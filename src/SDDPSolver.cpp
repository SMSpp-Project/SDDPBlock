/*--------------------------------------------------------------------------*/
/*-------------------------- File SDDPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPSolver class.
 *
 * \version 0.10
 *
 * \date 24 - 02 - 2021
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BendersBlock.h"
#include "FRealObjective.h"
#include "Objective.h"
#include "SDDPBlock.h"
#include "SDDPSolver.h"
#include "StochasticBlock.h"

#include <Eigen/Core>

#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"

#include "StOpt/sddp/backwardForwardSDDP.h"
#include "StOpt/sddp/LocalConstRegressionForSDDP.h"
#include "StOpt/sddp/LocalLinearRegressionForSDDP.h"

#include "StOpt/sddp/LocalConstRegressionForSDDPGeners.h"
#include "StOpt/sddp/LocalLinearRegressionForSDDPGeners.h"

#define BENDERSBFUNCTION_DEBUG

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( SDDPSolver );

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS of SDDPSolver -------------------------*/
/*--------------------------------------------------------------------------*/

int SDDPSolver::compute( bool changedvars ) {

 if( ! f_Block )
  return( kBlockLocked );

 // Possibly lock the SDDPBlock

 auto owned = f_Block->is_owned_by( f_id );        // check if already locked
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )  // if not try to lock
  return( kBlockLocked );                          // return error on failure

 process_outstanding_Modification();

 // ostream for StOpt output
 boost::iostreams::stream< boost::iostreams::null_sink >
  null_sink( ( boost::iostreams::null_sink() ) );
 std::ostream * output_stream = & null_sink;
 if( f_log && log_verbosity >= 2 )
  output_stream = f_log;

 const auto time_horizon = get_time_horizon();

 /* "dates" must be an array with size T + 1, where T is the time_horizon,
  * such that dates[ t ] contains the t-th time step (in our case it is simply
  * t) for each t in {0, ..., T-1}. The element in this array, with index T,
  * is associated with the cut to be used at the last time step. */

 Eigen::ArrayXd dates =
  Eigen::ArrayXd::LinSpaced( time_horizon + 1 , 0 , time_horizon );

 /* As input to the StOpt SDDP solver, it contains the maximum number of
  * iterations that the solver should perform. As output, it contains the
  * number of iterations performed by the StOpt SDDP solver. */
 number_iterations_performed = maximum_number_iterations;

 /* As input to the StOpt SDDP solver, it contains the desired accuracy that
  * the method should seek. As output, it contains the accuracy achieved by
  * the StOpt SDDP solver, which is given by
  *
  * | backwardValue - forwardValue | / forwardValue
  *
  * where backwardValue is the value of the last backward pass and
  * forwardValueForConv is the value obtained during the forward pass when
  * checking for convergence. */

 auto accuracy_achieved_stopt = accuracy;

 /*****************/
 /* INITIAL STATE */
 /*****************/

 if( ! initial_state.size() ) // TODO add a parameter to set initial_state
  initial_state = sddp_optimizer->oneAdmissibleState( 0 );

 /***************************/
 /* CUTS FOR THE LAST STAGE */
 /***************************/

 /* The cuts to be used at the last time instant. If any cuts are provided,
  * they are added to the last stage right now and do not need to be passed to
  * StOpt. Notice that any cut that is currently at the last stage are kept
  * there. */

 PolyhedralFunction::MultiVector A;
 PolyhedralFunction::RealVector b;
 const auto number_of_states = initial_state.size();

 if( ! last_stage_cuts.empty() ) {

  if( last_stage_cuts.size() % ( number_of_states + 1 ) != 0 )
   throw( std::logic_error( "SDDPSolver::compute: Invalid size of the cuts "
                            "for the last stage." ) );

  const auto number_of_cuts = last_stage_cuts.size() / ( number_of_states + 1 );

  // Pack the cuts

  A.resize( number_of_cuts );
  b.resize( number_of_cuts );

  for( Index i = 0 ; i < number_of_cuts ; ++i ) {
   const auto begin = i * ( number_of_states + 1 );
   b[ i ] = last_stage_cuts[ begin + number_of_states ];
   A[ i ].resize( number_of_states );
   for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j )
    A[ i ][ j ] = last_stage_cuts[ begin + j ];
  }
 }
 else {
  /* No cut for the last stage has been provided. We check whether the
   * PolyhedralFunction at the last stage contains a finite bound or at least
   * one cut (row). */

  const auto polyhedral_function =
   static_cast< SDDPBlock * >( f_Block )->get_polyhedral_functions().back();

  if( ( ! polyhedral_function->is_bound_set() ) &&
      ( polyhedral_function->get_nrows() == 0 ) ) {
   if( f_log )
    *f_log << "Warning: SDDPSolver::compute: No cut for the last stage has been"
           << " provided and the PolyhedralFunction at the last stage has no"
           << " bound and no row (cut). By default, the all-zero cut will then "
           << " be used for the last stage." << std::endl;
   b.resize( 1 , 0 );
   A.resize( 1 );
   A.front().resize( number_of_states , 0 );
  }
 }

 if( ! A.empty() ) {
  // Add the cuts to the last stage
  static_cast< SDDPBlock * >( f_Block )->add_cuts
   ( std::move( A ) , std::move( b ) , get_time_horizon() - 1 );
 }

 /* The cuts for the problem at the last stage have just been added. StOpt
  * requires a StOpt::SDDPFinalCut as argument, representing the cuts for the
  * last stage. So, we create a dummy one, which will be ignored within
  * oneStepBackward() and oneStepForward() when the problem at the last stage
  * is being solved. */

 StOpt::SDDPFinalCut final_cut
  ( Eigen::ArrayXXd::Zero( initial_state.size() + 1 , 1 ) );

 /* In some situations, we can avoid adding the cuts every time in
  * oneStepBackward() and oneStepForward(). One of these situations, for
  * instance, happens when no mesh discretization is provided. In this case,
  * the set of cuts to be considered in oneStepBackward() at any given time is
  * always a superset of that provided at all previous iterations for that
  * same time. So, we can add only the new cuts. The following vector stores
  * the number of cuts currently present at each stage and help determine
  * which cuts were generated by StOpt. */

 number_initial_cuts.resize( get_time_horizon() );
 for( Index stage = 0 ; stage < get_time_horizon() ; ++stage )
  number_initial_cuts[ stage ] =
   static_cast< SDDPBlock * >( f_Block )->get_number_cuts( stage );

 /***********************/
 /* MESH DISCRETIZATION */
 /***********************/

 // Meshes for regression
 Eigen::ArrayXi mesh_discretization_array;
 if( mesh_provided() ) {
  auto simulator = std::static_pointer_cast< ScenarioSimulator >
   ( sddp_optimizer->getSimulatorBackward() );
  const auto particle_length = simulator->get_particle_length();
  if( mesh_discretization.size() != particle_length )
   throw( std::logic_error
         ( "SDDPSolver::compute: simulation particle has dimension " +
           std::to_string( particle_length ) + " but given mesh discretization "
           "has dimension " + std::to_string( mesh_discretization.size() ) ) );
  mesh_discretization_array.resize( mesh_discretization.size() );
  for( Index i = 0 ; i < mesh_discretization.size() ; ++i )
   mesh_discretization_array( i ) = mesh_discretization[ i ];
 }

 /*********************/
 /* CONTROL VARIABLES */
 /*********************/

 current_iteration = 0;
 previous_pass_was_backward = false;

 /***********************/
 /* SOLVING THE PROBLEM */
 /***********************/

 // Invoke the StOpt SDDP solver
 auto backward_forward_values =
  StOpt::backwardForwardSDDP<StOpt::LocalLinearRegressionForSDDP>
  ( sddp_optimizer , number_simulations_for_convergence , initial_state ,
    final_cut , dates , mesh_discretization_array , regressors_filename ,
    cuts_filename , visited_states_filename , number_iterations_performed ,
    accuracy_achieved_stopt , convergence_frequency , *output_stream ,
    print_cpu_time );

 // Retrieve the backward and forward values

 backward_value = backward_forward_values.first;
 forward_value = backward_forward_values.second;

 // Log

 if( f_log && log_verbosity > 0 ) {
  *f_log << "Backward value: " << std::setprecision( 20 )
         << backward_value << std::endl;
  *f_log << "Forward value:  " << std::setprecision( 20 )
         << forward_value << std::endl;
 }

 // Unlock the SDDPBlock

 if( ! owned )              // if the Block was actually locked
  f_Block->unlock( f_id );  // unlock it

 // Compute the accuracy achieved

 if( forward_value != 0.0 )
  accuracy_achieved = std::abs( ( backward_value - forward_value ) /
                                forward_value );
 else
  accuracy_achieved = std::abs( backward_value );

 // Determine the status of SDDPSolver

 if( accuracy_achieved_stopt == 0.0 && accuracy_achieved != 0.0 )
  status = kCurveCross;
 else if( accuracy_achieved_stopt <= accuracy )
  status = kOK;
 else if( number_iterations_performed == maximum_number_iterations )
  status = kStopIter;
 else
  status = kError; // TODO

 // Possibly output the future cost functions

 if( output_frequency > 0 )
  output_future_cost_functions( f_output_filename );

 return status;
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::get_lb( void ) {
 if( ! f_Block )
  return - Inf< double >();

 const auto minimization =
  ( f_Block->get_objective_sense() == Objective::eMin );

 if( status == kOK || status == kLowPrecision ) {
  if( minimization )
   return backward_value;
  return forward_value;
 }

 if( status == kInfeasible ) {
  if( minimization )
   return Inf< double >();
  return - Inf< double >();
 }

 return - Inf< double >();
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::get_ub( void ) {
 if( ! f_Block )
  return Inf< double >();

 const auto minimization =
  ( f_Block->get_objective_sense() == Objective::eMin );

 if( status == kOK || status == kLowPrecision ) {
  if( minimization )
   return forward_value;
  return backward_value;
 }

 if( status == kInfeasible ) {
  if( minimization )
   return Inf< double >();
  return - Inf< double >();
 }

 return Inf< double >();
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::process_outstanding_Modification() {
 v_mod.clear();
}

/*--------------------------------------------------------------------------*/

Eigen::ArrayXd SDDPSolver::SDDPOptimizer::oneStepBackward
( const StOpt::SDDPCutOptBase & sddp_cut ,
  const std::tuple< std::shared_ptr< Eigen::ArrayXd > , int , int > & state ,
  const Eigen::ArrayXd & particle , const int & simulation_id ) const {

 const auto current_stage = date_next;

 /* The last argument is the simulation id indicating in which scenario the
  * resolution will be done. */

 Index scenario_index = 0;
 if( current_stage > 0 )
  scenario_index = simulator_backward->get_scenario_index( simulation_id );

 // Log

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  auto log = sddp_solver->f_log;
  *log << "***** SDDPSolver::SDDPOptimizer::oneStepBackward *****" << std::endl;
  *log << "  Stage:          " << current_stage << std::endl;
  *log << "  Scenario index: " << scenario_index << std::endl;
  if( current_stage > 0 ) {
   *log << "  Simulation id:  " << simulation_id << std::endl;
   if( sddp_solver->log_verbosity >= 4 ) {
    *log << "  Particle:       (";
    for( decltype( particle.size() ) i = 0 ; i < particle.size() ; ++i ) {
     if( i > 0 ) *log << ", ";
     *log << particle( i );
    }
    *log << ")" << std::endl;
   }
  }

  if( sddp_solver->log_verbosity >= 10 ) {
   *log << "  State:          (";
   const auto & state_variables = * std::get<0>( state );
   for( decltype( state_variables.size() ) i = 0 ;
        i < state_variables.size() ; ++i ) {
    if( i > 0 ) *log << ", ";
    *log << state_variables( i );
   }
   *log << ")" << std::endl;
  }
 }

 /***************/
 /* ADDING CUTS */
 /***************/

 if( current_stage < sddp_solver->get_time_horizon() - 1 ) {
  /* The cuts for the last stage are added only once in the beginning of
   * compute(). */

  const auto cuts = sddp_cut.getCutsAssociatedToTheParticle
   ( std::get<1>( state ) );
  sddp_solver->add_cuts( cuts , current_stage );
 }

 /*******************/
 /* STATE VARIABLES */
 /*******************/

 sddp_solver->set_state( * std::get<0>( state ).get() , current_stage );

 /***************/
 /* RANDOM DATA */
 /***************/

 sddp_solver->set_scenario( scenario_index , current_stage );

 /**************************/
 /* SOLVING THE SUBPROBLEM */
 /**************************/

 auto objective_value = sddp_solver->solve( current_stage );

 /**********************************/
 /* CONSTRUCTING THE LINEARIZATION */
 /**********************************/

 /* The oneStepBackard function returns a one-dimensional array whose size is
  * the number of state variables plus one and that contains a linearization
  * of the BendersBFunction. The first component contains the value of the
  * BendersBFunction and the remaining components contain the coefficients of
  * the linearization of the BendersBFunction. For i in {1, ...,
  * number_state_variables}, linearization( i ) contains the coefficient of
  * the linearization of the BendersBFunction associated with the i-th state
  * variable. */

 const auto number_state_variables = std::get<0>( state )->size();
 Eigen::ArrayXd linearization( number_state_variables + 1 );

 auto benders_function = sddp_solver->get_benders_function( current_stage );

 if( benders_function->has_linearization( true ) ) {
  linearization( 0 ) = objective_value;
  benders_function->get_linearization_coefficients( linearization.data() + 1 );

  if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
   *( sddp_solver->f_log ) << "  Objective:      " << objective_value
                           << std::endl;
  }

  // Debugging the BendersBFunction

#ifdef BENDERSBFUNCTION_DEBUG
  {
   const auto alpha = benders_function->get_linearization_constant();
   if( sddp_solver->f_log && sddp_solver->log_verbosity >= 20 ) {
    *( sddp_solver->f_log ) << "  Linearization: " << std::endl;
    *( sddp_solver->f_log ) << "    alpha:        " << alpha << std::endl;
    if( sddp_solver->log_verbosity >= 30 ) {
     *( sddp_solver->f_log ) << "    coefficients: (";
     for( decltype( linearization.size() ) i = 1 ;
          i < linearization.size() ; ++i ) {
      if( i > 1 ) *( sddp_solver->f_log ) << ", ";
      *( sddp_solver->f_log ) << linearization( i );
     }
     *( sddp_solver->f_log ) << ")" << std::endl;
    }
   }

   double gy = 0;
   const auto & state_variables = * std::get<0>( state ).get();
   for( decltype( state_variables.size() ) j = 0 ;
        j < state_variables.size() ; ++j )
    gy += linearization( j + 1 ) * state_variables( j );
   const double epsilon = 1.0e-4;
   const auto scale =
    std::max( 1.0 , std::min( abs( objective_value ) , abs( alpha + gy ) ) );
   const auto diff = std::abs( objective_value - ( alpha + gy ) );
   if( diff > epsilon * scale ) {
    std::cerr << "SDDPOptimizer::oneStepBackward: linearization precision "
              << "was not achieved:" << std::endl;
    std::cerr << "  precision required: " << std::setprecision( 20 )
              << epsilon << std::endl;
    std::cerr << "  precision achieved: " << std::setprecision( 20 )
              << ( diff / scale ) << std::endl;
    std::cerr << "  objective: " << std::setprecision( 20 )
              << objective_value << std::endl;
    std::cerr << "  alpha:     " << std::setprecision( 20 )
              << alpha << std::endl;
    std::cerr << "  g'y:       " << std::setprecision( 20 ) << gy << std::endl;
   }
  }
#endif

 }
 else if( benders_function->has_linearization( false ) ) {
  const auto alpha = benders_function->get_linearization_constant();
  linearization( 0 ) = alpha;
  benders_function->get_linearization_coefficients( linearization.data() + 1 );
 }
 else {
  // No linearization is available
  throw( std::logic_error( "SDDPOptimizer::oneStepBackward: no "
                           "linearization is available." ) );
 }

 sddp_solver->previous_pass_was_backward = true;

 return linearization;
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::SDDPOptimizer::oneStepForward
( const Eigen::ArrayXd & particle , Eigen::ArrayXd & state ,
  Eigen::ArrayXd & state_to_store , const StOpt::SDDPCutOptBase & sddp_cut ,
  const int & simulation_id ) const {

 const auto current_stage = date;

 // Index of the scenario to be considered

 Index scenario_index = 0;
 if( current_stage > 0 )
  scenario_index = simulator_forward->get_scenario_index( simulation_id );

 // Log

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  auto log = sddp_solver->f_log;
  *log << "***** SDDPSolver::SDDPOptimizer::oneStepForward *****"
            << std::endl;
  *log << "  Stage:          " << current_stage << std::endl;
  *log << "  Scenario index: " << scenario_index << std::endl;
  if( current_stage > 0 ) {
   *log << "  Simulation id:  " << simulation_id << std::endl;
   if( sddp_solver->log_verbosity >= 4 ) {
    *log << "  Particle:       (";
    for( decltype( particle.size() ) i = 0 ; i < particle.size() ; ++i ) {
     if( i > 0 ) *log << ", ";
     *log << particle( i );
    }
    *log << ")" << std::endl;
   }
  }

  if( sddp_solver->log_verbosity >= 10 ) {
   *log << "  State:          (";
   for( decltype( state.size() ) i = 0 ; i < state.size() ; ++i ) {
    if( i > 0 ) *log << ", ";
    *log << state( i );
   }
   *log << ")" << std::endl;
  }
 }

 /***************/
 /* ADDING CUTS */
 /***************/

 if( sddp_solver->mesh_provided() && ( current_stage > 0 )  &&
     ( current_stage < sddp_solver->get_time_horizon() - 1 ) ) {
  /* Cuts are added if and only if meshes were provided and the current stage
   * is not the first or the last one. There is a single mesh associated with
   * the first stage and the cuts have already been added during the backward
   * pass. The cuts for the last stage are added once in the beginning of
   * compute(). The array "particle" contains the random quantities on which
   * the regression over the expectation of the value function will be
   * based. */
  const auto cuts = sddp_cut.getCutsAssociatedToAParticle( particle );
  sddp_solver->add_cuts( cuts , current_stage );
 }

 /***************/
 /* RANDOM DATA */
 /***************/

 sddp_solver->set_scenario( scenario_index , current_stage );

 /*********************************/
 /* VARIABLES FROM PREVIOUS STAGE */
 /*********************************/

 sddp_solver->set_state( state , current_stage );

 /**************************/
 /* SOLVING THE SUBPROBLEM */
 /**************************/

 auto objective_value = sddp_solver->solve( current_stage );

 /* The objective_value takes into account the value of the future cost
  * function. For all stages other than the last one, we subtract the value of
  * the future cost function from objective_value. The problem at the last
  * stage, however, has a fixed future cost and it is kept as it is considered
  * part of the cost of that stage. */

 if( current_stage < sddp_solver->get_time_horizon() - 1 )
  objective_value -= static_cast< SDDPBlock * >( sddp_solver->f_Block )->
   get_future_cost( current_stage );

 /**************************/
 /* RETRIVING THE SOLUTION */
 /**************************/

 // Retrieve the solution x_t of the Block associated with the current stage.

 auto solution = sddp_solver->get_solution( current_stage );

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  *( sddp_solver->f_log ) << "  Objective:      " << objective_value << std::endl;
  if( sddp_solver->log_verbosity >= 10 ) {
   *( sddp_solver->f_log ) << "  Solution:       (";
   for( decltype( solution.size() ) i = 0 ; i < solution.size() ; ++i ) {
    if( i > 0 ) *( sddp_solver->f_log ) << ", ";
    *( sddp_solver->f_log ) << solution( i );
   }
   *( sddp_solver->f_log ) << ")" << std::endl;
  }
 }

 // Store in state the current state, i.e, ( x_t, w_t^{dep} ).

 state.resize( solution.size() );
 state << solution;

 // Store in state_to_store the vector ( x_t, w_{t-1}^{dep} ).

 state_to_store.resize( solution.size() );
 state_to_store << solution;

 /****************************/
 /* UPDATE CONTROL VARIABLES */
 /****************************/

 if( sddp_solver->previous_pass_was_backward ) {
  if( sddp_solver->output_frequency > 0 &&
      ( sddp_solver->current_iteration % sddp_solver->output_frequency == 0 ) )
   sddp_solver->output_future_cost_functions( sddp_solver->f_output_filename );
  sddp_solver->current_iteration++;
  sddp_solver->previous_pass_was_backward = false;
 }

 /*************************/
 /* RETURN SOLUTION VALUE */
 /*************************/

 return objective_value;
}

/*--------------------------------------------------------------------------*/

SDDPBlock::Index SDDPSolver::get_time_horizon( void ) const {
 return static_cast< SDDPBlock * >( f_Block )->get_time_horizon();
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::add_cuts( const Eigen::ArrayXXd & cuts ,
                           SDDPBlock::Index stage ) const {

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::add_cuts: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 // Number of cuts that will be added. In principle, all given cuts are added.

 auto number_cuts_to_be_added = cuts.cols();

 /* If a mesh has been provided, all cuts are replaced by the given ones, but
  * the initial cuts (those there were present at the time StOpt was called)
  * are kept. */
 auto number_cuts_to_keep = number_initial_cuts[ stage ];

 if( ! mesh_provided() ) {

  // Total number of cuts previously added by StOpt
  auto number_cuts_previously_added =
   static_cast< SDDPBlock * >( f_Block )->get_number_cuts( stage ) -
   number_initial_cuts[ stage ];

  assert( cuts.cols() >= number_cuts_previously_added );

  /* StOpt provides all cuts that were ever generated. Since no mesh has been
   * provided, we can consider only the most recent cuts (since the previous
   * ones have already been added before). We are assuming that the new cuts
   * provided by StOpt appear in the last columns of the matrix "cuts". */

  if( cuts.cols() == number_cuts_previously_added )
   return; // all cuts are already there

  // Add only the most recent cuts
  number_cuts_to_be_added = cuts.cols() - number_cuts_previously_added;

  // and keep the current ones
  number_cuts_to_keep = Inf< Index >();
 }

 // Store the given cuts in A and b

 PolyhedralFunction::MultiVector A( number_cuts_to_be_added );
 PolyhedralFunction::RealVector b( number_cuts_to_be_added );

 /* Each column in "cuts" contains a cut. The first element is the constant
  * term of the cut, and all the other elements are the coefficients. */

 for( Index k = 0 ; k < number_cuts_to_be_added ; ++k ) {
  const auto i = number_cuts_to_be_added - k - 1;
  const auto col = cuts.cols() - k - 1;
  A[ i ].resize( cuts.rows() - 1 );
  b[ i ] = cuts( 0 , col );
  for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j ) {
   A[ i ][ j ] = cuts( j + 1 , col );
  }
 }

 // Log

 if( f_log && log_verbosity >= 30 ) {
  *f_log << "  Adding the following cuts:" << std::endl;
  for( decltype( b )::size_type i = 0 ; i < b.size() ; ++i ) {
   *f_log << "    (" << b[ i ];
   for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j )
    *f_log << ", " << A[ i ][ j ];
   *f_log << ")" << std::endl;
  }
 }

 // Finally add the cuts

 static_cast< SDDPBlock * >( f_Block )->add_cuts
  ( std::move( A ) , std::move( b ) , stage , number_cuts_to_keep );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::set_state( const Eigen::ArrayXd & state ,
                            SDDPBlock::Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::set_state: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_state( state , stage );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::set_scenario( SDDPBlock::Index scenario_id ,
                               SDDPBlock::Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::set_scenario: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_scenario( scenario_id , stage );
}

/*--------------------------------------------------------------------------*/

BendersBFunction *
SDDPSolver::get_benders_function( SDDPBlock::Index stage ) const {
 auto benders_block = static_cast< BendersBlock * >
  ( static_cast< StochasticBlock * >( f_Block->get_nested_Blocks()[ stage ] )->
    get_nested_Blocks().front() );

 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );

 return static_cast< BendersBFunction * >( objective->get_function() );
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::solve( SDDPBlock::Index stage ) {

 /* Solving the subproblem consists in evaluating the Objective of the
  * BendersBFunction associated with the problem of the given stage. */

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::solve: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 auto benders_block = static_cast< BendersBlock * >
  ( static_cast< StochasticBlock * >( f_Block->get_nested_Blocks()[ stage ] )->
    get_nested_Blocks().front() );

 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );

 auto benders_function = static_cast< BendersBFunction * >
  ( objective->get_function() );

 auto status = benders_function->compute();

 if( status != kOK && status != kLowPrecision ) {
  // No feasible solution has been found

  std::string message;
  if( status == kUnbounded )
   message = " The sub-problem is unbounded.";
  else if( status == kInfeasible )
   message = " The sub-problem is infeasible.";
  else if( status == kError )
   message = " An error occurred while solving the sub-problem.";
  else if( status == kStopTime )
   message = " A time limit has been reached while solving the sub-problem.";
  else if( status == kStopIter )
   message = " A maximum number of iterations has been reached while solving "
    "the sub-problem.";

  throw( std::logic_error( "SDDPSolver::solve: the sub-problem at stage " +
                           std::to_string( stage ) + " was not solved." +
                           message ) );
 }

 auto solver = benders_function->get_solver();
 if( ! solver->has_var_solution() ) {
  throw( std::logic_error( "SDDPSolver::solve: the sub-problem at stage " +
                           std::to_string( stage ) + " has no solution." ) );
 }

 solver->get_var_solution(); // TODO Use Configuration to request only the
                             // active Variables of the PolyhedralFunction

 return benders_function->get_value();
}

/*--------------------------------------------------------------------------*/

template<class T>
T SDDPSolver::get_solution( SDDPBlock::Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::get_solution: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 const auto polyhedral_function =
  static_cast< SDDPBlock * >( f_Block )->get_polyhedral_functions()[ stage ];

 T solution( polyhedral_function->get_num_active_var() );

 auto data = solution.data();
 for( const auto & variable : * polyhedral_function ) {
  *data = static_cast< const ColVariable & >( variable ).get_value();
  data++;
 }
 return solution;
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::output_future_cost_functions( const std::string & filename )
 const {

 const auto & functions =
  static_cast< SDDPBlock *>( f_Block )->get_polyhedral_functions();
 if( functions.empty() )
  return;

 std::ofstream output( filename , std::ios::out );

 const char separator_character = ',';
 const auto num_var = functions.front()->get_num_active_var();

 output << "Timestep";
 for( Index i = 0 ; i < num_var ; ++i ) {
  output << separator_character << "a_" << std::to_string( i );
 }
 output << separator_character << "b" << std::endl;

 for( Index stage = 0 ; stage < get_time_horizon() ; ++stage ) {

  const auto & b = functions[ stage ]->get_b();
  const auto & A = functions[ stage ]->get_A();

  assert( b.size() == A.size() );

  for( Index i = 0 ; i < b.size() ; ++i ) {
   output << stage;
   for( Index j = 0 ; j < A[ i ].size() ; ++j )
    output << separator_character << std::setprecision( 20 ) << A[ i ][ j ];
   output << separator_character << std::setprecision( 20 ) << b[ i ]
          << std::endl;
  }
 }

 output.close();

}

/*--------------------------------------------------------------------------*/

StochasticBlock * SDDPSolver::SDDPOptimizer::get_block
( const double & stage ) const {

  /* Make sure the given stage is integer and belongs to the interval
   * [0, T-1], where T is the time horizon. */

  double integral_part;
  assert( std::modf( stage , &integral_part ) == 0.0 );

  assert( 0 <= stage && stage < sddp_solver->get_time_horizon() );

  // The number of sub-Blocks must be at least the time horizon.

  assert( sddp_solver->f_Block->get_nested_Blocks().size() >=
          sddp_solver->get_time_horizon() );

  // Return the Block associated with the given stage.

  return static_cast<StochasticBlock *>
    ( sddp_solver->f_Block->get_nested_Blocks()[ stage ] );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::SDDPOptimizer::updateDates
( const double & date , const double & date_next ) {

  // We assume that the given arguments correspond to consecutive stages
  assert( date_next == date + 1.0 );

  // We assume that -1 <= date < T
  assert( - 1.0 <= date && date < sddp_solver->get_time_horizon() );

  this->date = date;
  this->date_next = date_next;
}

/*--------------------------------------------------------------------------*/

Eigen::ArrayXd
SDDPSolver::SDDPOptimizer::oneAdmissibleState( const double & stage ) {

 const auto sddp_block = static_cast< SDDPBlock * >( sddp_solver->f_Block );

 auto state_size = sddp_block->get_admissible_state_size( stage );
 Eigen::ArrayXd state( state_size );

 auto state_iterator = sddp_block->get_admissible_state( stage );

 auto data = state.data();
 for( decltype( state_size ) i = 0 ; i < state_size ;
      ++i , ++data , ++state_iterator )
  *data = *state_iterator;

 return state;
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File SDDPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

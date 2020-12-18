/*--------------------------------------------------------------------------*/
/*-------------------------- File SDDPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPSolver class.
 *
 * \version 0.10
 *
 * \date 18 - 12 - 2020
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

 auto owned = f_Block->is_owned_by( f_id );        // check if already locked
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )  // if not try to lock
  return( kBlockLocked );                          // return error on failure

 process_outstanding_Modification();

 // ostringstream in which the text output of StOpt will be stored
 std::ostringstream output_stream;
 /*
 boost::iostreams::stream< boost::iostreams::null_sink >
  null_sink( ( boost::iostreams::null_sink() ) );
 std::ostream * output_stream = & null_sink;
 if( f_log && log_verbosity >= 5 )
  output_stream = f_log;
 */

 auto time_horizon = get_time_horizon();

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

 accuracy_achieved = accuracy;

 if( ! initial_state.size() ) // TODO add a parameter to set initial_state
  initial_state = sddp_optimizer->oneAdmissibleState( 0 );

 /*
 if( ! number_meshes.size() ) {
  // TODO add to SDDPSolver parameters
  auto simulator = std::static_pointer_cast< ScenarioSimulator >
   ( sddp_optimizer->getSimulatorBackward() );
  const auto particle_length = simulator->get_particle_length();
  number_meshes = Eigen::ArrayXi::Constant
   ( particle_length , std::min( 10 , simulator->getNbSimul() ) );
 }
 */

 // The cuts used at the last time step. TODO add a parameter to set final_cut
 StOpt::SDDPFinalCut final_cut
  ( Eigen::ArrayXXd::Zero( initial_state.size() + 1 , 1 ) );

 cut_controller.reset();

 // Invoke the StOpt SDDP solver
 auto backward_forward_values =
  StOpt::backwardForwardSDDP<StOpt::LocalLinearRegressionForSDDP>
  ( sddp_optimizer , number_simulations_for_convergence , initial_state ,
    final_cut , dates ,
    number_meshes ,
    regressors_filename , cuts_filename , visited_states_filename ,
    number_iterations_performed , accuracy_achieved , convergence_frequency ,
    //*output_stream , print_cpu_time );
    output_stream , print_cpu_time );

 backward_value = backward_forward_values.first;
 forward_value = backward_forward_values.second;

 if( f_log && log_verbosity > 0 ) {
  *f_log << "Backward value: " << backward_value << std::endl;
  *f_log << "Forward value:  " << forward_value << std::endl;
 }

 if( ! owned )              // if the Block was actually locked
  f_Block->unlock( f_id );  // unlock it

 if( accuracy_achieved <= accuracy )
  return( kOK );
 else if( number_iterations_performed == maximum_number_iterations )
  return( kStopIter );
 else
  return( kError ); // TODO
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::process_outstanding_Modification() {
 v_mod.clear();
}

/*--------------------------------------------------------------------------*/

Eigen::ArrayXd SDDPSolver::SDDPOptimizer::oneStepBackward
( const StOpt::SDDPCutOptBase &p_linCut,
  const std::tuple< std::shared_ptr<Eigen::ArrayXd>, int, int > & state,
  const Eigen::ArrayXd & particle, const int & simulation_id ) const {

 auto current_stage = date_next;

 if( sddp_solver->f_log && sddp_solver->log_verbosity ) {
  auto log = sddp_solver->f_log;
  *log << "***** SDDPSolver::SDDPOptimizer::oneStepBackward *****"
            << std::endl;
  *log << "  stage:         " << current_stage << std::endl;
  *log << "  simulation_id: " << simulation_id << std::endl;
  *log << "  particle:      (";
  for( decltype( particle.size() ) i = 0 ; i < particle.size() ; ++i ) {
   if( i > 0 ) *log << ", ";
   *log << particle( i );
  }
  *log << ")" << std::endl;

  if( sddp_solver->f_log && sddp_solver->log_verbosity >= 10 )
   *log << "  state:         " << *std::get<0>( state ) << std::endl;
 }

 /* "particle" contains the random quantities in which the regression over the
  * expectation of the value function will be based */

 auto cuts = p_linCut.getCutsAssociatedToTheParticle( std::get<1>( state ) );

 /* For each scenario j used in the forward pass (for j+1 in
  * {1, ..., G, ..., (n+1)G}), cuts(:, j) is a cut such that
  *
  * cuts(0, j) = alpha_{t+1}^j
  *
  * cuts(i, j) = beta_{i-1, t+1}^j, for i in {1, ..., nbstate}
  *
  * Add these cuts to the Block associated with time p_dateNext.
  */

 if( current_stage == sddp_solver->get_time_horizon() - 1 )
  sddp_solver->add_cuts( cuts , current_stage , true );
 else {
  if( cuts.cols() > 1 )
   assert( cuts.cols() >= simulator_forward->getNbSimul() );
  auto first_cut = cuts.cols() - 1;
  if( cuts.cols() >= simulator_forward->getNbSimul() )
   first_cut = cuts.cols() - simulator_forward->getNbSimul();
  sddp_solver->add_cuts( cuts , current_stage , true ,
                         Block::Range( first_cut , cuts.cols() ) );
  //   ( cuts( Eigen::all , Eigen::lastN( simulator_forward->getNbSimul() ) ) ,
  //     current_stage );
 }

 /*******************/
 /* STATE VARIABLES */
 /*******************/

 sddp_solver->set_state( * std::get<0>( state ).get() , current_stage );

 /***************/
 /* RANDOM DATA */
 /***************/

 /* The last argument is the simulation id indicating in which scenario the
  * resolution will be done. */

 auto scenario_index = simulator_backward->get_scenario_index( simulation_id );
 sddp_solver->set_scenario( scenario_index , current_stage );

 // Solve the problem

 auto objective_value = sddp_solver->solve( current_stage );

 /* This function returns a one-dimensional array whose size is the number of
  * state variables plus one and that contains a linearization of the
  * BendersBFunction. The first component contains the value of the
  * BendersBFunction and the remaining components contain the coefficients of
  * the linearization of the bendersBFunction. For i in {1, ...,
  * number_state_variables}, linearization( i ) contains the coefficient of
  * the linearization of the BendersBFunction associated with the i-th state
  * variable. */

 auto number_state_variables = std::get<0>( state )->size();
 Eigen::ArrayXd linearization( number_state_variables + 1 );

 auto benders_function = sddp_solver->get_benders_function( current_stage );

 if( benders_function->has_linearization( true ) ) {
  linearization( 0 ) = objective_value;
  benders_function->get_linearization_coefficients( linearization.data() + 1 );

#ifdef BENDERSBFUNCTION_DEBUG
  {
   const auto alpha = benders_function->get_linearization_constant();
   if( sddp_solver->f_log && sddp_solver->log_verbosity ) {
    *( sddp_solver->f_log ) << "  alpha:           " << alpha << std::endl;
    *( sddp_solver->f_log ) << "  objective value: " << objective_value
                            << std::endl;
    if( sddp_solver->log_verbosity >= 0 ) {
     *( sddp_solver->f_log ) << "  linearization coefficients:\n    (";
     for( decltype( linearization.size() ) i = 1 ;
          i < linearization.size() ; ++i ) {
      if( i > 1 ) *( sddp_solver->f_log ) << ", ";
      *( sddp_solver->f_log ) << linearization( i );
     }
     *( sddp_solver->f_log ) << ")" << std::endl;
    }
   }

   double gy = 0;
   auto state_variables = std::get<0>( state ).get();
   for( decltype( ( *state_variables ).size() ) j = 0 ;
        j < ( *state_variables ).size() ; ++j )
    gy += linearization( j + 1 ) * ( *state_variables )( j );
   const double epsilon = 1.0e-8;
   const auto max_diff = std::max( epsilon , epsilon *
                                   std::min( abs( objective_value ),
                                             abs( alpha + gy ) ) );
   if( std::abs( objective_value - ( alpha + gy ) ) > max_diff ) {
    std::cerr << "Wrong linearization in SDDPSolver:" << std::endl;
    std::cerr << "obj   = " << std::setprecision( 20 )
              << objective_value << std::endl;
    std::cerr << "alpha = " << std::setprecision( 20 ) << alpha << std::endl;
    std::cerr << "gy    = " << std::setprecision( 20 ) << gy << std::endl;
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

 sddp_solver->cut_controller.backward_pass( current_stage );

 return linearization;
}

/*--------------------------------------------------------------------------*/

SDDPBlock::Index SDDPSolver::get_time_horizon( void ) const {
 return static_cast< SDDPBlock * >( f_Block )->get_time_horizon();
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::add_cuts( const Eigen::ArrayXXd & cuts ,
                           SDDPBlock::Index stage , bool backward ,
                           Block::Range range ) const {

 if( ! cut_controller.add_cuts( stage , get_time_horizon() ) )
  return;

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::add_cuts: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 /* The first element is the function value, and all the other elements are
  * the coefficients. */

 PolyhedralFunction::MultiVector A;
 PolyhedralFunction::RealVector b;

 range.second = std::min< decltype( cuts.cols() ) >( range.second ,
                                                     cuts.cols() );

 if( range.second <= range.first )
  return;

 auto num_cuts = range.second - range.first;

 A.resize( num_cuts );
 b.resize( num_cuts );

 for( decltype( num_cuts ) cut = range.first ; cut < range.second ; ++cut ) {
  auto i = cut - range.first;
  A[ i ].resize( cuts.rows() - 1 );
  b[ i ] = cuts( 0 , cut );
  for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j ) {
   A[ i ][ j ] = cuts( j + 1 , cut );
  }
 }

 if( f_log && log_verbosity >= 5 ) {
  *f_log << "  Adding the following cuts:" << std::endl;
  for( decltype( b )::size_type i = 0 ; i < b.size() ; ++i ) {
   *f_log << "    (" << b[ i ];
   for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j )
    *f_log << ", " << A[ i ][ j ];
   *f_log << ")" << std::endl;
  }
 }

 auto replace_last_cuts = cut_controller.remove_cuts
  ( stage , get_time_horizon() , backward );

 static_cast< SDDPBlock * >( f_Block )->add_cuts
  ( std::move( A ) , std::move( b ) , stage , replace_last_cuts );
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
  throw( std::logic_error( "SDDPSolver::solve: the sub-problem at stage " +
                           std::to_string( stage ) + " was not solved." ) );
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

double SDDPSolver::SDDPOptimizer::oneStepForward
( const Eigen::ArrayXd & particle , Eigen::ArrayXd &state ,
  Eigen::ArrayXd & state_to_store ,
  const StOpt::SDDPCutOptBase & sddp_cut ,
  const int & simulation_id ) const {

 auto current_stage = date;

 if( sddp_solver->f_log && sddp_solver->log_verbosity ) {
  auto log = sddp_solver->f_log;
  *log << "***** SDDPSolver::SDDPOptimizer::oneStepForward *****"
            << std::endl;
  *log << "  stage:         " << current_stage << std::endl;
  *log << "  simulation_id: " << simulation_id << std::endl;
  *log << "  particle:      (";
  for( decltype( particle.size() ) i = 0 ; i < particle.size() ; ++i ) {
   if( i > 0 ) *log << ", ";
   *log << particle( i );
  }
  *log << ")" << std::endl;
  if( sddp_solver->f_log && sddp_solver->log_verbosity >= 10 )
   *log << "  state:         " << state << std::endl;
 }

 /********/
 /* CUTS */
 /********/

 // Retrieve the cuts.
 auto cuts = sddp_cut.getCutsAssociatedToAParticle( particle );

 /* StOpt provides all cuts that were ever generated. We consider only
  * the last ones when updating the cuts. */

 // Update the cuts in the Block associated with the current stage.

 if( current_stage == sddp_solver->get_time_horizon() - 1 )
  sddp_solver->add_cuts( cuts , current_stage , false );
 else {
  if( cuts.cols() > 1 )
   assert( cuts.cols() >= simulator_forward->getNbSimul() );
  auto first_cut = cuts.cols() - 1;
  if( cuts.cols() >= simulator_forward->getNbSimul() )
   first_cut = cuts.cols() - simulator_forward->getNbSimul();
  sddp_solver->add_cuts( cuts , current_stage , false ,
                         Block::Range( first_cut , cuts.cols() ) );
 }

 /***************/
 /* RANDOM DATA */
 /***************/

 auto scenario_index = simulator_forward->get_scenario_index( simulation_id );
 sddp_solver->set_scenario( scenario_index , current_stage );

 /*********************************/
 /* VARIABLES FROM PREVIOUS STAGE */
 /*********************************/

 sddp_solver->set_state( state , current_stage );

 /**************************/
 /* SOLVING THE SUBPROBLEM */
 /**************************/

 auto objective_value = sddp_solver->solve( current_stage ) -
  static_cast< SDDPBlock * >( sddp_solver->f_Block )->
  get_future_cost( current_stage );

 /**************************/
 /* RETRIVING THE SOLUTION */
 /**************************/

 // Retrieve the solution x_t of the Block associated with the current stage.

 auto solution = sddp_solver->get_solution( current_stage );

 if( sddp_solver->f_log && sddp_solver->log_verbosity ) {
  *( sddp_solver->f_log ) << "  objective: " << objective_value << std::endl;
  if( sddp_solver->f_log && sddp_solver->log_verbosity >= 10 )
   *( sddp_solver->f_log ) << "  solution:  " << solution << std::endl;
 }

 // Store in state the current state, i.e, ( x_t, w_t^{dep} ).

 state.resize( solution.size() );
 state << solution;

 // Store in state_to_store the vector ( x_t, w_{t-1}^{dep} ).

 state_to_store.resize( solution.size() );
 state_to_store << solution;

 /*************************/
 /* RETURN SOLUTION VALUE */
 /*************************/

 return objective_value;
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

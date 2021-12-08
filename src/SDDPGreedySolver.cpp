/*--------------------------------------------------------------------------*/
/*----------------------- File SDDPGreedySolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPGreedySolver class.
 *
 * \version 0.10
 *
 * \date 04 - 12 - 2021
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
#include "BlockSolverConfig.h"
#include "CDASolver.h"
#include "FRealObjective.h"
#include "SDDPBlock.h"
#include "SDDPGreedySolver.h"
#include "StochasticBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( SDDPGreedySolver );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of SDDPGreedySolver ----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_Block( Block * block ) {

 if( f_Block == block )
  return;

 if( f_Block ) {
  // TODO clean
  v_inner_block_configured.clear();
  v_inner_solver_configured.clear();
 }

 Solver::set_Block( block );

 if( ! f_Block )
  return;

 SDDPBlock * sddp_block;
 if( ! ( sddp_block = dynamic_cast< SDDPBlock * >( block ) ) )
  throw( std::invalid_argument( "SDDPGreedySolver::set_Block: given Block "
                                "is not an SDDPBlock." ) );

 v_inner_block_configured.assign( sddp_block->get_time_horizon() , false );
 v_inner_solver_configured.assign( sddp_block->get_time_horizon() , false );

}  // end( SDDPGreedySolver::set_Block )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int SDDPGreedySolver::compute( bool changedvars ) {

 if( ! f_Block )
  return kBlockLocked;

 // Possibly lock the SDDPBlock

 auto owned = f_Block->is_owned_by( f_id );        // check if already locked
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )  // if not try to lock
  return( kBlockLocked );                          // return error on failure

 process_outstanding_Modification();

 const auto time_horizon = get_time_horizon();
 status_compute = Solver::kLowPrecision;
 fault_stage = Inf<Index>();
 solution_value = 0.0;
 f_has_var_solution = false;

 // Initial state for the first stage problem
 if( time_horizon > 0 ) {
  if( ! v_initial_state.empty() )
   // Use the initial state given as parameter to the SDDPGreedySolver
   set_state( v_initial_state , 0 );
  else {
   const auto & state =
    static_cast< SDDPBlock * >( f_Block )->get_initial_state();
   if( ! state.empty() )
    // Use the initial state given by SDDPBlock
    set_state( state , 0 );
  }
 }

 // If required, load the random cuts
 if( ! f_random_cuts_filename.empty() )
  static_cast< SDDPBlock * >( f_Block )->
   deserialize_random_cuts( f_random_cuts_filename );

 for( Index stage = 0 ; stage < time_horizon ; ++stage ) {

  if( f_log && log_verbosity )
   *f_log << "Solving problem at stage " << stage << std::endl;

  if( stage > 0 ) {
   set_state( get_solution( stage - 1 ) , stage );
  }

  set_scenario( stage );

  if( callback ) callback( stage );

  configure_inner_block( stage );

  auto sub_status = solve( stage , true );

  if( sub_status == Solver::kInfeasible ) {
   fault_stage = stage;
   if( stage == 0 ) status_compute = kInfeasible;
   else status_compute = kSubproblemInfeasible;
   break;
  }
  else if( sub_status == Solver::kUnbounded ) {
   fault_stage = stage;
   status_compute = Solver::kUnbounded;
   break;
  }
  else if( sub_status >= Solver::kError ) {
   fault_stage = stage;
   status_compute = kError;
   break;
  }
  else if( ! get_sub_solver( stage )->has_var_solution() ) {
   fault_stage = stage;
   status_compute = kSolutionNotFound;
   break;
  }
  else if( sub_status == Solver::kStopTime || sub_status == Solver::kStopIter ) {
   solution_value += get_sub_solution_value( stage );
   if( fault_stage == Inf<Index>() ) {
    fault_stage = stage;
    status_compute = sub_status;
   }
  }
  else {
   solution_value += get_sub_solution_value( stage );

   if( ! f_subgradients_filename.empty() )
    store_subgradients( stage , scenario_id );
  }

  if( f_unregister_solver )
   unregister_solver_inner_block( stage );
 }

 // Unlock the SDDPBlock

 if( ! owned )              // if the Block was actually locked
  f_Block->unlock( f_id );  // unlock it

 f_has_var_solution =
  ( status_compute == Solver::kOK ) ||
  ( status_compute == Solver::kLowPrecision ) ||
  ( status_compute == Solver::kStopIter ) ||
  ( status_compute == Solver::kStopTime );

 // Output the subgradients if required.

 output_subgradients( f_subgradients_filename );

 return status_compute;
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING THE DATA ----------------------*/
/*--------------------------------------------------------------------------*/

SDDPGreedySolver::Index SDDPGreedySolver::get_time_horizon( void ) const {
 return static_cast< SDDPBlock * >( f_Block )->get_time_horizon();
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::get_var_solution( Configuration *solc ) {
 /* During the call to compute(), every Block associated with a stage in {0,
  * ..., T-2} has its solutions written in it. Therefore, we only need to
  * write the solution associated with the last stage. */

 auto solver = get_sub_solver( get_time_horizon() - 1 );

 if( ! solver )
  return; // The Solver must have been unregistered (but the Solution should
          // have already been written into the Block)

 if( ! solver->has_var_solution() )
  throw( std::logic_error( "SDDPGreedySolver::get_var_solution: subproblem "
                           "at the last stage does not have a solution." ) );
 else {
  solver->get_var_solution();

  // TODO make SDDPGreedySolver a CDASolver?
  for( Index t = 0 ; t < get_time_horizon() ; ++t ) {
   auto solver = get_sub_solver( t );
   if( auto cda_solver = dynamic_cast< CDASolver * >( solver ) ) {
    //assert( cda_solver->has_dual_solution() );
    if( cda_solver->has_dual_solution() )
     cda_solver->get_dual_solution();
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::configure_inner_block( Index stage ) {

 if( v_inner_block_configured[ stage ] && v_inner_solver_configured[ stage ] )
  return;

 auto benders_function = get_benders_function( stage );
 auto inner_block = benders_function->get_inner_block();

 // BlockConfig

 if( ! v_inner_block_configured[ stage ] ) {

  if( ( ! f_inner_block_config ) &&
      ( ! f_inner_block_config_filename.empty() ) ) {
   auto c = Configuration::deserialize( f_inner_block_config_filename );
   if( ! ( f_inner_block_config = dynamic_cast< BlockConfig * >( c ) ) ) {
    delete c;
    throw( std::invalid_argument
           ( "SDDPGreedySolver::configure_inner_block: file " +
             f_inner_block_config_filename + " is not a BlockConfig." ) );
   }
  }

  if( f_inner_block_config ) {
   f_inner_block_config->apply( inner_block );
   v_inner_block_configured[ stage ] = true;
  }
 }

 // BlockSolverConfig

 if( ! v_inner_solver_configured[ stage ] ) {

  if( ( ! f_inner_block_solver_config ) &&
      ( ! f_inner_block_solver_config_filename.empty() ) ) {
   auto c = Configuration::deserialize( f_inner_block_solver_config_filename );
   if( ! ( f_inner_block_solver_config =
           dynamic_cast< BlockSolverConfig * >( c ) ) ) {
    delete c;
    throw( std::invalid_argument
           ( "SDDPGreedySolver::configure_inner_block: file " +
             f_inner_block_solver_config_filename +
             " is not a BlockSolverConfig." ) );
   }
  }

  if( f_inner_block_solver_config ) {
   f_inner_block_solver_config->apply( inner_block );
   v_inner_solver_configured[ stage ] = true;
  }
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::unregister_solver_inner_block( Index stage ) {

 auto benders_function = get_benders_function( stage );
 auto inner_block = benders_function->get_inner_block();

 // BlockSolverConfig

 BlockSolverConfig * inner_block_solver_config = nullptr;

 if( v_BSC.size() > stage && v_BSC[ stage ] )
  inner_block_solver_config = v_BSC[ stage ]->clone();
 else if( f_inner_block_solver_config )
  inner_block_solver_config = f_inner_block_solver_config->clone();

 if( inner_block_solver_config ) {
  inner_block_solver_config->clear();
  inner_block_solver_config->apply( inner_block );
  delete inner_block_solver_config;
 }
 else {
  inner_block->unregister_Solvers();
 }

 v_inner_solver_configured[ stage ] = false;
}

/*--------------------------------------------------------------------------*/

int SDDPGreedySolver::solve( Index stage , bool write_solution ) {

 auto benders_function = get_benders_function( stage );

 auto status = benders_function->compute();

 auto solver = benders_function->get_solver();

 if( write_solution ) {
  if( solver->has_var_solution() )
   solver->get_var_solution();
  if( auto cda_solver = dynamic_cast< CDASolver * >( solver ) )
   if( cda_solver->has_dual_solution() )
    cda_solver->get_dual_solution();
 }

 return status;
}

/*--------------------------------------------------------------------------*/

Solver * SDDPGreedySolver::get_sub_solver( Index stage ) const {
 auto benders_function = get_benders_function( stage );
 return benders_function->get_solver();
}

/*--------------------------------------------------------------------------*/

BendersBFunction * SDDPGreedySolver::get_benders_function( Index stage ) const {

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPGreedySolver::get_benders_function: "
                                "invalid stage index: " +
                                std::to_string( stage ) ) );

 auto benders_block = static_cast< BendersBlock * >
  ( static_cast< SDDPBlock * >( f_Block )->
    get_sub_Block( stage )->get_inner_block() );

 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );

 return static_cast< BendersBFunction * >( objective->get_function() );
}

/*--------------------------------------------------------------------------*/

std::vector<double> SDDPGreedySolver::get_solution
( SDDPBlock::Index stage ) const {

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPGreedySolver::get_solution: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );

 Index solution_size = 0;
 for( Index i = 0 ;
      i < sddp_block->get_num_polyhedral_function_per_sub_block() ; ++i ) {
  solution_size +=
   sddp_block->get_polyhedral_function( stage , i )->get_num_active_var();
 }

 std::vector<double> solution;
 solution.reserve( solution_size );

 for( Index i = 0 ;
      i < sddp_block->get_num_polyhedral_function_per_sub_block() ; ++i ) {
  const auto polyhedral_function =
   sddp_block->get_polyhedral_function( stage , i );

  for( const auto & variable : * polyhedral_function ) {
   solution.push_back
    ( static_cast< const ColVariable & >( variable ).get_value() );
  }
 }
 return solution;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_state( const std::vector<double> & state ,
                                  Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPGreedySolver::set_state: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_state( state , stage , 0 );
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::process_outstanding_Modification( void ) {
 while( ! v_mod.empty() ) {
  auto mod = v_mod.front();  // pick (a reference to) the first Modification
  v_mod.pop_front();
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_scenario( Index stage ) {
 static_cast< SDDPBlock * >( f_Block )->set_scenario( scenario_id , stage );
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::store_subgradients( Index stage ,
                                           Index scenario_index ) {

 // First, we store the initial state.

 std::vector< double > initial_state;

 if( stage == 0 ) {
  // Store the initial state of the first stage problem.

  if( ! v_initial_state.empty() )
   // Use the initial state given as parameter to the SDDPGreedySolver.
   initial_state = v_initial_state;
  else {
   // Use the initial state given by SDDPBlock.
   initial_state = static_cast< SDDPBlock * >( f_Block )->get_initial_state();
  }

  if( ! initial_state.empty() )
   f_subgradients.store_initial_state( initial_state , 0 );
 }

 // Store the solution at the given stage as the initial state of the next
 // stage.

 f_subgradients.store_initial_state( get_solution( stage ) , stage + 1 );

 // Store the subgradient with respect to the final state.
 store_subgradient_final_state( stage );

 if( stage > 0 ) {
  // Store the subgradient with respect to the initial state.
  store_subgradient_initial_state( stage , scenario_index ,
                                   get_solution( stage - 1 ) );
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::store_subgradient_final_state( Index stage ) {

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );

 // There must be a single PolyhedralFunction per sub-Block.
 assert( sddp_block->get_num_polyhedral_function_per_sub_block() == 1 );

 // The solution is already written into the Block, so there is no need to set
 // the values of the active Variables of the PolyhedralFunction.

 // Find the index of the active row.

 auto polyhedral_function = sddp_block->get_polyhedral_function( stage );
 polyhedral_function->compute();

 std::vector< double > subgradient;

 if( polyhedral_function->has_linearization() ) {
  // The PolyhedralFunction has a linearization and, therefore, there is a row
  // that is active. By optimality conditions, a subgradient of the objective
  // function if given by *minus* the cut of the future cost function that is
  // active at the final state.
  subgradient.resize( polyhedral_function->get_num_active_var() );
  polyhedral_function->get_linearization_coefficients( subgradient.data() );
 }
 else if( polyhedral_function->get_value() ==
          polyhedral_function->get_global_bound() ) {
  // The bound is active, so the subgradient is zero.
  subgradient.assign( polyhedral_function->get_num_active_var() , 0.0 );
 }

 if( ! subgradient.empty() ) {
  // Store the subgradient.
  f_subgradients.store_subgradient_final_state( std::move( subgradient ) ,
                                                stage );
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::store_subgradient_initial_state
( Index stage , Index scenario_index ,
  const std::vector< double > & initial_state ) {

 if( stage == 0 )
  return;

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );

 // Get the random cut associated with the previous stage and the given
 // scenario index.

 auto & polyhedral_function = sddp_block->get_random_cut
  ( stage - 1 , scenario_index );

 if( polyhedral_function.get_num_active_var() == 0 ) {
  // The PolyhedralFunction representing the random cut has no
  // Variables. Thus, collect the active Variables of the PolyhedralFunction
  // representing the future cost function at the given stage and make them
  // active Variables of the random cut.
  const auto future_cost_function =
   sddp_block->get_polyhedral_function( stage );
  PolyhedralFunction::VarVector active_variables
   ( future_cost_function->get_num_active_var() );
  for( Index i = 0 ; i < future_cost_function->get_num_active_var() ; ++i )
   active_variables[ i ] = static_cast< ColVariable * >
    ( future_cost_function->get_active_var( i ) );

  polyhedral_function.set_variables( std::move( active_variables ) );
 }

 /* In order to evaluate the PolyhedralFunction at the initial state, we must
  * set the values of its active Variables to be equal to the initial
  * state. Thus, we save the current values of the active Variables in order
  * to undo this modification at the end. */

 std::vector< double > original_variable_values;
 original_variable_values.reserve( polyhedral_function.get_num_active_var() );

 auto initial_state_it = std::cbegin( initial_state );

 for( auto & variable : polyhedral_function ) {
  auto & col_variable = static_cast< ColVariable & >( variable );

  // Save the current value of the active Variable.
  original_variable_values.push_back( col_variable.get_value() );

  // Change the value of the active Variable.
  col_variable.set_value( * ( initial_state_it++ ) );
 }

 // Compute the PolyhedralFunction and obtain the index of the active row.

 polyhedral_function.compute();

 std::vector< double > subgradient;

 if( polyhedral_function.has_linearization() ) {
  // The PolyhedralFunction has a linearization and, therefore, there is a row
  // that is active. By optimality conditions, a subgradient of the objective
  // function if given by *minus* the cut of the future cost function that is
  // active at the final state.
  subgradient.resize( polyhedral_function.get_num_active_var() );
  polyhedral_function.get_linearization_coefficients( subgradient.data() );
 }
 else if( polyhedral_function.get_value() ==
          polyhedral_function.get_global_bound() ) {
  // The bound is active, so the subgradient is zero.
  subgradient.assign( polyhedral_function.get_num_active_var() , 0.0 );
 }

 if( ! subgradient.empty() ) {
  // Store the subgradient.
  f_subgradients.store_subgradient_initial_state
   ( std::move( subgradient ) , stage );
 }

 // Put back the original values of the active Variables.

 auto original_variable_values_it = std::cbegin( original_variable_values );
 for( auto & variable : polyhedral_function ) {
  auto & col_variable = static_cast< ColVariable & >( variable );
  col_variable.set_value( * ( original_variable_values_it++ ) );
 }

}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::output_subgradients
( const std::string & filename ) const {

 if( filename.empty() )
  return;

 std::ofstream file( filename );

 if( ! file.is_open() )
  throw( std::runtime_error( "SDDPGreedySolver::output_subgradients: it was "
                             "not possible to open the file \"" +
                             filename + "\"." ) );

 const auto separator = ",";

 file << get_time_horizon() << separator
      << f_subgradients.initial_states.size() << std::endl;

 // Initial states

 for( const auto & state : f_subgradients.initial_states ) {
  file << state.first;
  for( const auto & component : state.second )
   file << separator << component;
  file << std::endl;
 }

 // Subgradients with respect to the initial state

 for( const auto & subgradient : f_subgradients.subgradients_initial_state ) {
  file << subgradient.first << separator << "I";
  for( const auto & component : subgradient.second )
   file << separator << component;
  file << std::endl;
 }

 // Subgradients with respect to the final state

 for( const auto & subgradient : f_subgradients.subgradients_final_state ) {
  file << subgradient.first << separator << "F";
  for( const auto & component : subgradient.second )
   file << separator << component;
  file << std::endl;
 }

 // Finally, we output the scenarios

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );
 const auto & scenario_set = sddp_block->get_scenario_set();

 for( Index stage = 0 ; stage < get_time_horizon() ; ++stage ) {

  auto scenario_begin = scenario_set.sub_scenario_begin( scenario_id , stage );
  auto scenario_end = scenario_set.sub_scenario_end( scenario_id , stage );

  file << stage;

  for( auto scenario_it = scenario_begin ; scenario_it != scenario_end ;
       ++scenario_it )
   file << separator << *scenario_it;
  file << std::endl;
 }

 file.close();
}

/*--------------------------------------------------------------------------*/
/*------------------- End File SDDPGreedySolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/

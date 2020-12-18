/*--------------------------------------------------------------------------*/
/*----------------------- File SDDPGreedySolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPGreedySolver class.
 *
 * \version 0.10
 *
 * \date 04 - 10 - 2020
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
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int SDDPGreedySolver::compute( bool changedvars ) {

 process_outstanding_Modification();
 set_scenario();
 set_initial_state();

 auto time_horizon = get_time_horizon();
 status_compute = Solver::kLowPrecision;
 fault_stage = Inf<Index>();
 solution_value = 0.0;

 for( Index stage = 0 ; stage < time_horizon ; ++stage ) {

  if( stage > 0 ) {
   set_state( get_solution( stage - 1 ) , stage );
  }

  auto sub_status = solve( stage , true );

  if( sub_status == Solver::kInfeasible ) {
   fault_stage = stage;
   status_compute = ( stage == 0 ) ? kInfeasible : kSubproblemInfeasible;
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
  }
 }

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

int SDDPGreedySolver::solve( Index stage , bool write_solution ) {

 auto benders_function = get_benders_function( stage );

 auto status = benders_function->compute();

 auto solver = benders_function->get_solver();

 if( solver->has_var_solution() && write_solution ) {
  solver->get_var_solution();
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

 const auto polyhedral_function =
  static_cast< SDDPBlock * >( f_Block )->get_polyhedral_functions()[ stage ];

 std::vector<double> solution( polyhedral_function->get_num_active_var() );

 Index i = 0;
 for( const auto & variable : * polyhedral_function ) {
  solution[ i++ ] = static_cast< const ColVariable & >( variable ).get_value();
 }
 return solution;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_state( const std::vector<double> & state ,
                                  Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPGreedySolver::set_state: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_state( state , stage );
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::process_outstanding_Modification( void ) {
 while( ! v_mod.empty() ) {
  auto mod = v_mod.front();  // pick (a reference to) the first Modification
  v_mod.pop_front();
  scenario_is_set = false;
  initial_state_is_set = false;
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_scenario( void ) {
 if( ! scenario_is_set ) {
  static_cast< SDDPBlock * >( f_Block )->set_scenario( scenario_id );
  scenario_is_set = true;
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_initial_state( void ) {
 if( ! initial_state_is_set ) {
  static_cast< SDDPBlock * >( f_Block )->set_admissible_state( 0 );
  initial_state_is_set = true;
 }
}

/*--------------------------------------------------------------------------*/
/*------------------- End File SDDPGreedySolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------- File SDDPSimulationSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SDDPSimulationSolver class, a CDASolver for the
 * SDDPBlock that trains the policy of the SDDPBlock and simulates it on
 * every scenario, its value being the mean of the costs of the simulations.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato, Antonio Frangioni,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SDDPSimulationSolver
 #define __SDDPSimulationSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"
#include "SDDPBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class SDDPGreedySolver;  // forward declaration of SDDPGreedySolver
 class SDDPSolver;        // forward declaration of SDDPSolver

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SDDPSimulationSolver_CLASSES Classes in SDDPSimulationSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS SDDPSimulationSolver -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// trains the policy of an SDDPBlock and simulates it on every scenario
/** The SDDPSimulationSolver estimates the expected cost of the policy of an
 * SDDPBlock: it trains the policy with the SDDPSolver registered with the
 * same SDDPBlock (if any; without one the cuts the SDDPBlock already has
 * are the policy), then simulates it on each scenario of the SDDPBlock
 * with the SDDPGreedySolver registered there, and its value [see
 * get_var_value()] is the mean of the costs of the simulations. It does not
 * replace the two Solver, which it uses as they are configured: it has to
 * be registered with the SDDPBlock before them, so that whoever takes the
 * first Solver of the SDDPBlock (say, the InvestmentFunction whose inner
 * Block is the SDDPBlock) takes it.
 *
 * Each scenario is simulated on a copy of the SDDPBlock, made out of the
 * image of the SDDPBlock as the training leaves it (cuts included) and of
 * the BlockConfig of the Block of each of its stages, by a new
 * SDDPGreedySolver with the ComputeConfig and the callback [see
 * SDDPGreedySolver::set_callback()] of the one registered with the
 * SDDPBlock: the SDDPGreedySolver passes to each stage the final state of
 * the previous one, which a unit whose Variable are generated cannot always
 * take, while the units of a copy are generated only when the simulation
 * gets to their stage. The copies are simulated in parallel (OpenMP), and
 * the SDDPBlock has to have a single sub-Block per stage, the copies making
 * more than one of no use.
 *
 * After each scenario the solution of the copy, primal and dual, is written
 * in it, and the events of type eEverySolution [see
 * ThinComputeInterface::event_type] are called, one simulation at a time:
 * meanwhile get_solution_Block() returns the copy, a twin of the SDDPBlock
 * whose data are those of the scenario, so that an event handler can read
 * what it needs out of the solution of each scenario (say, the
 * linearization of the cost with respect to some data of the stages) and
 * combine it as the value is.
 *
 * With MPI, the training runs on every process, and the cuts that the
 * process of rank 0 has at the end are given to every other one, so that
 * every process simulates all the scenarios on the same policy and calls
 * the same events: whoever uses the SDDPSimulationSolver finds the same
 * results on every process, with no MPI of its own.
 *
 * The SDDPSimulationSolver has no solution of its own to write in the
 * SDDPBlock: get_var_solution() and get_dual_solution() do nothing, the
 * solutions being those of the scenarios, given by the events. */

class SDDPSimulationSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: it does nothing but making room for the events
 SDDPSimulationSolver( void ) : CDASolver() {
  v_events.resize( max_event_number() );
  }

/*--------------------------------------------------------------------------*/
 /// the events supported: eEverySolution, and those before it

 [[nodiscard]] EventID max_event_number( void ) const override {
  return( eEverySolution + 1 );
  }

/*--------------------------------------------------------------------------*/
 /// destructor: it does nothing

 ~SDDPSimulationSolver() override = default;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// trains the policy and simulates it on every scenario
 /** If an SDDPSolver is registered with the SDDPBlock, the policy is
  * trained by it first, after removing the cuts of the SDDPBlock if
  * \p changedvars is true (the data of the stages, the investment say, have
  * changed, and the cuts with them). Then every scenario is simulated [see
  * the class description], the events of type eEverySolution being called
  * after each of them.
  *
  * @return kOK if the training succeeded (an iteration limit is a failure:
  *         there is no policy to simulate) and every scenario has been
  *         simulated; otherwise the status of the first Solver that failed,
  *         or that of the event handler that asked to stop, kError if it
  *         was not one of the "stop" ones. */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
 /// the weights of the Block of the stages are applied
 /** As SDDPGreedySolver::handles_weight_change(), the stages being solved by
  * their own Solver. */

 [[nodiscard]] bool handles_weight_change( void ) const override {
  return( true );
  }

 [[nodiscard]] bool handles_variable_weight( void ) const override {
  return( true );
  }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// true if every scenario has been simulated
 bool has_var_solution( void ) override { return( f_has_solution ); }

 /// does nothing: the solutions are those of the scenarios [see compute()]
 void get_var_solution( Configuration * solc = nullptr ) override {}

 /// true if every scenario has been simulated
 bool has_dual_solution( void ) override { return( f_has_solution ); }

 /// true if every scenario has been simulated
 bool is_dual_feasible( void ) override { return( f_has_solution ); }

 /// does nothing: the solutions are those of the scenarios [see compute()]
 void get_dual_solution( Configuration * solc = nullptr ) override {}

 /// the mean of the costs of the simulations of the scenarios
 OFValue get_var_value( void ) override { return( f_value ); }

 /// the mean of the costs of the simulations, for a maximization
 OFValue get_lb( void ) override {
  if( f_has_solution && ( f_sense == Objective::eMax ) )
   return( f_value );
  return( - Inf< OFValue >() );
  }

 /// the mean of the costs of the simulations, for a minimization
 OFValue get_ub( void ) override {
  if( f_has_solution && ( f_sense == Objective::eMin ) )
   return( f_value );
  return( Inf< OFValue >() );
  }

/*--------------------------------------------------------------------------*/
 /// the copy of the SDDPBlock a scenario has just been simulated on
 /** While an event of type eEverySolution is being called, the copy of the
  * SDDPBlock holding the solution of the scenario just simulated [see the
  * class description]; get_Block() otherwise. */

 [[nodiscard]] Block * get_solution_Block( void ) const override {
  return( f_solution_block ? f_solution_block : get_Block() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// the first registered Solver of the SDDPBlock of class T, nullptr if none
 template< class T >
 T * registered( void ) const {
  for( auto solver : f_Block->get_registered_solvers() )
   if( auto s = dynamic_cast< T * >( solver ) )
    return( s );
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// gives the cuts of the process of rank 0 to every other one (MPI only)

 void share_cuts( SDDPSolver * trainer );

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED FIELDS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 OFValue f_value = 0;            ///< the mean cost of the simulations

 bool f_has_solution = false;    ///< whether every scenario was simulated

 int f_sense = Objective::eMin;  ///< the sense of the problems of the stages

 Block * f_solution_block = nullptr;
 ///< the copy simulated last, while its eEverySolution events are called

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class SDDPSimulationSolver )

/** @} end( group( SDDPSimulationSolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* SDDPSimulationSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File SDDPSimulationSolver.h --------------------*/
/*--------------------------------------------------------------------------*/

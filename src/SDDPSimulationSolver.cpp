/*--------------------------------------------------------------------------*/
/*--------------------- File SDDPSimulationSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPSimulationSolver class.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SDDPSimulationSolver.h"

#include "BendersBFunction.h"
#include "BendersBlock.h"
#include "ColRowSolution.h"
#include "FRealObjective.h"
#include "RBlockConfig.h"
#include "SDDPGreedySolver.h"
#include "SDDPSolver.h"
#include "SMSnetCDF.h"

#include <memory>
#include <string>
#include <vector>

#ifdef USE_MPI
 #include <boost/mpi/collectives.hpp>
 #include <boost/mpi/communicator.hpp>
 #include <boost/serialization/vector.hpp>
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( SDDPSimulationSolver );

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL HELPERS -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

// the Block of stage t of sddp_block, the inner Block of the BendersBFunction
// of the Objective of its sub-Block

Block * stage_block( SDDPBlock * sddp_block , Block::Index t )
{
 auto benders_block = sddp_block->get_sub_Block( t , 0 )->get_inner_block();
 return( static_cast< BendersBFunction * >( static_cast< FRealObjective * >(
          benders_block->get_objective() )->get_function() )->
         get_inner_block() );
 }

/*--------------------------------------------------------------------------*/

/* The copies of an SDDPBlock on which the scenarios are simulated: the image
 * of the SDDPBlock and the BlockConfig of the Block of each of its stages,
 * as the training leaves them. A copy has that BlockConfig and a single
 * Solver, an SDDPGreedySolver with the ComputeConfig and the callback of the
 * given one, which configures the Block of each stage when it gets there and
 * takes its Solver away when it leaves the copy. */

struct SimulationImage {
 std::vector< unsigned char > image;
 std::vector< CRBlockConfig > stage_bc;

 explicit SimulationImage( SDDPBlock * sddp_block ) {
  sddp_block->Block::serialize( image );
  const auto num_stages = sddp_block->get_time_horizon();
  stage_bc.reserve( num_stages );
  for( Block::Index t = 0 ; t < num_stages ; ++t )
   stage_bc.emplace_back( stage_block( sddp_block , t ) );
  }

 // a new copy, with the SDDPGreedySolver that simulates it in solver
 SDDPBlock * copy( const SDDPGreedySolver * original ,
                   SDDPGreedySolver * & solver ) const {
  auto c = static_cast< SDDPBlock * >( Block::deserialize( image ) );
  if( ! c )
   throw( std::logic_error( "SDDPSimulationSolver::compute: the image of the "
                            "SDDPBlock cannot be read back" ) );
  solver = nullptr;
  try {
   for( Block::Index t = 0 ; t < stage_bc.size() ; ++t ) {
    // apply() takes the content of what it applies, hence a clone
    auto bc = stage_bc[ t ].clone();
    bc->apply( stage_block( c , t ) );
    delete bc;
    }
   solver = static_cast< SDDPGreedySolver * >(
                              Solver::new_Solver( original->classname() ) );
   c->register_Solver( solver );
   auto cc = original->get_ComputeConfig();
   solver->set_ComputeConfig( cc );
   delete cc;
   solver->set_callback( original->get_callback() );
   }
  catch( ... ) {
   if( solver )
    release( c , solver );
   else
    delete c;
   solver = nullptr;
   throw;
   }
  return( c );
  }

 static void release( SDDPBlock * c , SDDPGreedySolver * solver ) {
  c->unregister_Solver( solver , true );
  delete c;
  }
 };

}  // namespace

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF SDDPSimulationSolver --------------------*/
/*--------------------------------------------------------------------------*/

int SDDPSimulationSolver::compute( bool changedvars )
{
 f_has_solution = false;
 f_value = 0;

 const auto sddp_block = dynamic_cast< SDDPBlock * >( f_Block );
 if( ! sddp_block )
  throw( std::logic_error( "SDDPSimulationSolver::compute: the Block is not "
                           "an SDDPBlock" ) );

 if( sddp_block->get_num_sub_blocks_per_stage() != 1 )
  throw( std::logic_error( "SDDPSimulationSolver::compute: the SDDPBlock "
                           "has " + std::to_string(
                            sddp_block->get_num_sub_blocks_per_stage() ) +
                           " sub-Blocks per stage, and the scenarios are "
                           "simulated on copies of an SDDPBlock with one" ) );

 const auto simulator = registered< SDDPGreedySolver >();
 if( ! simulator )
  throw( std::logic_error( "SDDPSimulationSolver::compute: no "
                           "SDDPGreedySolver is registered with the "
                           "SDDPBlock" ) );
 const auto trainer = registered< SDDPSolver >();

 f_sense = stage_block( sddp_block , 0 )->get_objective_sense();

 // the Solver used lock the SDDPBlock: this Solver lends them its identity,
 // and they get theirs back on every way out of here
 struct Lend {
  Solver * solver;
  void * id;
  ~Lend() { if( solver ) solver->set_id( id ); }
  };
 Lend lend_simulator{ simulator , simulator->id() };
 simulator->set_id( id() );
 Lend lend_trainer{ trainer , trainer ? trainer->id() : nullptr };
 if( trainer )
  trainer->set_id( id() );

 // the training: the cuts of other data are no use for the new ones
 if( trainer ) {
  if( changedvars )
   sddp_block->remove_cuts();
  if( const auto status = trainer->compute() ; status != kOK )
   // whatever the reason (an iteration limit included), there is no
   // policy to simulate
   return( kError );
  share_cuts( trainer );
  }

 const Index num_scenarios = sddp_block->get_scenario_set().size();
 if( ! num_scenarios )
  throw( std::logic_error( "SDDPSimulationSolver::compute: the SDDPBlock "
                           "has no scenario" ) );

 const SimulationImage simulation_image( sddp_block );

 int status = kOK;
 bool interrupt = false;
 OFValue total = 0;

 #pragma omp parallel for reduction( + : total )
 for( int scenario = 0 ; scenario < int( num_scenarios ) ; ++scenario ) {

  if( interrupt )
   continue;

  // an exception cannot leave the parallel loop, nor a critical section:
  // it is a failure of the scenario, as a compute() that finds nothing
  SDDPBlock * copy = nullptr;
  SDDPGreedySolver * solver = nullptr;
  std::string failure;
  #pragma omp critical( SDDPSimulationSolver_copy )
  try {
   copy = simulation_image.copy( simulator , solver );
   }
  catch( const std::exception & e ) {
   failure = e.what();
   }

  int scenario_status = kError;
  if( failure.empty() )
   try {
    solver->set_par( SDDPGreedySolver::intScenarioId , scenario );
    scenario_status = solver->compute( true );
    if( solver->has_var_solution() ) {
     // the solution of each stage, primal and dual, into the copy
     solver->get_var_solution();
     solver->get_dual_solution();
     }
    }
   catch( const std::exception & e ) {
    failure = e.what();
    scenario_status = kError;
    }

  if( ! failure.empty() )
   std::cerr << "SDDPSimulationSolver::compute: an error occurred while "
             << "simulating scenario " << scenario << ": '" << failure << "'"
             << std::endl;

  if( ( ! failure.empty() ) || ( ! solver->has_var_solution() ) ) {
   #pragma omp critical( SDDPSimulationSolver )
   if( ! interrupt ) {
    interrupt = true;
    // a provably infeasible scenario says so, any other failure is an error
    status = ( scenario_status == kInfeasible ) ? kInfeasible : kError;
    }
   }
  else {
   const auto value = solver->get_var_value();
   // the events of the scenario, one simulation at a time, with its
   // solution copied into the SDDPBlock
   #pragma omp critical( SDDPSimulationSolver )
   if( ( ! interrupt ) && ( ! v_events[ eEverySolution ].empty() ) ) {
    try {
     ColRowSolution solution;
     solution.read( copy );
     solution.write( get_Block() );
     }
    catch( const std::exception & e ) {
     std::cerr << "SDDPSimulationSolver::compute: the solution of scenario "
               << scenario << " could not be written in the SDDPBlock: '"
               << e.what() << "'" << std::endl;
     interrupt = true;
     status = kError;
     }
    if( ! interrupt )
     for( auto & event : v_events[ eEverySolution ] ) {
      const auto action = event();
      if( action > eContinue ) {
       interrupt = true;
       status = ( action <= eStopOK ) ? action : kError;
       break;
       }
      }
    }
   total += value;
   }

  if( copy ) {
   #pragma omp critical( SDDPSimulationSolver_copy )
   SimulationImage::release( copy , solver );
   }
  }  // end( for each scenario )

 if( interrupt )
  return( status );

 f_value = total / num_scenarios;
 f_has_solution = true;
 return( kOK );
 }

/*--------------------------------------------------------------------------*/

void SDDPSimulationSolver::share_cuts( SDDPSolver * trainer )
{
#ifdef USE_MPI
 /* With MPI the cuts reach the SDDPBlock of a process through the problems
  * that process solves in the training, and the process of rank 0 is the
  * only one sure to solve some at every stage: its cuts, its State, are
  * given to all the others, which then simulate on the same policy. */

 boost::mpi::communicator world;
 if( world.size() < 2 )
  return;

 NcMemFile::Image image;
 if( ! world.rank() ) {
  std::unique_ptr< State > state( trainer->get_State() );
  NcMemFile file;
  state->serialize( file );
  file.close( image );
  }
 boost::mpi::broadcast( world , image , 0 );
 if( world.rank() ) {
  NcMemFile file( image );
  std::unique_ptr< State > state( State::new_State( file ) );
  trainer->put_State( std::move( *state ) );
  }
#endif
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File SDDPSimulationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/

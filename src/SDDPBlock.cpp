/*--------------------------------------------------------------------------*/
/*--------------------------- File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPBlock class.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claude Opus 4.7 \n
 *         Antrophic \n
 *
 * \copyright &copy; by Rafael Durbano Lobato, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractPath.h"
#include "BendersBlock.h"
#include "SDDPBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( SDDPBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS of SDDPBlock -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------- CONSTRUCTING AND DESTRUCTING SDDPBlock ------------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::deserialize_random_cuts( const std::string & filename )
{
 if( filename.empty() )
  return;

 netCDF::NcFile file( filename.c_str() , netCDF::NcFile::read );

 const auto TimeHorizon = file.getDim( "TimeHorizon" );
 if( TimeHorizon.isNull() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the dimension TimeHorizon "
           "was not provided." ) );

 const auto time_horizon = TimeHorizon.getSize();

 if( time_horizon != get_time_horizon() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the expected TimeHorizon "
           "dimension is " + std::to_string( get_time_horizon() ) +
           ", but " + std::to_string( time_horizon ) + " was given." ) );

 const auto NumberScenarios = file.getDim( "NumberScenarios" );

 if( NumberScenarios.isNull() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the dimension "
           "NumberScenarios was not provided." ) );

 const auto number_scenarios = NumberScenarios.getSize();

 if( number_scenarios != size() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the expected NumberScenarios "
           "dimension is " + std::to_string( size() ) +
           ", but " + std::to_string( number_scenarios ) + " was given." ) );

 // Possibly clear the previous random cuts
 random_cuts.resize( boost::extents[ 0 ][ 0 ] );

 // Create the random cuts
 random_cuts.resize( boost::extents[ time_horizon ][ number_scenarios ] );

 for( Index t = 0 ; t < time_horizon ; ++t ) {

  // Collect the active Variables of the PolyhedralFunction at stage t
  const auto polyhedral_function = get_polyhedral_function( t );
  PolyhedralFunction::VarVector active_variables
   ( polyhedral_function->get_num_active_var() );
  for( Index i = 0 ; i < polyhedral_function->get_num_active_var() ; ++i )
   active_variables[ i ] = static_cast< ColVariable * >
    ( polyhedral_function->get_active_var( i ) );

  for( Index s = 0 ; s < number_scenarios ; ++s ) {
   // Set the active Variables of the PolyhedralFunction
   auto variables = active_variables;
   random_cuts[ t ][ s ].set_variables( std::move( variables ) );

   // Deserialize the PolyhedralFunction (if provided)
   auto group_name = "PolyhedralFunction_" +
    std::to_string( t ) + "_" + std::to_string( s );
   auto group = file.getGroup( group_name );
   if( group.isNull() )
    continue;
   random_cuts[ t ][ s ].deserialize( group );
  }
 }
}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_Modification( sp_Mod mod , Observer::ChnlName chnl )
{
 // TODO
 if( anyone_there() )
  Block::add_Modification( std::make_shared< NBModification >( this ) ,
			   chnl );
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/

int SDDPBlock::get_objective_sense( void ) const
{
 try {
  auto sub_Block = get_sub_Block( 0 );
  if( sub_Block )
   return( sub_Block->get_objective_sense() );
  }
 catch( ... ) {}
 return( Objective::eUndef );
 }

/*--------------------------------------------------------------------------*/

StochasticBlock * SDDPBlock::get_sub_Block( Index stage ,
					    Index sub_block_index ) const
{
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::get_sub_Block: invalid stage " +
                                std::to_string( stage ) ) );
 if( sub_block_index >= num_sub_blocks_per_stage )
  throw( std::invalid_argument(
			"SDDPBlock::get_sub_Block: invalid sub-Block index "
			+ std::to_string( sub_block_index ) ) );
 const auto index = stage * num_sub_blocks_per_stage + sub_block_index;
 return( static_cast< StochasticBlock * >( v_Block[ index ] ) );
 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_cuts( PolyhedralFunction::MultiVector && A ,
                          PolyhedralFunction::RealVector && b , Index stage ,
                          Index sub_block_index ,
                          Index number_cuts_to_keep )
{
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::add_cuts: invalid stage index: " +
                                std::to_string( stage ) ) );

 auto polyhedral_function =
  get_polyhedral_function( stage , 0 , sub_block_index );

 const auto num_rows = polyhedral_function->get_nrows();

 // Possibly remove the last (num_rows - number_cuts_to_keep) cuts
 if( number_cuts_to_keep < num_rows )
  polyhedral_function->delete_rows( Range( number_cuts_to_keep , num_rows) );

 // Add the given cuts
 polyhedral_function->add_rows( std::move( A ) , b );
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::store_random_cut( std::vector< double > && coefficients ,
                                  double alpha , Index stage ,
                                  Index scenario_index )
{
 assert( stage < get_time_horizon() );
 assert( scenario_index < size() );

 if( ! f_random_cuts_initialized ) {
  // Create the PolyhedralFunctions that will store the random cuts.

  // Ensure the random cuts are initialized by only one thread.
#pragma omp critical (SDDPBlock_random_cut)
  {
   if( ! f_random_cuts_initialized ) {
    initialize_random_cuts();
    f_random_cuts_initialized = true;
   }
  }
 }

 // Store the given random cut.

 random_cuts[ stage ][ scenario_index ].add_row( std::move( coefficients ) ,
                                                 alpha );
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::initialize_random_cuts( void )
{
 const auto time_horizon = get_time_horizon();
 const auto number_scenarios = size();
 random_cuts.resize( boost::extents[ 0 ][ 0 ] );
 random_cuts.resize( boost::extents[ time_horizon ][ number_scenarios ] );

 for( Index t = 0 ; t < time_horizon ; ++t ) {

  const auto polyhedral_function = get_polyhedral_function( t );
  PolyhedralFunction::VarVector active_variables
   ( polyhedral_function->get_num_active_var() );
  for( Index i = 0 ; i < polyhedral_function->get_num_active_var() ; ++i )
   active_variables[ i ] = static_cast< ColVariable * >
    ( polyhedral_function->get_active_var( i ) );

  for( Index s = 0 ; s < number_scenarios ; ++s ) {
   auto variables = active_variables;
   random_cuts[ t ][ s ].set_variables( std::move( variables ) );
   random_cuts[ t ][ s ].set_is_convex( polyhedral_function->is_convex() );
  }
 }
}

/*--------------------------------------------------------------------------*/

double SDDPBlock::get_future_cost( Index stage , Index sub_block_index )
 const
{
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::get_future_cost: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 auto function = get_polyhedral_function( stage , 0 , sub_block_index );
 if( ! function )
  // If this SDDPBlock has no PolyhedralFunction, the future cost must be 0
  return( 0 );
 function->compute();
 return( function->get_value() );
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const Eigen::ArrayXd & values , Index stage ,
                           Index sub_block_index )
{
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 benders_block->set_variable_values( values );
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const std::vector< double > & values ,
			   Index stage , Index sub_block_index )
{
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 benders_block->set_variable_values( values );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > SDDPBlock::get_state( Index stage ,
                                            Index sub_block_index ) const
{
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >(
	  get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 return( benders_block->get_variable_values() );
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_admissible_state( Index stage , Index sub_block_index )
{
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );

 auto benders_block = static_cast< BendersBlock * >(
	   get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );

 auto admissible_state = get_admissible_state( stage );
 benders_block->set_variable_values( admissible_state );
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR PRINTING & SAVING THE SDDPBlock ---------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::print( std::ostream & output , char vlvl ) const
{
 output << std::endl << "SDDPBlock with ";

 if( v_Block.empty() )
  output << "no inner Block";
 else
  output << v_Block.size() << " sub-Blocks" << std::endl;
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::prepare_generator_pool( void )
{
 if( ! f_scenario_generator )
  throw( std::logic_error( "SDDPBlock::prepare_generator_pool: "
                           "no ScenarioGenerator is attached to this "
                           "SDDPBlock." ) );

 if( ! f_scenario_generator->is_pool_initialized() )
  throw( std::logic_error( "SDDPBlock::prepare_generator_pool: the "
                           "ScenarioGenerator's pool has not been "
                           "initialized; call init_random_pool() or "
                           "init_representative_pool() on the generator "
                           "before invoking this method." ) );

 // make sure we are at the start of the pool
 f_scenario_generator->reset_pool();

 // sanity check on the per-stage decomposition: scenario_set must have its
 // structural metadata already set (either from deserialize() in the
 // legacy path or from set_structural_metadata() in the new ScenarioGenerator
 // path); get_scenario_size() must match the generator's get_scenario_size()
 if( scenario_set.get_time_horizon() == 0 )
  throw( std::logic_error( "SDDPBlock::prepare_generator_pool: the "
                           "structural metadata of scenario_set is not "
                           "set; this is normally done at deserialize "
                           "time. Cannot proceed." ) );

 const auto expected_scenario_size = scenario_set.get_scenario_size();
 if( f_scenario_generator->get_scenario_size() != expected_scenario_size )
  throw( std::logic_error( "SDDPBlock::prepare_generator_pool: the "
                           "ScenarioGenerator's scenario_size (" +
                           std::to_string(
                             f_scenario_generator->get_scenario_size() ) +
                           ") does not match the expected SDDPBlock "
                           "scenario_size (" +
                           std::to_string( expected_scenario_size ) + ")." ) );

 // iterate the pool, collecting scenarios and probabilities directly
 // into the SDDPBlock-side cache (not into scenario_set.scenarios)
 f_generator_pool_cache.clear();
 std::vector< double > probabilities;

 while( true ) {
  const auto s = f_scenario_generator->get_current_scenario();
  f_generator_pool_cache.emplace_back( s.begin() , s.end() );
  probabilities.push_back(
                f_scenario_generator->get_current_scenario_probability() );
  if( ! f_scenario_generator->next_scenario() )
   break;
  }

 const auto pool_size = f_generator_pool_cache.size();
 if( pool_size == 0 )
  throw( std::logic_error( "SDDPBlock::prepare_generator_pool: the "
                           "ScenarioGenerator returned an empty pool." ) );

 // verify uniform probabilities: v2 step 1 still only supports the case
 // where every scenario has weight 1 / pool_size (within a small
 // tolerance). Non-uniform weights would require weighted-average cuts
 // in the backward pass of SDDPSolver, which is deferred (cf. v2+ TODO).
 const double uniform = 1.0 / double( pool_size );
 constexpr double tol = 1e-9;
 for( std::size_t i = 0 ; i < probabilities.size() ; ++i ) {
  if( std::abs( probabilities[ i ] - uniform ) > tol ) {
   // wipe the partial cache so further calls see an unprepared state
   f_generator_pool_cache.clear();
   throw( std::logic_error(
    "SDDPBlock::prepare_generator_pool: non-uniform scenario "
    "probabilities are not yet supported (scenario " + std::to_string( i ) +
    " has probability " + std::to_string( probabilities[ i ] ) +
    ", while a uniform distribution would assign " +
    std::to_string( uniform ) + " to every scenario). Weighted-average "
    "cuts in the SDDP backward pass are planned for a future version." ) );
   }
  }

 // mirror the pool size into scenario_set so legacy callers that read
 // get_scenario_set().size() observe the right count (the .scenarios
 // storage of scenario_set is intentionally left empty in this path)
 scenario_set.set_num_scenarios( static_cast< Index >( pool_size ) );

 // leave the generator at the beginning of the pool
 f_scenario_generator->reset_pool();

 }  // end( SDDPBlock::prepare_generator_pool )

/*--------------------------------------------------------------------------*/

void SDDPBlock::prepare_multi_stage_generator_pool( void )
{
 if( ! f_scenario_generator )
  throw( std::logic_error( "SDDPBlock::prepare_multi_stage_generator_pool: "
                           "no ScenarioGenerator is attached to this "
                           "SDDPBlock." ) );

 auto * mgen =
  dynamic_cast< MultiStageScenarioGenerator * >( f_scenario_generator );
 if( ! mgen )
  throw( std::logic_error( "SDDPBlock::prepare_multi_stage_generator_pool: "
                           "the attached ScenarioGenerator is not a "
                           "MultiStageScenarioGenerator." ) );

 if( ! mgen->is_pool_initialized() )
  throw( std::logic_error( "SDDPBlock::prepare_multi_stage_generator_pool: "
                           "the MultiStageScenarioGenerator's pool has not "
                           "been initialized; call init_random_pool() or "
                           "init_representative_pool() on the generator "
                           "before invoking this method." ) );

 // Cross-check that the generator's stage count matches the SDDPBlock's
 // time horizon (the latter being derived from the number of sub-Blocks
 // and the SubScenarioSize metadata in scenario_set).
 const auto T = scenario_set.get_time_horizon();
 if( T == 0 )
  throw( std::logic_error( "SDDPBlock::prepare_multi_stage_generator_pool: "
                           "the structural metadata of scenario_set is not "
                           "set; this is normally done at deserialize "
                           "time. Cannot proceed." ) );

 const auto gen_T = mgen->get_stage_number();
 if( static_cast< Index >( gen_T ) != T )
  throw( std::logic_error( "SDDPBlock::prepare_multi_stage_generator_pool: "
                           "the MultiStageScenarioGenerator's stage number ("
                           + std::to_string( gen_T ) +
                           ") does not match the SDDPBlock's time horizon ("
                           + std::to_string( T ) + ")." ) );

 // Reset the stage cursor (and the iteration of stage 0) and prepare
 // the cache. reset_pool() would only rewind the current stage's
 // iteration; to also rewind the stage cursor we use the
 // previous_stage( INFStage ) idiom — see the comments on
 // MultiStageScenarioGenerator in ScenarioGenerator.h.
 mgen->previous_stage( MultiStageScenarioGenerator::INFStage );
 f_multi_stage_pool_cache.assign( T , {} );

 // Walk strategy (v2 step 2 framework, stage-independent assumption):
 //
 // For each stage t in 0..T-1, navigate from a freshly-reset pool down
 // to stage t (taking the default-history descent), then enumerate
 // next_scenario() until exhausted to collect every x_t at this node.
 //
 // This is correct under the assumption that the per-stage pool does
 // not depend on the path history H_t — exactly the case the planned
 // MultiStageScenarioSet (one DiscreteScenarioSet per stage) will
 // implement. Generalising to history-dependent multi-stage trees would
 // require a full depth-first walk and a different cache layout (a
 // forest rather than a per-stage list); we defer that to a follow-up
 // when an actual non-independent subclass exists to motivate the
 // design.
 //
 // Uniform probabilities are required *within each stage's pool*: the
 // weights returned by get_current_scenario_probability() must all be
 // equal to 1 / N_t (with a small tolerance). The N_t values may
 // differ across stages.

 for( Index t = 0 ; t < T ; ++t ) {
  // Navigate to stage t from the start: previous_stage( INFStage )
  // rewinds the cursor to stage 0 (and resets stage 0's iteration);
  // a plain reset_pool() would only touch the current stage.
  mgen->previous_stage( MultiStageScenarioGenerator::INFStage );
  for( Index u = 0 ; u < t ; ++u ) {
   if( ! mgen->next_stage() )
    throw( std::logic_error(
     "SDDPBlock::prepare_multi_stage_generator_pool: "
     "MultiStageScenarioGenerator refused to advance to stage "
     + std::to_string( u + 1 ) + " while preparing stage "
     + std::to_string( t ) + "." ) );
   }

  // Sanity-check the per-stage scenario_size matches sub_scenario_size[t].
  const auto expected_sz = scenario_set.get_sub_scenario_size( t );
  if( mgen->get_scenario_size() != expected_sz )
   throw( std::logic_error(
    "SDDPBlock::prepare_multi_stage_generator_pool: "
    "MultiStageScenarioGenerator's scenario_size at stage "
    + std::to_string( t ) + " ("
    + std::to_string( mgen->get_scenario_size() )
    + ") does not match the SDDPBlock-side sub_scenario_size ("
    + std::to_string( expected_sz ) + ")." ) );

  // Enumerate x_t at this node.
  std::vector< double > probabilities;
  while( true ) {
   const auto s = mgen->get_current_scenario();
   f_multi_stage_pool_cache[ t ].emplace_back( s.begin() , s.end() );
   probabilities.push_back( mgen->get_current_scenario_probability() );
   if( ! mgen->next_scenario() )
    break;
   }

  const auto N_t = f_multi_stage_pool_cache[ t ].size();
  if( N_t == 0 ) {
   f_multi_stage_pool_cache.clear();
   throw( std::logic_error(
    "SDDPBlock::prepare_multi_stage_generator_pool: "
    "MultiStageScenarioGenerator returned an empty pool at stage "
    + std::to_string( t ) + "." ) );
   }

  const double uniform = 1.0 / double( N_t );
  constexpr double tol = 1e-9;
  for( std::size_t k = 0 ; k < probabilities.size() ; ++k )
   if( std::abs( probabilities[ k ] - uniform ) > tol ) {
    f_multi_stage_pool_cache.clear();
    throw( std::logic_error(
     "SDDPBlock::prepare_multi_stage_generator_pool: non-uniform "
     "scenario probabilities are not yet supported (at stage "
     + std::to_string( t ) + ", scenario " + std::to_string( k ) +
     " has probability " + std::to_string( probabilities[ k ] ) +
     ", while a uniform distribution over " + std::to_string( N_t ) +
     " scenarios would assign " + std::to_string( uniform ) + " each).") );
    }
  }

 // Mirror an aggregate scenario count into scenario_set so that legacy
 // callers that still read get_scenario_set().size() observe something
 // sensible. In the stage-independent multi-stage case there is no
 // single "number of full scenarios" (it would be the *product* of the
 // per-stage N_t's), so we report stage 0's pool size: this matches the
 // value returned by size() (without the stage parameter) and avoids
 // surprising callers that conflate the two.
 scenario_set.set_num_scenarios(
             static_cast< Index >( f_multi_stage_pool_cache[ 0 ].size() ) );

 // leave the generator at the beginning of the pool (stage 0, scenario
 // 0): previous_stage( INFStage ) rewinds both the stage cursor and
 // stage 0's iteration; reset_pool() alone would only rewind the
 // iteration of whatever stage we ended up on.
 mgen->previous_stage( MultiStageScenarioGenerator::INFStage );

 }  // end( SDDPBlock::prepare_multi_stage_generator_pool )

/*--------------------------------------------------------------------------*/

void SDDPBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 // type

 group.putAtt( "type" , "SDDPBlock" );

 // TimeHorizon

 const auto time_horizon = get_time_horizon();
 auto TimeHorizon_dim = group.addDim( "TimeHorizon" , time_horizon );

 // NumSubBlocksPerStage

 group.addDim( "NumSubBlocksPerStage" , num_sub_blocks_per_stage );

 // StochasticBlock_i

 for( Index i = 0 ; i < get_time_horizon() ; ++i ) {
  auto sub_group = group.addGroup( "StochasticBlock_" + std::to_string( i ) );
  get_sub_Block( i )->serialize( sub_group );
 }

 // AbstractPaths to PolyhedralFunctions

 std::vector< AbstractPath > paths;
 paths.reserve( get_time_horizon() * num_polyhedral_per_sub_block );

 for( Index t = 0 ; t < get_time_horizon() ; ++t ) {
  for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i ) {
   auto reference_block = get_sub_Block( t )->get_nested_Block( 0 );
   assert( reference_block );
   paths.emplace_back( get_polyhedral_function( t , i ) , reference_block );
  }
 }

 AbstractPath::serialize( paths , group );

 if( num_polyhedral_per_sub_block != 1 )
  group.addDim( "NumPolyhedralFunctionsPerSubBlock" ,
                num_polyhedral_per_sub_block );

 // Scenarios
 //
 // In the legacy path (f_scenario_generator == nullptr) the scenarios live
 // inside scenario_set and are serialized via ScenarioSet::serialize().
 //
 // In the new path (f_scenario_generator != nullptr) we serialize the
 // ScenarioGenerator under a dedicated "ScenarioGenerator" sub-group and
 // the structural metadata (SubScenarioSize / NumberRandomDataGroups /
 // SizeRandomDataGroups) at the SDDPBlock level for symmetry with the
 // legacy layout. We do NOT serialize the (possibly empty / stale)
 // contents of scenario_set in this case: scenario_set is recreated from
 // the generator at the next deserialize+load cycle.

 if( f_scenario_generator ) {
  auto gen_group = group.addGroup( "ScenarioGenerator" );
  f_scenario_generator->serialize( gen_group );

  // mirror ScenarioSet's structural metadata at SDDPBlock level (these
  // are problem-structure metadata owned by SDDPBlock, not by the
  // generator)
  std::vector< Index > sss( scenario_set.get_time_horizon() );
  for( Index t = 0 ; t < scenario_set.get_time_horizon() ; ++t )
   sss[ t ] = scenario_set.get_sub_scenario_size( t );
  ::SMSpp_di_unipi_it::serialize( group , "SubScenarioSize" ,
                                  netCDF::NcUint() , TimeHorizon_dim ,
                                  sss , false );

  const auto & rd = scenario_set.get_size_random_data_groups();
  if( ! rd.empty() ) {
   auto NRDG_dim = group.addDim( "NumberRandomDataGroups" , rd.size() );
   ::SMSpp_di_unipi_it::serialize( group , "SizeRandomDataGroups" ,
                                   netCDF::NcUint() , NRDG_dim , rd ,
                                   false );
   }
  }
 else
  scenario_set.serialize( group );

 // Initial state

 auto InitialStateSize = group.addDim( "InitialStateSize" ,
                                       initial_state.size() );

 ::SMSpp_di_unipi_it::serialize( group , "InitialState" ,
                                 netCDF::NcDouble() , InitialStateSize ,
                                 initial_state , false );

 // StateSize

 std::vector< Index > state_size( time_horizon );
 for( Index t = 0 ; t < time_horizon - 1 ; ++t )
  state_size[ t ] = admissible_state_begin[ t + 1 ] -
                    admissible_state_begin[ t ];
 if( time_horizon > 0 )
  state_size.back() = admissible_states.size() - admissible_state_begin.back();

 ::SMSpp_di_unipi_it::serialize( group , "StateSize" , netCDF::NcUint64() ,
                                 TimeHorizon_dim , state_size , false );

 // AdmissibleState

 auto AdmissibleState_dim = group.addDim( "AdmissibleState_dim" ,
                                          admissible_states.size() );

 ::SMSpp_di_unipi_it::serialize( group , "AdmissibleState" ,
                                 netCDF::NcDouble() , AdmissibleState_dim ,
                                 admissible_states , false );
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::serialize_random_cuts( const std::string & filename ) const
{
 if( filename.empty() || ( random_cuts.num_elements() == 0 ) )
  return;

 netCDF::NcFile file( filename , netCDF::NcFile::replace );

 const auto time_horizon = get_time_horizon();
 file.addDim( "TimeHorizon" , time_horizon );

 const auto number_scenarios = random_cuts[ 0 ].size();
 file.addDim( "NumberScenarios" , number_scenarios );

 for( Index t = 0 ; t < time_horizon ; ++t )
  for( Index s = 0 ; s < number_scenarios ; ++s ) {
   auto group_name = "PolyhedralFunction_" +
    std::to_string( t ) + "_" + std::to_string( s );
   auto group = file.addGroup( group_name );
   random_cuts[ t ][ s ].serialize( group );
   }
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

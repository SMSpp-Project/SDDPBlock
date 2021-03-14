/*--------------------------------------------------------------------------*/
/*--------------------------- File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPBlock class.
 *
 * \version 0.10
 *
 * \date 09 - 03 - 2021
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractPath.h"
#include "BendersBlock.h"
#include "SDDPBlock.h"
#include "StochasticBlock.h"

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

Block * SDDPBlock::deserialize_sub_Block( const netCDF::NcGroup & group ,
                                          Index i ) {

 std::string sub_group_name = "StochasticBlock_" + std::to_string( i );
 auto sub_group = group.getGroup( sub_group_name );

 if( sub_group.isNull() ) {
  sub_group = group.getGroup( "StochasticBlock" );
  if( sub_group.isNull() )
   throw std::logic_error( "SDDPBlock::deserialize: neither group '" +
                           sub_group_name + "' nor 'StochasticBlock' "
                           "was found." );
  sub_group_name = "StochasticBlock";
 }

 auto type = sub_group.getAtt( "type" );
 if( type.isNull() )
  throw std::logic_error( "SDDPBlock::deserialize: attribute 'type' of '" +
                          sub_group_name + "' must be present." );

 std::string type_name;
 type.getValues( type_name );

 if( type_name != "StochasticBlock" )
  throw std::logic_error( "SDDPBlock::deserialize: attribute 'type' of '" +
                          sub_group_name + "' must contain "
                          "'StochasticBlock'." );

 auto sub_Block = new_Block( sub_group , this );

 if( ! sub_Block )
  throw std::logic_error( "SDDPBlock::deserialize: sub-group '" +
                          sub_group_name + "' is incomplete." );

 if( sub_group_name != "StochasticBlock" ) {

  // If StochasticBlock_i does not have sub-group "Block" then
  // "StochasticBlock" must have one.
  if( sub_group.getGroup( "Block" ).isNull() ) {

   auto StochasticBlock_group = group.getGroup( "StochasticBlock" );
   if( StochasticBlock_group.isNull() )
    throw std::logic_error( "SDDPBlock::deserialize: sub-group 'Block' was not "
                            "provided neither in '" + sub_group_name +
                            "' nor in 'StochasticBlock'" );


   auto Block_group = StochasticBlock_group.getGroup( "Block" );
   if( Block_group.isNull() )
    throw std::logic_error( "SDDPBlock::deserialize: sub-group 'Block' was not "
                            "provided neither in '" + sub_group_name +
                            "' nor in 'StochasticBlock'" );

   auto inner_block = new_Block( Block_group, this );
   if( ! inner_block )
    throw std::logic_error( "SDDPBlock::deserialize: the 'Block' sub-group of "
                            "the 'StochasticBlock' group has an invalid or "
                            "incomplete description." );

   static_cast< StochasticBlock * >( sub_Block )->
    set_inner_block( inner_block );
  }

  // If StochasticBlock_i does not have the description of vector of
  // "DataMapping" then, if "StochasticBlock" has one, we use it.

  Index num_data_mappings;
  if( ! ::SMSpp_di_unipi_it::deserialize_dim( sub_group , "NumberDataMappings" ,
                                              num_data_mappings , true ) ) {

   auto StochasticBlock_group = group.getGroup( "StochasticBlock" );
   if( ! StochasticBlock_group.isNull() ) {

    if( ::SMSpp_di_unipi_it::deserialize_dim( StochasticBlock_group ,
                                              "NumberDataMappings" ,
                                              num_data_mappings , true ) ) {

     std::vector< std::unique_ptr< SimpleDataMappingBase > > data_mappings;
     data_mappings.reserve( num_data_mappings );
     SimpleDataMappingBase::deserialize
      ( group , data_mappings ,
        static_cast< StochasticBlock *>( sub_Block )->get_inner_block() );

     static_cast< StochasticBlock * >( sub_Block )->
      set_data_mappings( std::move( data_mappings ) );
    }
   }
  }
 }

 return sub_Block;
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::deserialize( const netCDF::NcGroup & group ) {

 // TimeHorizon

 Index time_horizon;
 ::SMSpp_di_unipi_it::deserialize_dim( group , "TimeHorizon" ,
                                       time_horizon , false );

 ::SMSpp_di_unipi_it::deserialize_dim( group , "TimeHorizon" ,
                                       num_sub_blocks_per_stage , true );

 // NumSubBlocksPerStage

 if( ! ::SMSpp_di_unipi_it::deserialize_dim
     ( group , "NumSubBlocksPerStage" , num_sub_blocks_per_stage ) ) {
  num_sub_blocks_per_stage = 1;
 }

 // StochasticBlock

 v_Block.reserve( time_horizon * num_sub_blocks_per_stage );

 for( Index i = 0 ; i < time_horizon ; ++i )
  for( Index j = 0 ; j < num_sub_blocks_per_stage ; ++j )
   v_Block.push_back( deserialize_sub_Block( group , i ) );

 // PolyhedralFunctions

 auto path_group = group.getGroup( "AbstractPath" );

 auto paths = AbstractPath::vector_deserialize( path_group );

 if( ! ::SMSpp_di_unipi_it::deserialize_dim
     ( group , "NumPolyhedralFunctionsPerSubBlock" ,
       num_polyhedral_per_sub_block ) ) {
  num_polyhedral_per_sub_block = 1;
 }

 if( paths.size() != num_polyhedral_per_sub_block * time_horizon &&
     ! ( paths.size() == num_polyhedral_per_sub_block && time_horizon > 1 ) ) {
  if( num_polyhedral_per_sub_block == 1 )
   throw ( std::invalid_argument
           ( "SDDPBlock::deserialize: The number of AbstractPath to "
             "PolyhedralFunction must be either equal to 1 or equal to "
             "the time horizon." ) );
  else
   throw ( std::invalid_argument
           ( "SDDPBlock::deserialize: The number of AbstractPath to "
             "PolyhedralFunction must be either equal to K or equal to K "
             "times the time horizon, where K is the number of "
             "PolyhedralFunction per sub-Block." ) );
 }

 v_polyhedral_functions.clear();
 v_polyhedral_functions.reserve
  ( num_sub_blocks_per_stage * num_polyhedral_per_sub_block * time_horizon );

 for( Index t = 0 ; t < time_horizon ; ++t ) {
  for( Index j = 0 ; j < num_sub_blocks_per_stage ; ++j ) {
   auto reference_block = get_sub_Block( t , j )->get_nested_Block( 0 );
   assert( reference_block );
   for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i ) {
    Index path_index = num_polyhedral_per_sub_block * t + i;
    if( paths.size() == num_polyhedral_per_sub_block )
     path_index = i;
    auto polyhedral_function = dynamic_cast< PolyhedralFunction * >
     ( paths[ path_index ].get_element< Function >( reference_block ) );
    if( ! polyhedral_function )
     throw ( std::invalid_argument
             ( "SDDPBlock::deserialize: PolyhedralFunction for stage "
               + std::to_string( t ) + " was not found." ) );
    v_polyhedral_functions.push_back( polyhedral_function );
   }
  }
 }

 // Scenarios

 scenario_set.deserialize( group );

 // Initial state

 ::SMSpp_di_unipi_it::deserialize( group , "InitialState" ,
                                   initial_state , false );

 // StateSize

 std::vector< Index > state_size;

 ::SMSpp_di_unipi_it::deserialize( group , "StateSize" , { time_horizon } ,
                                   state_size , false , true );

 bool state_size_is_scalar = ( state_size.size() == 1 );
 if( state_size.size() == 1 )
  state_size.resize( time_horizon , state_size[ 0 ] );
 else if( state_size.size() != time_horizon )
  throw ( std::logic_error( "SDDPBlock::deserialize: 'StateSize' must be "
                            "either a scalar or an array with size "
                            "'TimeHorizon'." ) );

 // AdmissibleState

 ::SMSpp_di_unipi_it::deserialize( group , "AdmissibleState" ,
                                   admissible_states , false );

 if( state_size_is_scalar ) {
  if( admissible_states.size() != state_size[ 0 ] &&
      admissible_states.size() != time_horizon * state_size[ 0 ] )
   throw ( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                             "array has an invalid size." ) );

  if( admissible_states.size() != time_horizon * state_size[ 0 ] ) {
   std::vector<double> state = admissible_states;
   admissible_states.reserve( time_horizon * state_size[ 0 ] );
   for( Index t = 1 ; t < time_horizon ; ++t )
    admissible_states.insert( admissible_states.cend() ,
                              state.cbegin() , state.cend() );
  }
 }
 else if( admissible_states.size() !=
          std::accumulate( state_size.begin() , state_size.end() ,
                           decltype( state_size )::value_type( 0 ) ) ) {
  throw ( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                            "array has an invalid size." ) );
 }

 // Construct the vector admissible_state_begin

 admissible_state_begin.resize( time_horizon );
 if( time_horizon > 0 )
  admissible_state_begin.front() = 0;
 for( Index t = 1 ; t < time_horizon ; ++t )
  admissible_state_begin[ t ] =
   admissible_state_begin[ t - 1 ] + state_size[ t - 1 ];

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_Modification( sp_Mod mod , Observer::ChnlName chnl ) {
 // TODO
 if( anyone_there() )
  Block::add_Modification( std::make_shared<NBModification>( this ) , chnl );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/

StochasticBlock * SDDPBlock::get_sub_Block
( Index stage , Index sub_block_index ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::get_sub_Block: invalid stage " +
                                std::to_string( stage ) ) );
 if( sub_block_index >= num_sub_blocks_per_stage )
  throw( std::invalid_argument( "SDDPBlock::get_sub_Block: invalid sub-Block in"
                                "dex " + std::to_string( sub_block_index ) ) );
 const auto index = stage * num_sub_blocks_per_stage + sub_block_index;
 return static_cast< StochasticBlock * >( v_Block[ index ] );
}

/*--------------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_cuts( PolyhedralFunction::MultiVector && A ,
                          PolyhedralFunction::RealVector && b , Index stage ,
                          Index number_cuts_to_keep ) {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::add_cuts: invalid stage index: " +
                                std::to_string( stage ) ) );

 const auto num_rows = v_polyhedral_functions[ stage ]->get_nrows();

 // Possibly remove the last (num_rows - number_cuts_to_keep) cuts
 if( number_cuts_to_keep < num_rows )
  v_polyhedral_functions[ stage ]->
   delete_rows( Range( number_cuts_to_keep , num_rows) );

 // Add the given cuts
 v_polyhedral_functions[ stage ]->add_rows( std::move( A ) , b );
}

/*--------------------------------------------------------------------------*/

double SDDPBlock::get_future_cost( Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::get_future_cost: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 v_polyhedral_functions[ stage ]->compute();
 return v_polyhedral_functions[ stage ]->get_value();
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const Eigen::ArrayXd & values , Index stage ,
                           Index sub_block_index ) {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 benders_block->set_variable_values( values );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const std::vector<double> & values , Index stage ,
                           Index sub_block_index ) {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 benders_block->set_variable_values( values );
}

/*--------------------------------------------------------------------------*/

std::vector< double > SDDPBlock::get_state( Index stage ,
                                            Index sub_block_index ) const {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 return benders_block->get_variable_values();
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_admissible_state( Index stage , Index sub_block_index ) {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );

 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );

 auto admissible_state = get_admissible_state( stage );
 benders_block->set_variable_values( admissible_state );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_scenario( Index scenario_id , Index stage ,
                              Index sub_block_index ) {
 auto sub_scenario_begin = scenario_set.
  sub_scenario_begin( scenario_id , stage );

 try {
  get_sub_Block( stage , sub_block_index )->set_data( sub_scenario_begin );
 }
 catch( const std::exception & e ) {
  std::cout << "SDDPBlock::set_scenario: exception while setting scenario "
            << scenario_id << " of stage " << stage << ".\n"
            << e.what() << std::endl;
  std::exit( EXIT_FAILURE );
 }
}

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE SDDPBlock ----------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::serialize( netCDF::NcGroup & group ) const {

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
  auto sub_group = group.addGroup( "StochasticBlock_" +
                                   std::to_string( i ) );
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
  state_size[ t ] = admissible_state_begin[ t+1 ] - admissible_state_begin[ t ];
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

void SDDPBlock::print( std::ostream &output ) const {
 output << std::endl << "SDDPBlock with ";

 if( v_Block.empty() )
  output << "no inner Block";
 else
  output << v_Block.size() << " sub-Blocks" << std::endl;
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

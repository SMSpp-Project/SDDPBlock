/*--------------------------------------------------------------------------*/
/*--------------------------- File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPBlock class.
 *
 * \version 0.10
 *
 * \date 06 - 03 - 2020
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

Block * SDDPBlock::deserialize_sub_Block( netCDF::NcGroup & group , Index i ) {

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

 auto type = group.getAtt( "type" );
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

void SDDPBlock::deserialize( netCDF::NcGroup & group ) {

 // TimeHorizon

 Index time_horizon;
 ::SMSpp_di_unipi_it::deserialize_dim( group , "TimeHorizon" ,
                                       time_horizon , false );

 // StochasticBlock

 v_Block.reserve( time_horizon );

 for( Index i = 0 ; i < time_horizon ; ++i )
  v_Block.push_back( deserialize_sub_Block( group , i ) );

 // PolyhedralFunctions

 auto path_group = group.getGroup( "AbstractPath" );

 auto paths = AbstractPath::vector_deserialize( path_group );

 if( paths.size() != time_horizon &&
     ! ( paths.size() == 1 && time_horizon > 1 ) ) {
  throw ( std::invalid_argument
          ( "SDDPBlock::deserialize: The number of AbstractPath to "
            "PolyhedralFunction must be either equal to 1 or equal to "
            "the time horizon." ) );
 }

 v_polyhedral_functions.clear();
 v_polyhedral_functions.reserve( time_horizon );

 for( Index i = 0 ; i < time_horizon ; ++i ) {
  auto reference_block = static_cast< StochasticBlock * >( v_Block[ i ] )->
   get_nested_Blocks().front();
  assert( reference_block );
  auto path_index = ( paths.size() == 1 ) ? 0 : i;
  auto polyhedral_function = dynamic_cast< PolyhedralFunction * >
   ( AbstractPath::get_element< Function >( paths[ path_index ] ,
                                            reference_block ) );
  if( ! polyhedral_function )
   throw ( std::invalid_argument( "SDDPBlock::deserialize: PolyhedralFunction "
                                  + std::to_string( i ) + " was not found." ) );
  v_polyhedral_functions.push_back( polyhedral_function );
 }

 // NumberScenarios

 ::SMSpp_di_unipi_it::deserialize_dim( group , "NumberScenarios" ,
                                       num_scenarios , false );

 // ScenarioSize

 ::SMSpp_di_unipi_it::deserialize_dim( group , "ScenarioSize" ,
                                       scenario_size , false );

 // SubScenarioSize

 if( ::SMSpp_di_unipi_it::deserialize( group , "SubScenarioSize" ,
                                       time_horizon , sub_scenario_size ,
                                       true , false ) ) {

  if( scenario_size !=
      std::accumulate( sub_scenario_size.begin() ,
                       sub_scenario_size.end() , 0 ) )
   throw ( std::logic_error( "SDDPBlock::deserialize: The sum of the "
                             "elements in 'SubScenarioSize' must be equal to "
                             "'ScenarioSize'" ) );
 }
 else {
  // SubScenarioSize was not provided

  if( scenario_size % time_horizon != 0 )
   throw ( std::logic_error( "SDDPBlock::deserialize: 'SubScenarioSize' was "
                             "not provided. Thus, 'ScenarioSize' must be a "
                             "multiple of 'TimeHorizon'." ) );

  sub_scenario_size.resize( time_horizon , scenario_size / time_horizon );
 }

 // Scenarios

 ::SMSpp_di_unipi_it::deserialize( group , "Scenarios" , scenarios ,
                                   false , false );

 if( scenarios.shape()[ 0 ] != num_scenarios ||
     scenarios.shape()[ 1 ] != scenario_size ) {

  throw ( std::logic_error( "SDDPBlock::deserialize: 'Scenarios' must be a "
                            "two-dimensional array whose first and second "
                            "dimensions have sizes 'NumberScenarios' and "
                            "'ScenarioSize', respectively." ) );
 }

 // NumberRandomDataGroups and SizeRandomDataGroups

 if( std::adjacent_find( sub_scenario_size.begin() , sub_scenario_size.end() ,
                         std::not_equal_to<>() ) != sub_scenario_size.end() ) {
  // Not all sub-scenarios have the same size. In this case, we consider a
  // single random data group.
  num_random_data_groups = 1;
 }
 else {

  // All sub-scenarios have the same size. Thus, we check
  // NumberRandomDataGroups and SizeRandomDataGroups.

  if( ! ::SMSpp_di_unipi_it::deserialize_dim
      ( group , "NumberRandomDataGroups" , num_random_data_groups , true ) ) {
   // NumberRandomDataGroups was not provided. Hence, there must be a single
   // random data group.
   num_random_data_groups = 1;
  }
  else {
   // NumberRandomDataGroups was provided. Now, we check SizeRandomDataGroups.

   if( ::SMSpp_di_unipi_it::deserialize
       ( group , "SizeRandomDataGroups" , num_random_data_groups ,
         size_random_data_groups , true , false ) ) {

    // SizeRandomDataGroups was provided.

    if( ( scenario_size / time_horizon ) !=
        std::accumulate( size_random_data_groups.begin() ,
                         size_random_data_groups.end() , 0 ) )
     throw ( std::logic_error( "SDDPBlock::deserialize: The sum of the sizes "
                               "in 'SizeRandomDataGroups' must be equal to "
                               "'ScenarioSize' / 'TimeHorizon'." ) );
   }
   else {
    // SizeRandomDataGroups was not provided.
    size_random_data_groups = { ( scenario_size / time_horizon ) };
    if( num_random_data_groups > 1 ) {
     throw ( std::logic_error( "SDDPBlock::deserialize: 'NumberRandomData"
                               "Groups' must be provided since "
                               "'NumberRandomDataGroups' > 1." ) );
    }
   }
  }
 }

 // StateSize

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
  if( admissible_states.size() != state_size[ 0 ] ||
      admissible_states.size() != time_horizon * state_size[ 0 ] )
   throw ( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                             "array has an invalid size." ) );
 }
 else if( admissible_states.size() !=
          std::accumulate( state_size.begin() , state_size.end() , 0 ) ) {
  throw ( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                            "array has an invalid size." ) );
 }

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_Modification( sp_Mod mod , Observer::ChnlName chnl ) {
 // TODO
 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/

StochasticBlock * SDDPBlock::get_sub_Block( Index i ) const {
 if( i >= v_Block.size() )
  throw( std::invalid_argument( "SDDPBlock::get_sub_Block: invalid sub-Block "
                                "index: " + std::to_string( i ) ) );
 return static_cast< StochasticBlock * >( v_Block[ i ] );
}

/*--------------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::update_cuts( PolyhedralFunction::MultiVector && A ,
                             PolyhedralFunction::RealVector & b ,
                             Index stage ) {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::update_cuts: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 v_polyhedral_functions[ stage ]->add_rows( std::move( A ) , b );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const Eigen::ArrayXd & values , Index stage ) {
 assert( stage < get_time_horizon() );
 auto benders_block = static_cast< BendersBlock * >
  ( static_cast< StochasticBlock * >( v_Block[ stage ] )->
    get_nested_Blocks().front() );
 benders_block->set_variable_values( values );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_scenario( const Eigen::ArrayXd & scenario , Index stage ) {
 assert( stage < get_time_horizon() );
 static_cast< StochasticBlock * >( v_Block[ stage ] )->set_data( scenario );
}

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE SDDPBlock ----------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::serialize( netCDF::NcGroup & group ) const {

 Block::serialize( group );

 // type

 group.putAtt( "type" , "SDDPBlock" );

 // TimeHorizon

 auto TimeHorizon_dim = group.addDim( "TimeHorizon" , get_time_horizon() );

 // StochasticBlock_i

 for( Index i = 0 ; i < v_Block.size() ; ++i ) {
  auto sub_group = group.addGroup( "StochasticBlock_" +
                                   std::to_string( i ) );
  v_Block[ 0 ]->serialize( sub_group );
 }

 // AbstractPaths to PolyhedralFunctions

 std::vector< AbstractPath > paths;
 paths.reserve( v_polyhedral_functions.size() );

 for( Index i = 0 ; i < paths.size() ; ++i ) {
  auto reference_block = static_cast< StochasticBlock * >( v_Block[ i ] )->
   get_nested_Blocks().front();
  assert( reference_block );
  paths.push_back( AbstractPath::build_path< PolyhedralFunction >
                   ( v_polyhedral_functions[ i ] , reference_block ) );
 }

 AbstractPath::serialize( paths , group );

 /////////////////////////////////////////////////////////

 // NumberScenarios

 auto NumberScenarios_dim = group.addDim( "NumberScenarios" , num_scenarios );

 // ScenarioSize

 auto ScenarioSize_dim = group.addDim( "ScenarioSize" , scenario_size );

 // SubScenarioSize

 ::SMSpp_di_unipi_it::serialize( group , "SubScenarioSize" ,
                                 netCDF::NcUint64() , TimeHorizon_dim ,
                                 sub_scenario_size , false );

 // Scenarios

 ::SMSpp_di_unipi_it::serialize( group , "Scenarios" , netCDF::NcDouble() ,
                                 { NumberScenarios_dim , ScenarioSize_dim } ,
                                 scenarios , false , false );

 // NumberRandomDataGroups and SizeRandomDataGroups

 auto NumberRandomDataGroups_dim = group.addDim( "NumberRandomDataGroups" ,
                                                 num_random_data_groups );

 if( size_random_data_groups.size() > 0 )
  ::SMSpp_di_unipi_it::serialize( group , "SizeRandomDataGroups" , netCDF::NcUint64() ,
                                  NumberRandomDataGroups_dim ,
                                  size_random_data_groups , false );

 // StateSize

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

/*--------------------------------------------------------------------------*/
/*--------------------------- File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPBlock class.
 *
 * \version 0.10
 *
 * \date 08 - 01 - 2020
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

  // If StochasticBlock_i does not have description of vector of
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
 ::SMSpp_di_unipi_it::deserialize_dim( group, "TimeHorizon",
                                       time_horizon, false );

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

 // TODO
 // NumberScenarios

 // ScenarioSize

 // SubScenarioSize

 // Scenarios

 // NumberRandomDataGroups

 // SizeRandomDataGroups

 // StateSize

 // AdmissibleState
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

 // type

 group.putAtt( "type" , "SDDPBlock" );

 // TimeHorizon

 group.addDim( "TimeHorizon" , get_time_horizon() );

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

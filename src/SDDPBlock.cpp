/*--------------------------------------------------------------------------*/
/*--------------------------- File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPBlock class.
 *
 * \version 0.10
 *
 * \date 23 - 12 - 2019
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
#include "SDDPSolver.h"
#include "StochasticBlock.h"
#include "StOpt/sddp/backwardForwardSDDP.h"
#include "StOpt/sddp/LocalConstRegressionForSDDP.h"
#include "StOpt/sddp/LocalLinearRegressionForSDDP.h"

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

 std::string sub_group_name = "StochasticBlock";
 if( i < Inf<Index>() )
  sub_group_name += "_" + std::to_string( i );

 auto sub_group = group.getGroup( sub_group_name );

 if( sub_group.isNull() )
  throw std::logic_error( "SDDPBlock::deserialize: the '" +
                          sub_group_name + "' was not found." );

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

 return sub_Block;
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::deserialize( netCDF::NcGroup & group ) {

 // TimeHorizon

 Index time_horizon;
 ::deserialize_dim( group, "TimeHorizon", time_horizon, false );

 // StochasticBlock

 v_Block.reserve( time_horizon );

 auto stochastic_block_group = group.getGroup( "StochasticBlock" );

 if( ! stochastic_block_group.isNull() ) {
  // The "StochasticBlock" group is present. Therefore, the sub-Blocks are all
  // identical.
  for( Index i = 0 ; i < time_horizon ; ++i )
   v_Block.push_back( deserialize_sub_Block( group ) );
 }
 else {
  for( Index i = 0 ; i < time_horizon ; ++i )
   v_Block.push_back( deserialize_sub_Block( group , i ) );
 }

 // PolyhedralFunctions

 auto paths = AbstractPath::vector_deserialize( group );

 if( paths.size() != time_horizon )
  throw ( std::invalid_argument( "SDDPBlock::deserialize: The number of Abstrac"
                                 "tPath to PolyhedralFunction must be equal to "
                                 "the time horizon." ) );

 v_polyhedral_functions.clear();
 v_polyhedral_functions.reserve( time_horizon );

 for( Index i = 0 ; i < time_horizon ; ++i ) {
  auto reference_block = static_cast< StochasticBlock * >( v_Block[ i ] )->
   get_inner_block();
  assert( reference_block );
  auto polyhedral_function = dynamic_cast< PolyhedralFunction * >
   ( AbstractPath::get_element< Function >( paths[ i ] , reference_block ) );
  if( ! polyhedral_function )
   throw ( std::invalid_argument( "SDDPBlock::deserialize: PolyhedralFunction "
                                  + std::to_string( i ) + " was not found." ) );
  v_polyhedral_functions.push_back( polyhedral_function );
 }
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
  ( static_cast< StochasticBlock * >( v_Block[ stage ] )->get_inner_block() );
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
  auto reference_block =
   static_cast< StochasticBlock * >( v_Block[ i ] )->get_inner_block();
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

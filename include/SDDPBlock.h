/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of SDDPBlock, a class for representing a multistage stochastic
 * programming problem specifically designed to be solved by an SDDP solver.
 *
 * \version 0.1
 *
 * \date 23 - 12 - 2019
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SDDPBlock
#define __SDDPBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "PolyhedralFunction.h"
#include "ScenarioSimulator.h"
#include "StOpt/sddp/SimulatorSDDPBase.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class StochasticBlock;      // forward declaration of StochasticBlock

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SDDPBlock_CLASSES Classes in SDDPBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SDDPBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// SDDPBlock, representing a multistage stochastic programming problem
/** The SDDPBlock is a class that represents a multistage stochastic
 * programming problem of the form
 *
 * \f[
 *    \min_{x_1 \in X_1} f_1(x_1) +
 *    \mathbb{E} \left \lbrack \min_{x_2 \in X_2(x_1, \xi_2)} f_2(x_2; \xi_2) +
 *    \mathbb{E} \left \lbrack \dots +
 *    \mathbb{E} \left \lbrack \min_{x_{T} \in X_{T}(x_{T-1}, \xi_{T})}
 *                             f_{T}(x_{T}; \xi_{T})
 *    \right\rbrack \right\rbrack\right\rbrack,
 * \f]
 *
 * where \f$ f_t: \mathbb{R}^{n_t} \rightarrow \mathbb{R} \f$ for \f$
 * t \in \{1, \dots, T\} \f$, \f$ \{\xi_{t}\}_{t \in \{2, \dots, T\}}
 * \f$ is a stochastic process, and \f$ T \f$ is the time horizon. For
 * each \f$ t \in \{2, \dots, T\}\f$, the function \f$ f_t \f$ may
 * depend on \f$ \xi_t \f$, while the set \f$ X_{t} \f$ may depend on
 * both \f$ x_{t-1} \f$ and \f$ \xi_t \f$. For each \f$ t \in \{1,
 * \dots, T\} \f$, we call
 *
 * \f[
 *    \min_{x_t \in X_t(x_{t-1}, \xi_t)} f_t(x_t; \xi_t) +
 *         \mathcal{V}_{t+1}(x_t)
 * \f]
 *
 * the problem associated with stage \f$ t \f$, where
 *
 * \f[
 *    \mathcal{V}_{t+1}(x_t) =
 *          \mathbb{E} \left\lbrack V_{t+1}(x_t, \xi_{t+1}) \right\rbrack
 *  \f]
 *
 * is the (expected value) cost-to-go function (also called value
 * function, future value function, future cost function), with \f$
 * \mathcal{V}_{T+1} \equiv 0 \f$,
 *
 * \f[
 *    V_{t}(x_{t-1}, \xi_{t}) =
 *             \min_{x_t \in X_t(x_{t-1}, \xi_{t})}
 *             f_{t}(x_{t}; \xi_{t}) + \mathcal{V}_{t+1}(x_t)
 * \f]
 *
 * and
 *
 * \f[
 *    f_1(x_1) \doteq f_1(x_1; \xi_1) \quad \mbox{ and } \quad
 *    X_1 \doteq X_1(x_{0}, \xi_1)
 * \f]
 *
 * with given \f$ x_0 \in \mathbb{R}^{n_0} \f$ and deterministic \f$
 * \xi_1\f$. We consider an approximation to the problem associated
 * with stage \f$ t \in \{1, \dots, T\} \f$ as the problem
 *
 * \f{equation}{
 *    \min_{x_t \in X_t(x_{t-1}, \xi_t)}
 *    f_t(x_t; \xi_t) + \mathcal{P}_{t+1}(x_t)        \quad       (1)
 * \f}
 *
 * where \f$ \mathcal{P}_{t+1}(x_t) \f$ is a polyhedral
 * function, i.e., it is a function of the form
 *
 * \f[
 *    \mathcal{P}_{t+1}(x_t) = \max_{i \in \{1,\dots,k_t\}}
 *                                     \{ a_{t,i}^{\top}x_t + b_{t,i} \}
 * \f]
 *
 * with \f$ a_{t,i} \in \mathbb{R}^{n_t} \f$ and \f$ b_{t,i} \in
 * \mathbb{R} \f$ for each \f$ i \in \{1,\dots,k_t\} \f$.
 *
 * The SDDPBlock has \f$ T \f$ sub-Blocks, each one being a
 * StochasticBlock. The \f$t\f$-th sub-Block represents an approximation to
 * the problem associated with stage \f$ t \f$ as defined in
 * (1).
 */

class SDDPBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------- CONSTRUCTING AND DESTRUCTING SDDPBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SDDPBlock
 *  @{ */

 /// constructor
 /** Constructs an SDDPBlock with the given \p father Block. The input
  * parameter has a default value (\c nullptr), so that this can be used as
  * the void constructor.
  *
  * @param father A pointer to the father Block of this SDDPBlock.
  */
 SDDPBlock( Block * father = nullptr ) : Block( father ) { }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~SDDPBlock() {
  for( auto & block : v_Block )
   delete block;
  v_Block.clear();
 }

/*--------------------------------------------------------------------------*/
 /// de-serialize an SDDPBlock out of netCDF::NcGroup
 /** The method takes a netCDF::NcGroup supposedly containing all the
  * information required to de-serialize the SDDPBlock. Besides the mandatory
  * "type" attribute of any :Block, the group must contain the following:
  *
  * - The "TimeHorizon" dimension, containing the time horizon.
  *
  * - The description of the sub-Blocks. If the sub-group "StochasticBlock" is
  *   present, then all sub-Blocks are assumed to be identical and the
  *   description in this sub-group is used to construct each of the
  *   sub-Blocks of this SDDPBlock. If "StochasticBlock" is not present, then
  *   the sub-groups "StochasticBlock_i", for i in {0, ..., TimeHorizon - 1}
  *   must be present, where sub-group "StochasticBlock_i" contains the
  *   description of the i-th sub-Block of this SDDPBlock.
  *
  * - All the dimensions and variables necessary to describe a vector of
  *   AbstractPath as described in the comments of
  *   AbstractPath::deserialize(). The number of AbstractPath must be equal to
  *   TimeHorizon and, therefore, the dimension associated with the number of
  *   AbstractPath is TimeHorizon. The i-th AbstractPath in this vector must
  *   be the path to the PolyhedralFunction associated with the i-th sub-Block
  *   of this SDDPBlock. The i-th AbstractPath is taken with respect to the inner
  *   Block of the i-th sub-Block of this SDDPBlock.
  *
  * @param group A netCDF::NcGroup holding the data of the SDDPBlock.
  */

 void deserialize( netCDF::NcGroup & group ) override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR Saving THE DATA OF THE SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Saving the data of the SDDPBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// serialize an SDDPBlock into a netCDF::NcGroup
 /** Serialize an SDDPBlock into a netCDF::NcGroup with the format
  * explained in the comments of the deserialize() function.
  *
  * @param group The NcGroup in which this SDDPBlock will be serialized.
  */

 virtual void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the SDDPBlock
    @{ */

 /// returns the time horizon
 /** This function returns the time horizon associated with this SDDPBlock.
  *
  * @return The time horizon.
  */
 inline virtual std::size_t get_time_horizon() const {
   return v_Block.size();
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of PolyhedralFunction
 /** This function returns the vector of PolyhedralFunction associated with
  * this SDDPBlock.
  *
  * @return The vector of PolyhedralFunction.
  */
 const std::vector< PolyhedralFunction * > &
 get_polyhedral_functions( ) const {
  return v_polyhedral_functions;
 }

/*--------------------------------------------------------------------------*/

 /// returns the i-th sub-Block of this SDDPBlock
 /** This function returns the i-th sub-Block of this SDDPBlock. The given
  * index \p i must be between 0 and get_time_horizon() - 1. If \p i is an
  * invalid index, an exception is thrown.
  *
  * @param i The index of the desired sub-Block. It must be a number between 0
  *        and get_time_horizon() - 1.
  *
  * @return The i-th sub-Block of this SDDPBlock.
  */
 virtual StochasticBlock * get_sub_Block( Index i ) const;

/**@} ----------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods describing the behavior of an SDDPBlock
 * @{ */

 /// update the cuts of the problem at the given stage
 /** This function updates the cuts of the problem at the given \p stage. The
  * parameters must satisfy the following requirements:
  *
  * 1. \p cuts must be a matrix with as many columns as there are cuts to be
  *    added and the number of rows must be equal to the number of Variable
  *    defined in the BendersBlock associated with stage \p stage.
  *
  * 2. \p stage must be an integer between 0 and get_time_horizon() - 1.
  *
  * @param cuts An Eigen::ArrayXXd containing the cuts to be added.
  *
  * @param stage The stage whose cuts should be updated.
  */
 void update_cuts( PolyhedralFunction::MultiVector && A ,
                   PolyhedralFunction::RealVector & b , Index stage );

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem at the
  * given \p stage. The size of the \p values array parameter must be equal to
  * the number N of state Variable of the problem at the given \p stage, so
  * that the value of the i-th state Variable will be values( i ), for each i
  * in {0, ..., N-1}.
  *
  * @param values The Eigen::ArrayXd containing the values of the Variable.
  */
 void set_state( const Eigen::ArrayXd & values , Index stage );

/*--------------------------------------------------------------------------*/

 /// updates the sub-Block at the given stage for the given scenario
 /** This function updates the sub-Block at the given \p stage for the given
  * \p scenario.
  *
  * @param scenario The scenario that must be set.
  */
 void set_scenario( const Eigen::ArrayXd & scenario , Index stage );

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 virtual void print( std::ostream &output ) const override;

/*--------------------------------------------------------------------------*/

 virtual void load( std::istream &input ) override {}

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 /// Pointers to the PolyhedralFunction of each sub-Block
 std::vector< PolyhedralFunction * > v_polyhedral_functions;

 /// Simulator for the forward step of the SDDP method
 std::shared_ptr< ScenarioSimulator > simulator_forward;

 /// Simulator for the backward step of the SDDP method
 std::shared_ptr< ScenarioSimulator > simulator_backward;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// deserializes the i-th sub-Block out of the given group
 /** This auxiliary function deserializes the \p i-th sub-Block out of the
  * given \p group.
  *
  * @param group The netCDF::NcGroup containing the description of the
  *        sub-Block.
  *
  * @param i The index of the sub-Block to be deserialized. If i is
  *        Inf<Index>(), then the sub-Block is deserialized out of the
  *        sub-group named "StochasticBlock". If i < Inf<Index>(), then the
  *        sub-Block is deserialized out of the sub-group named
  *        "StochasticBlock_i".
  *
  * @return A pointer to the Block that was deserialized.
  */
 Block * deserialize_sub_Block( netCDF::NcGroup & group ,
                                Index i = Inf<Index>() );

/*--------------------------------------------------------------------------*/

};   // end( class SDDPBlock )

/** @} end( group( SDDPBlock_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SDDPBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SDDPBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/

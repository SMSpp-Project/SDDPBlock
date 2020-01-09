/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of SDDPBlock, a class for representing a multistage stochastic
 * programming problem specifically designed to be solved by an SDDP solver.
 *
 * \version 0.1
 *
 * \date 09 - 01 - 2020
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
/** The SDDPBlock is a class that derives from Block and represents a
 * multistage stochastic programming problem of the form
 *
 * \f[
 *    \min_{\substack{x_0 \in \mathbb{R}^{n_0} \\ A_0 x_0 + B_0 x_{-1} = b_0\\
 *                    x_0 \ge 0}} c_0^{\top}x_0 +
 *    \mathbb{E} \left \lbrack
 *    \min_{\substack{x_1 \in \mathbb{R}^{n_1} \\ A_1 x_1 + B_1 x_0 = b_1\\
 *                    x_1 \ge 0}} c_1^{\top}x_1 +
 *    \mathbb{E} \left \lbrack \dots +
 *    \mathbb{E} \left \lbrack
 *    \min_{\substack{x_{T-1} \in \mathbb{R}^{n_{T-1}} \\
 *          A_{T-1} x_{T-1} + B_{T-1} x_{T-2} = b_{T-1}\\
 *                    x_{T-1} \ge 0}} c_{T-1}^{\top}x_{T-1}
 *    \right\rbrack \right\rbrack\right\rbrack,
 * \f]
 *
 * where \f$ T \f$ is called the time horizon and \f$ \xi = \{ (b_t,
 * c_t, A_t, B_t) \}_{t \in \{1, \dots, T-1\}} \f$ is a stochastic
 * process. This means that some (or all) the components of the
 * matrices \f$ A_t \f$ and \f$ B_t \f$ and the vectors \f$ b_t \f$
 * and \f$ c_t \f$ may be random variables. Notice that \f$ x_{-1} \f$
 * and \f$ (b_0, c_0, A_0, B_0) \f$, which we denote by \f$ \xi_0 \f$,
 * are deterministic. The term \f$ B_0 x_{-1} \f$ in the first stage
 * problem could be disregarded (i.e., we could have \f$ B_0 = 0 \f$
 * or \f$ x_{-1} = 0 \f$ without loss of generality), but we keep them
 * in order to have all subproblems with the same structure, which
 * will facilitate our approach.
 *
 * For each \f$ t \in \{0, \dots, T-1\}\f$, we call
 *
 * \f[
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t}\\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{V}_{t+1}(x_t, \xi_t)
 * \f]
 *
 * the problem associated with stage \f$ t \f$, where
 *
 * \f[
 *    \mathcal{V}_{t+1}(x_t, \xi_t) =
 *      \mathbb{E}
 *        \left\lbrack
 *          V_{t+1}(x_t, \xi_{t+1}) \mid \xi_t
 *        \right\rbrack
 * \f]
 *
 * is the (expected value) cost-to-go function (also called value
 * function, future value function, future cost function), with \f$
 * \mathcal{V}_{T} \equiv 0 \f$ and
 *
 * \f[
 *    V_{t}(x_{t-1}, \xi_{t}) =
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t} \\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{V}_{t+1}(x_t, \xi_t)
 * \f]
 *
 * with given \f$ x_{-1} \f$ and (deterministic)
 * \f$ \xi_0\f$. We consider an approximation to the problem
 * associated with stage \f$ t \in \{0, \dots, T-1\} \f$ as the problem
 *
 * \f[
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t} \\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{P}_{t+1}(x_t)                            \qquad (1)
 * \f]
 *
 * where \f$ \mathcal{P}_{t+1}(x_t) \f$ is a polyhedral
 * function, i.e., it is a function of the form
 *
 * \f[
 *    \mathcal{P}_{t+1}(x_t) = \max_{i \in \{1,\dots,k_t\}}
 *                                     \{ d_{t,i}^{\top}x_t + e_{t,i} \}
 * \f]
 *
 * with \f$ d_{t,i} \in \mathbb{R}^{n_t} \f$ and \f$ e_{t,i} \in
 * \mathbb{R} \f$ for each \f$ i \in \{1,\dots,k_t\} \f$.
 *
 * An SDDPBlock is characterized by the following:
 *
 * 1) It has a time horizon \f$ T \f$.
 *
 * 2) It has \f$ T \f$ sub-Blocks, each one being a StochasticBlock. The
 *    \f$t\f$-th sub-Block represents an approximation to the problem
 *    associated with stage \f$ t \f$ as defined in (1).
 *
 * 3) It has pointers to \f$ T - 1 \f$ PolyhedralFunction. The \f$t\f$-th
 *    PolyhedralFunction represents the function \f$ \mathcal{P}_{t+1} \f$ in
 *    (1) and, therefore, must be defined in the \f$t\f$-th sub-Block of this
 *    SDDPBlock or in any of the sub-Blocks of that sub-Block, recursively.
 *
 * 4) It has a set of scenarios. Each scenario is represented by a vector of
 *    double and spans all the time horizon \f$ T \f$ . Each vector is divided
 *    into \f$ T \f$ parts, each one being associated with a stage of the
 *    multistage problem. Let \f$ S \f$ denote a vector representing a
 *    scenario. Then, \f$ S \f$ is defined as
 *
 *    \f[
 *        S = ( S_0 , ... , S_{T-1} )
 *    \f]
 *
 *    where \f$ S_t \f$ is a sub-vector of S with size \f$ s_t \f$, for each
 *    \f$ t \in \{ 0, ..., T-1 \} \f$, and is associated with the sub-problem
 *    at stage \f$ t \f$, i.e., it provides data for the \f$t\f$-th sub-Block
 *    of this SDDPBlock. We say that \f$ S_t \f$ represents the \f$t\f$-th
 *    sub-scenario of the scenario represented by \f$ S \f$.
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
  * - The description of the sub-Blocks. This is given by the sub-groups
  *   "StochasticBlock" and "StochasticBlock_t", for each t in {0, ...,
  *   TimeHorizon - 1}. These sub-groups are optional, but they cannot be all
  *   absent. If "StochasticBlock_t" is not provided by some t in {0, ...,
  *   TimeHorizon - 1}, then "StochasticBlock" must be provided and contain a
  *   complete description of the t-th sub-Block of this SDDPBlock. If
  *   "StochasticBlock" is not provided, then "StochasticBlock_t" must be
  *   provided for each t in {0, ..., TimeHorizon - 1} and contain a complete
  *   description of the t-th sub-Block of this SDDPBlock.
  *
  *   If "StochasticBlock_t" is provided but the description of its inner
  *   Block is not provided, then the sub-group "StochasticBlock" must be
  *   provided and contain the description of an inner Block of a
  *   StochasticBlock. In this case, the description of the inner Block
  *   provided in the sub-group "StochasticBlock" will be used to construct
  *   the inner Block of the StochasticBlock described by the
  *   "StochasticBlock_t" sub-group.
  *
  *   If "StochasticBlock_t" is provided but the description of its vector of
  *   DataMapping is not provided, then if the sub-group "StochasticBlock" is
  *   provided and contains a description of a vector of DataMapping, then it
  *   is used to construct the vector of DataMapping of the StochasticBlock
  *   decribed by the sub-group "StochasticBlock_t".
  *
  * - All the dimensions and variables necessary to describe a vector of
  *   AbstractPath as described in the comments of
  *   AbstractPath::deserialize(). The number of AbstractPath must be equal to
  *   "TimeHorizon - 1". The i-th AbstractPath in this vector must be the path
  *   to the PolyhedralFunction associated with the i-th sub-Block of this
  *   SDDPBlock. The i-th AbstractPath is taken with respect to the inner
  *   Block of the i-th sub-Block of this SDDPBlock.
  *
  * - The "SubScenarioSize" variable, of type netCDF::NcUint64 and indexed
  *   over dimension "TimeHorizon". This dimension is optional. If it is not
  *   provided, then all sub-scenarios are assumed to have the same size,
  *   i.e., \f$ s_i = s_j \f$ for all \f$ i,j \in \{ 0, ..., T-1 \}\f$. If it
  *   is present, then SubScenarioSize[t] is the size of the sub-scenario
  *   associated with stage t, i.e., SubScenarioSize[t] = \f$ s_t \f$, for
  *   each \f$ i \in \{ 0, ..., T-1 \}\f$.
  *
  * - The two-dimensional variable "Scenarios" of type netCDF::NcDouble,
  *   containing the scenarios. The i-th row of "Scenarios" contains the i-th
  *   scenario, so that Scenarios[ i ][ j ] is the j-th component of the i-th
  *   scenario.
  *
  * @param group A netCDF::NcGroup holding the data describing this SDDPBlock.
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

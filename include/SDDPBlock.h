/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of SDDPBlock, a class for representing a multistage stochastic
 * programming problem specifically designed to be solved by an SDDP solver.
 *
 * \version 0.1
 *
 * \date 27 - 02 - 2020
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
 *   \min_{x_0 \in \mathcal{X}^{n_0}} f_0(x_0) +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_1 \in \mathcal{X}^{n_1}} f_1(x_1) +
 *   \mathbb{E} \left \lbrack \dots +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{X}^{n_{T-1}}} f_{T-1}(x_{T-1})
 *   \right\rbrack \right\rbrack\right\rbrack,
 * \f]
 *
 * where T is called the time horizon, \f$\mathcal{X}^{n_t} \equiv
 * \mathcal{X}^{n_t}(x_{t-1}, \xi_t) \subseteq \mathbb{R}^{n_t}\f$ for each
 * \f$t \in \{0, \dots, T-1\}\f$, and \f$ \xi = \{ \xi_t \}_{t \in \{1, \dots,
 * T-1\}} \f$ is a stochastic process. Notice that \f$ x_{-1} \f$ and \f$
 * \xi_0 \f$ are deterministic. For each \f$ t \in \{0, \dots, T-1\}\f$, we
 * call
 *
 * \f[
 *   \min_{x_t \in \mathcal{X}^{n_t}} f_t(x_t) +
 *   \mathcal{V}_{t+1}(x_t, \xi_t)
 * \f]
 *
 * the problem associated with stage \f$ t \f$, where
 *
 * \f[
 *   \mathcal{V}_{t+1}(x_t, \xi_t) =
 *    \mathbb{E}
 *      \left\lbrack
 *        V_{t+1}(x_t, \xi_{t+1}) \mid \xi_t
 *      \right\rbrack
 * \f]
 *
 * is the (expected value) cost-to-go function (also called value function,
 * future value function, future cost function), with \f$ \mathcal{V}_{T}
 * \equiv 0 \f$ and
 *
 * \f[
 *
 *    V_{t}(x_{t-1}, \xi_{t}) =
 *    \min_{x_t \in \mathcal{X}^{n_t}} f_t(x_t) +
 *    \mathcal{V}_{t+1}(x_t, \xi_t)
 * \f]
 *
 * with given \f$ x_{-1} \f$ and (deterministic) \f$ \xi_0\f$. We consider an
 * approximation to the problem associated with stage \f$ t \in \{0, \dots,
 * T-1\} \f$ as the problem
 *
 * \f[
 *    \min_{x_t \in \mathcal{X}^{n_t}} f_t(x_t) +
 *    \mathcal{P}_{t+1}(x_t)
 *    \qquad (1)
 * \f]
 *
 * where \f$ \mathcal{P}_{t+1}(x_t) \f$ is a polyhedral function, i.e., it is
 * a function of the form
 *
 * \f[
 *    \mathcal{P}_{t+1}(x_t) = \max_{i \in \{1,\dots,k_t\}}
 *                                     \{ d_{t,i}^{\top}x_t + e_{t,i} \}
 * \f]
 *
 * with \f$ d_{t,i} \in \mathbb{R}^{n_t} \f$ and \f$ e_{t,i} \in \mathbb{R}
 * \f$ for each \f$ i \in \{1,\dots,k_t\} \f$.
 *
 * An SDDPBlock is then characterized by the following:
 *
 * - It has a time horizon T.
 *
 * - It has T sub-Blocks, each one being a StochasticBlock. The t-th sub-Block
 *   represents an approximation to the problem associated with stage t as
 *   defined in (1).
 *
 * - It has pointers to "T - 1" PolyhedralFunction. The t-th
 *   PolyhedralFunction represents the function \f$ \mathcal{P}_{t+1} \f$ in
 *   (1) and, therefore, must be defined in the t-th sub-Block of this
 *   SDDPBlock or in any of the sub-Blocks of that sub-Block, recursively.
 *
 * - It has a set of scenarios \f$\mathcal{S}\f$. Each scenario in
 *   \f$\mathcal{S}\f$ is represented by a vector of double and spans all the
 *   time horizon T. Each vector is divided into T parts, each one being
 *   associated with a stage of the multistage problem. Let \f$S\f$ denote a
 *   vector representing a scenario in \f$\mathcal{S}\f$. Then, \f$S\f$ is
 *   defined as
 *
 *   \f[
 *     S = ( S_0 , \dots, S_{T-1} )
 *   \f]
 *
 *   where \f$S_t\f$ is a sub-vector of \f$S\f$ with size \f$s_t\f$, for each
 *   \f$t \in \{ 0, \dots, T-1 \}\f$, and is associated with the sub-problem
 *   at stage \f$t\f$, i.e., it provides data for the \f$t\f$-th sub-Block of
 *   this SDDPBlock. We say that \f$S_t\f$ represents the \f$t\f$-th
 *   sub-scenario of the scenario represented by \f$S\f$.
 *
 *   We assume that the sub-scenarios are organized in such a way that related
 *   random data appear in contiguous areas of the sub-scenario. For instance,
 *   suppose that the random data is associated with demand, inflow, and wind
 *   power. In this case, the data related to demand should be a contiguous
 *   sub-vector \f$D_t\f$ of the sub-scenario associated with stage \f$t\f$,
 *   as well as that related to inflow (\f$F_t\f$) and wind power
 *   (\f$W_t\f$). In this example, the sub-scenario \f$S_t\f$ could be
 *   organized as
 *
 *   \f[
 *   S_t = ( D_t , F_t , W_t ).
 *   \f]
 *
 *   We say that this sub-scenario has three groups of related random
 *   data. The order in which the groups of related random data appear in
 *   \f$S_t\f$ is not relevant. We could have, for instance,
 *
 *   \f[
 *   S_t = ( W_t , D_t , F_t ).
 *   \f]
 *
 *   But the sub-scenario associated with stage \f$t\f$ must respect the same
 *   order for each scenario in \f$\mathcal{S}\f$.
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
  * - The description of the sub-Blocks of the SDDPBlock. This is given by the
  *   sub-groups "StochasticBlock" and "StochasticBlock_t", for each t in {0,
  *   ..., "TimeHorizon - 1"}. These sub-groups are optional, but they cannot
  *   be all absent. If "StochasticBlock_t" is not provided for some t in {0,
  *   ..., "TimeHorizon - 1"}, then the "StochasticBlock" group must be
  *   provided and contain a complete description of the t-th sub-Block of
  *   this SDDPBlock. If the "StochasticBlock" group is not provided, then
  *   "StochasticBlock_t" must be provided for each t in {0, ..., "TimeHorizon
  *   - 1"} and contain a complete description of the t-th sub-Block of this
  *   SDDPBlock.
  *
  *   If "StochasticBlock_t" is provided but the description of its inner
  *   Block is not provided, then the "StochasticBlock" group must be provided
  *   and contain the description of an inner Block of a StochasticBlock. In
  *   this case, the description of the inner Block provided in the
  *   "StochasticBlock" group will be used to construct the inner Block of the
  *   StochasticBlock described by the "StochasticBlock_t" group.
  *
  *   If "StochasticBlock_t" is provided but the description of its vector of
  *   DataMapping is not provided, then if the "StochasticBlock" group is
  *   provided and contains a description of a vector of DataMapping, then it
  *   is used to construct the vector of DataMapping of the StochasticBlock
  *   described by the "StochasticBlock_t" group.
  *
  * - The AbstractPath group containing the description of a vector of
  *   AbstractPath as described in the AbstractPath class. The number of
  *   AbstractPath must be equal to either 1 or "TimeHorizon - 1". If the number
  *   of AbstractPath is "TimeHorizon - 1" then the i-th AbstractPath in this
  *   vector must be the path to the PolyhedralFunction associated with the
  *   i-th sub-Block of this SDDPBlock. The i-th AbstractPath is taken with
  *   respect to the inner Block of the i-th sub-Block of this SDDPBlock. If
  *   the number of AbstractPath in this vector is 1, then all paths are
  *   assumed to be equal: for each i in {0, ..., TimeHorizon-1}, the provided
  *   AbstractPath will be the path to the PolyhedralFunction associated with
  *   the i-th sub-Block of this SDDPBlock (taken with respect to the inner
  *   Block of this i-th sub-Block).
  *
  * - The "NumberScenarios" dimension specifying the number of
  *   scenarios.
  *
  * - The "ScenarioSize" dimension containing the size of a single
  *   scenario, which spans all stages.
  *
  * - The "SubScenarioSize" variable, of type netCDF::NcUint64
  *   and indexed over dimension "TimeHorizon". This dimension is
  *   optional. If it is not provided, then all sub-scenarios are assumed to
  *   have the same size, i.e.,
  *
  *     s_t = ScenarioSize / TimeHorizon
  *
  *   for all t in {0, ..., "TimeHorizon - 1"}, and "ScenarioSize" is a multiple
  *   of "TimeHorizon". If this dimension is provided, then SubScenarioSize[t]
  *   is the size of the sub-scenario associated with stage t, i.e., s_t =
  *   SubScenarioSize[t], for each t in {0, ..., TimeHorizon-1}. In the latter
  *   case, the following must hold:
  *
  *   \f[
  *     \text{ScenarioSize} = \sum_{t = 0}^{\text{TimeHorizon} - 1}
  *                           \text{SubScenarioSize}[t].
  *   \f]
  *
  * - The two-dimensional variable "Scenarios" of type netCDF::NcDouble and
  *   indexed over the dimensions "NumberScenarios" and "ScenarioSize",
  *   containing the scenarios. The i-th row of "Scenarios" contains the i-th
  *   scenario, so that Scenarios[i][j] is the j-th component of the i-th
  *   scenario.
  *
  * - The "NumberRandomDataGroups" dimension containing the number of
  *   groups of related random data within each sub-scenario. This dimension
  *   is optional. If it is not provided, then we assume that there is a
  *   single group of related random data. Also, this dimension is meaningful
  *   only if all sub-scenarios have the same size.
  *
  * - The "SizeRandomDataGroups" variable, of type netCDF::Uint64 and indexed
  *   over the "NumberRandomDataGroups" dimension, containing the size of each
  *   group of related random data in each sub-scenario. For each i in {0,
  *   ..., "NumberRandomDataGroups - 1"}, NumberRandomDataGroups[i] is the size
  *   of the i-th group of a sub-scenario. This variable is optional. It is
  *   required only if "NumberRandomDataGroups" is provided and
  *   "NumberRandomDataGroups" > 1.
  *
  * - The "StateSize" variable, of type netCDF::Uint64 and being either a
  *   scalar or a one-dimensional array indexed over "TimeHorizon" dimension,
  *   specifying the sizes of the states at each stage. If this variable is a
  *   scalar, then all states are assumed to have the same size given by
  *   "StateSize". If it is an array then, for each t in {0, ...,
  *   TimeHorizon-1}, StateSize[t] contains the size of the state at stage
  *   t. The state being a vector, its size is the dimension of the space in
  *   which it lies.
  *
  * - The "AdmissibleState" variable, of type netCDF::NcDouble, containing an
  *   admissible state for each stage. An admissible state for a stage is an
  *   state which makes the problem at that stage feasible.  If "StateSize" is
  *   scalar and "AdmissibleState" has dimension "StateSize", then all stages
  *   are assumed to have the same admissible state given by
  *   "AdmissibleState". Otherwise, "AdmissibleState" contains the
  *   concatenation of the states for all stages as follows.
  *
  *   - If "StateSize" is scalar then, for each t in {0, ..., TimeHorizon -
  *     1}, an admissible state for stage t is given by
  *
  *       (AdmissibleState[s_t], ..., AdmissibleState[s_t + StateSize - 1]),
  *
  *     where s_t = t * StateSize. In this case, "AdmissibleState" must have
  *     size "TimeHorizon * StateSize".
  *
  *   - If "StateSize" is a one-dimensional array then, for each t in {0, ...,
  *     "TimeHorizon - 1"}, an admissible state for stage t is given by
  *
  *       (AdmissibleState[s_t], ..., AdmissibleState[s_t + StateSize[t] - 1]),
  *
  *     where s_t = \f$ \sum_{i=0}^{t-1} \f$ StateSize[i]. In this case,
  *     "AdmissibleState" must have size equal to
  *
  *     \f[
  *        \sum_{i=0}^{\text{TimeHorizon} - 1} \text{StateSize}[i].
  *     \f]
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

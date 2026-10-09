/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of SDDPBlock, a class for representing a multistage stochastic
 * programming problem specifically designed to be solved by an SDDP solver.
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

#ifndef __SDDPBlock
 #define __SDDPBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "Solution.h"
#include "Objective.h"
#include "PolyhedralFunction.h"
#include "ScenarioGenerator.h"
#include "ScenarioSimulator.h"
#include "ScenarioSet.h"
#include "StochasticBlock.h"
#include "StOpt/sddp/SimulatorSDDPBase.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SDDPBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// SDDPBlock, representing a multistage stochastic programming problem
/** SDDPBlock derives from Block and represents a multistage stochastic
 * programming problem of the form
 *
 * \f[
 *   \min_{x_0 \in \mathcal{X}_0} f_0(x_0) +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_1 \in \mathcal{X}_1} f_1(x_1) +
 *   \mathbb{E} \left \lbrack \dots +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{X}_{T-1}} f_{T-1}(x_{T-1})
 *   + \mathcal{V}_T( x_{T-1} )
 *   \right\rbrack \right\rbrack\right\rbrack \tag{1}
 * \f]
 *
 * where \f$ T \f$ is the number of stages (the "TimeHorizon" of the
 * SDDPBlock, not to be confused with the number of time instants of
 * the problem of a stage), \f$ \mathcal{X}_t \equiv \mathcal{X}_t( x_{t-1} ,
 * \xi_t ) \subseteq \mathbb{R}^{n_t} \f$ for each \f$ t \in \{ 0 , \dots ,
 * T - 1 \} \f$, \f$ \xi = \{ \xi_t \}_{t \in \{ 1 , \dots , T - 1 \}} \f$
 * is a stochastic process, \f$ x_{-1} \f$ and \f$ \xi_0 \f$ are
 * deterministic, and \f$ \mathcal{V}_T \f$ is the terminal value function
 * discussed below (identically zero in the simplest case). Writing
 * \f$ \xi_{[t]} = ( \xi_0 , \dots , \xi_t ) \f$ for the history of the
 * process up to stage \f$ t \f$, (1) is equivalent to the Bellman recursion
 *
 * \f[
 *   V_t( x_{t-1} , \xi_t ) = \min_{x_t \in \mathcal{X}_t( x_{t-1} , \xi_t )}
 *     f_t( x_t ) + \mathcal{V}_{t+1}( x_t ) \;,
 *   \qquad
 *   \mathcal{V}_{t+1}( x_t ) = \mathbb{E} \left\lbrack
 *     V_{t+1}( x_t , \xi_{t+1} ) \mid \xi_{[t]} \right\rbrack \tag{2}
 * \f]
 *
 * for \f$ t = T - 1 , \dots , 0 \f$, where the second part of (2) defines
 * \f$ \mathcal{V}_{t+1} \f$ for \f$ t < T - 1 \f$ while
 * \f$ \mathcal{V}_T \f$ is given, and \f$ V_0( x_{-1} , \xi_0 ) \f$ is the
 * optimal value of (1). The minimization in the first part of (2) is the
 * problem associated with stage \f$ t \f$, and \f$ \mathcal{V}_{t+1} \f$
 * is its (expected) cost-to-go function, also called value function, future
 * value function or future cost function. We consider an approximation to
 * the problem associated with stage \f$ t \f$ in which
 * \f$ \mathcal{V}_{t+1} \f$ is replaced by a polyhedral function
 * \f$ \mathcal{P}_{t+1} \f$,
 *
 * \f[
 *   \min_{x_t \in \mathcal{X}_t( x_{t-1} , \xi_t )} f_t(x_t) +
 *   \mathcal{P}_{t+1}( z_t ) \;,
 *   \qquad
 *   \mathcal{P}_{t+1}( z_t ) = \max \Big\{ \max_{j \in \mathcal{J}_t}
 *     \{ \alpha_{t,j} + \beta_{t,j}^\top z_t \} \,,\, L_t \Big\}
 *   \tag{3}
 * \f]
 *
 * where \f$ z_t \in \mathbb{R}^{m_t} \f$ is the part of \f$ x_t \f$ on which
 * the future depends (the state at the end of stage \f$ t \f$, see below),
 * each \f$ j \in \mathcal{J}_t \f$ is a row (a cut) of the PolyhedralFunction
 * with constant \f$ \alpha_{t,j} \f$ and coefficients \f$ \beta_{t,j} \f$,
 * and \f$ L_t \f$ is its global lower bound (\f$ -\infty \f$ when no bound is
 * set). For minimization problems, as written here, the PolyhedralFunction is
 * convex, and (3) underestimates the problem of stage \f$ t \f$ as long as
 * each cut and the bound underestimate \f$ \mathcal{V}_{t+1} \f$;
 * maximization problems, with a concave PolyhedralFunction (a minimum of
 * affine functions and an upper bound), are handled with all the inequalities
 * reversed.
 *
 * An SDDPBlock is then characterized by the following:
 *
 * - It has a time horizon T.
 *
 * - It has T sub-Blocks, each of which is a StochasticBlock. The t-th
 *   sub-Block represents an approximation to the problem associated with
 *   stage t as defined in (3). See the note below for the case in which it
 *   may have more than T sub-Blocks.
 *
 * - It has pointers to T PolyhedralFunction (to T times
 *   "NumPolyhedralFunctionsPerSubBlock" of them in general, see
 *   deserialize()). The t-th PolyhedralFunction represents the function
 *   \f$ \mathcal{P}_{t+1} \f$ in (3) and, therefore, must be defined in the
 *   t-th sub-Block of this SDDPBlock or in any of the sub-Blocks of that
 *   sub-Block, recursively; the last one, \f$ \mathcal{P}_T \f$, is the
 *   terminal value function \f$ \mathcal{V}_T \f$ of (1).
 *
 * - It has a set of scenarios \f$\mathcal{S}\f$. Each scenario in
 *   \f$\mathcal{S}\f$ is represented by a vector of double and spans all the
 *   time horizon T. Each vector is divided into T parts, each associated with
 *   a stage of the multistage problem. Let \f$S\f$ denote a vector
 *   representing a scenario in \f$\mathcal{S}\f$. Then, \f$S\f$ is defined as
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
 *
 * \par Stages and state
 *
 * For each stage \f$ t \f$, the sub-Block is a StochasticBlock whose inner
 * Block is a BendersBlock, which the methods that set and read the state
 * (set_state(), get_state(), set_admissible_state()) and the Solvers of this
 * module assume without checking it. The Variables of the BendersBlock are a
 * copy \f$ z_{t-1} \in \mathbb{R}^{m_{t-1}} \f$ of the state that stage
 * \f$ t \f$ receives. Its Objective is a BendersBFunction whose inner Block
 * \f$ B_t \f$ is the model of the stage (e.g., a UCBlock), and which, given
 * \f$ z_{t-1} \f$, writes
 *
 * \f[
 *   [ \, M_t z_{t-1} + m_t \, ]_i \tag{4}
 * \f]
 *
 * into one side (the left, the right or both) of the \f$ i \f$-th of a set of
 * RowConstraint of \f$ B_t \f$ (cf. BendersBFunction::set_mapping()). Its
 * value is the optimal value of \f$ B_t \f$, i.e.,
 * \f$ V_t( z_{t-1} , \xi_t ) \f$ of (2) with \f$ \mathcal{V}_{t+1} \f$
 * replaced by \f$ \mathcal{P}_{t+1} \f$, since the PolyhedralFunction of stage
 * \f$ t \f$ lies inside \f$ B_t \f$. Stage \f$ t \f$ passes on the vector
 * \f$ z_t \f$ of the values of the active Variables of its (first)
 * PolyhedralFunction, which a Solver copies into the BendersBlock of stage
 * \f$ t + 1 \f$ before solving it. Hence, \f$ m_t \f$ ("StateSize"[ t ]) must
 * be both the number of those active Variables and the number of Variables of
 * the BendersBlock of stage \f$ t + 1 \f$, and the two vectors must list the
 * same quantities in the same order (which is not checked). The state
 * \f$ z_{-1} \f$ of the first stage is "InitialState", unless the
 * vdblInitialState parameter of the Solver (SDDPSolver, SDDPGreedySolver)
 * is given.
 *
 * \par Cuts
 *
 * Let \f$ y_i \f$ be the dual value of the \f$ i \f$-th row of (4) in the last
 * solution of \f$ B_t \f$, with the sign convention of RowConstraint: for a
 * minimization, \f$ y_i \f$ is the coefficient of the row in the Lagrangian
 * function, i.e., minus the derivative of the optimal value of \f$ B_t \f$
 * w.r.t. the side of the row that is binding. The coefficients of the
 * linearization of the BendersBFunction at the current state
 * \f$ \bar z_{t-1} \f$ (cf.
 * BendersBFunction::get_linearization_coefficients()) are then
 *
 * \f[
 *   g_t = - \sum_{i \in \mathcal{I}_t} y_i \, M_{t,i}^\top \tag{5}
 * \f]
 *
 * where \f$ M_{t,i} \f$ is the \f$ i \f$-th row of \f$ M_t \f$ and
 * \f$ \mathcal{I}_t \f$ contains the rows whose dual value refers to the side
 * written by (4); a nonzero dual value of the other side of a ranged row does
 * not count. Also, the constant \f$ a_t \f$ is computed either as
 * \f$ \bar v - g_t^\top \bar z_{t-1} \f$, with \f$ \bar v \f$ the bound on the
 * optimal value of \f$ B_t \f$ given by its Solver, or from the dual values of
 * all the rows of \f$ B_t \f$ (cf. the intLinComp parameter of
 * BendersBFunction). If \f$ B_t \f$ is convex (e.g., a Linear Program) and the
 * dual values are optimal, then \f$ g_t \f$ is a subgradient of
 * \f$ V_t( \cdot , \xi_t ) \f$ at \f$ \bar z_{t-1} \f$ and
 *
 * \f[
 *   V_t( z , \xi_t ) \geq a_t + g_t^\top z
 *   \quad \text{for all } z \in \mathbb{R}^{m_{t-1}} \tag{6}
 * \f]
 *
 * holds for \f$ V_t \f$ computed with \f$ \mathcal{P}_{t+1} \f$, and therefore
 * also for \f$ V_t \f$ of (2), since
 * \f$ \mathcal{P}_{t+1} \leq \mathcal{V}_{t+1} \f$. In its backward pass a
 * Solver (see SDDPSolver) computes (6) at a trial state for the sub-scenarios
 * of stage \f$ t \f$. It then combines these linearizations into an estimate
 * \f$ \alpha + \beta^\top z_{t-1} \f$ of a linearization of
 * \f$ \mathcal{V}_t \f$ (an average over the sub-scenarios, or a regression on
 * the random data that depend on the previous stage), and adds it as a new row
 * of \f$ \mathcal{P}_t \f$, i.e., of the PolyhedralFunction of stage
 * \f$ t - 1 \f$ (cf. add_cuts()).
 *
 * The cuts underestimate \f$ \mathcal{V}_t \f$ only if the latter is convex
 * in the state, which holds by backward induction when each stage problem is
 * convex jointly in its state and its decisions (e.g., a linear program in
 * which the state enters only through (4)): \f$ \mathcal{V}_T \f$ is convex
 * (a polyhedral function); if \f$ \mathcal{V}_{t+1} \f$ is convex,
 * \f$ f_t( x_t ) + \mathcal{V}_{t+1}( x_t ) \f$ is convex in \f$ x_t \f$, the
 * feasible set \f$ \{ ( z_{t-1} , x_t ) : x_t \in \mathcal{X}_t \} \f$ is
 * convex, and the minimum over \f$ x_t \f$ of a function jointly convex in
 * \f$ ( z_{t-1} , x_t ) \f$ is convex in \f$ z_{t-1} \f$, and therefore each
 * \f$ V_t( \cdot , \xi_t ) \f$ is convex, and so is their expectation
 * \f$ \mathcal{V}_t \f$. When the random data of stage \f$ t \f$ depend on
 * those of stage \f$ t - 1 \f$ (the time-related random variables of
 * SDDPSolver), \f$ \mathcal{V}_t \f$ depends on these data as well, and so
 * does its model: SDDPSolver estimates the cuts by a regression on them, and
 * before solving a stage it writes in \f$ \mathcal{P}_{t+1} \f$ the cuts of
 * the current realization of these data, and therefore (3) holds with
 * \f$ \alpha_{t,j} \f$ and \f$ \beta_{t,j} \f$ depending on it. The
 * linearization (6) of one sub-scenario at the trial state is a cut of
 * \f$ V_t( \cdot , \xi_t ) \f$ for that sub-scenario only; SDDPSolver can
 * store these random cuts (see store_random_cut() and get_random_cut()),
 * which SDDPGreedySolver uses to give the derivative of the cost of a
 * scenario with respect to its initial state.
 *
 * \par Terminal condition
 *
 * At the last stage, the PolyhedralFunction \f$ \mathcal{P}_T \f$ is data of
 * the problem and stands for \f$ \mathcal{V}_T \f$ in (1); no Solver computes
 * cuts for it, and remove_cuts() leaves it untouched. SDDPSolver::compute()
 * keeps the rows and the bound it has and adds to it the cuts given in the
 * vdblLastStageCuts parameter of SDDPSolver; if it still has neither rows nor
 * a finite bound, it receives the all-zero cut (or the bound 0, if it has no
 * active Variable). Thus, \f$ \mathcal{V}_T \equiv 0 \f$ unless last-stage
 * cuts or a finite bound are given, either in that PolyhedralFunction or
 * (the cuts) through SDDPSolver; with a finite global bound \f$ L \f$ and
 * no rows, \f$ \mathcal{V}_T \equiv L \f$, and with rows (1) also values
 * the state left at the end of the horizon (e.g., the water left in the
 * reservoirs).
 *
 * \par Initial conditions that are not in the state
 *
 * Only the Variables of the BendersBlock of each stage form the state, and
 * only the rows written by (4) depend on it. Any other initial condition of
 * \f$ B_t \f$ (e.g., the initial flow rates of hydro units, the initial
 * commitment, up/down time and power of thermal units, the initial power and
 * storage of batteries) is data of \f$ B_t \f$. A Solver that computes the
 * cuts (SDDPSolver, ParallelSDDPSolver) does not update these data with the
 * solution of stage \f$ t - 1 \f$, which amounts to fixing each such
 * condition to a prescribed value; the cuts are then those of a problem in
 * which the stages are decoupled in these conditions. We remark that this
 * problem approximates the one in which the conditions are part of the state,
 * but it is a restriction of it (and its value an upper bound) only if the
 * final value of each condition in stage \f$ t - 1 \f$ is also required to
 * equal the value prescribed for stage \f$ t \f$, which is not done. A Solver
 * that simulates a policy may instead update these data between stages. In
 * fact, SDDPGreedySolver calls, right before solving each stage, the function
 * given to SDDPGreedySolver::set_callback(), which may copy into \f$ B_t \f$
 * the final values of stage \f$ t - 1 \f$ (say, the final flow rate of each
 * hydro unit as its initial flow rate, the final power of each thermal unit
 * and battery as its initial power, the final storage of each battery as its
 * initial storage). The values so written remain in the data of the stage,
 * unless the function given to SDDPGreedySolver::set_end_callback(), which
 * compute() calls at the end of the simulation, restores the values the
 * stage had before it (as the investmentblock_solver tool does); then a
 * Solver that later computes cuts on the same sub-Blocks starts from the
 * data of the stages, and otherwise from the values the simulation has
 * written. Furthermore, three
 * other treatments of these conditions are not modeled: (i) adding them to
 * the state, (ii) requiring each of them to have, at the end of each stage,
 * the value it has at the beginning, and (iii) relaxing them, i.e., letting
 * each stage choose its initial conditions, which gives a lower
 * approximation.
 *
 * \par Random data
 *
 * At stage \f$ t \f$, the random data are those that the DataMapping of its
 * StochasticBlock write when a sub-scenario \f$ S_t \f$ is set (cf.
 * set_scenario()). They may be any datum of \f$ B_t \f$ or of the
 * BendersBFunction for which a DataMapping can be built, e.g., the constants
 * \f$ m_t \f$ of (4) through BendersBFunction::modify_constants(). Note that
 * (4) is written again on the whole side of each mapped row whenever the
 * state changes; hence, a random datum that appears in that side must be
 * written into \f$ m_t \f$ (a change of it in \f$ B_t \f$ would be
 * overwritten). In a model of an electrical system the random data typically
 * are the demand, the inflows and the renewable generation. Outages of units
 * are not modeled as a random process, as no stage represents the failure of
 * a unit or carries it over to the following stages, although one can write a
 * different availability per scenario where the unit has a setter for it that
 * a DataMapping can call (e.g., ThermalUnitBlock::set_availability()). How
 * the expectation in (2) is estimated from the scenarios, and which random
 * data are allowed where, is described in SDDPSolver.
 *
 * \par Seasonal storage valuation
 *
 * In the seasonal valuation of the storages of an electrical system the
 * horizon is split into consecutive periods, one per stage, and \f$ B_t \f$
 * is a UCBlock over the \f$ \bar T_t \f$ instants
 * \f$ 0 , \dots , \bar T_t - 1 \f$ of period \f$ t \f$. In this setting the
 * seasonal storages are the reservoirs \f$ n \f$ of the HydroUnitBlock of a
 * HydroSystemUnitBlock, and \f$ z_{t-1} \f$ collects their initial volumes
 * \f$ V^0_n \f$. These enter the problem only through the right-hand side
 * \f$ V^0_n + A_{n,0} \f$ of the water balance rows at instant 0 (cf.
 * HydroUnitBlock), with \f$ A_{n,0} \f$ the inflow of instant 0, and (4)
 * writes them there; hence, the inflows of instant 0 belong to \f$ m_t \f$.
 * Conversely, the final volumes \f$ v^{hy}_{n , \bar T_t - 1} \f$ are the
 * active Variables of the PolyhedralFunction of the HydroSystemUnitBlock (the
 * future cost of the water), and \f$ z_t \f$ collects them in the same order
 * (it is the vector \f$ v^{f} \f$ of \ref ucbm_multi). Summing the water
 * balance over the instants of the period shows that the final volumes are an
 * affine function of the initial volumes, the flows and the inflows (see
 * HydroUnitBlock), and hence that the state evolves linearly from one stage
 * to the next; the cost of a stage can therefore be seen either as a function
 * of the flows (with the final volumes given by that affine map) or, as here,
 * as the value of the problem in which the final volumes are variables tied
 * to the flows by the water balances; of course, the two are the same
 * problem. By (5), when \f$ M_t \f$ is the identity, the coefficient of
 * \f$ V^0_n \f$ in \f$ g_t \f$ is minus the dual value of the water balance
 * row of reservoir \f$ n \f$ at instant 0, i.e., the derivative of the cost
 * of the stage w.r.t. the initial volume (usually nonpositive, as more water
 * lowers the cost). Its coefficient in a cut is then the estimate of the
 * expectation of this derivative over the sub-scenarios. A UCBlock has
 * integer variables, while (6) requires a convex \f$ B_t \f$. Thus, the dual
 * values in (5) must be those of a convex relaxation of \f$ B_t \f$, e.g.,
 * its continuous relaxation, or the Lagrangian relaxation of its linking
 * constraints solved by a LagrangianDualSolver, which (for linear objectives)
 * amounts to replacing the feasible set of each unit by its convex hull, and
 * in general its cost by the convex envelope of the cost on that convex hull
 * (the closed convex hull of the epigraph). The cuts are then valid for the
 * value function of that relaxation, and the duals of the water balances,
 * which are rows of a subproblem, are those described in \ref ucbm_multi.
 * Besides the three treatments of the non-state initial conditions listed
 * above, two things are not modeled. A UCBlock has a future value function on
 * the final volumes only, and therefore no other state is valued at the end
 * of a stage (e.g., neither the energy left in batteries nor the energy
 * curtailed by a load curtailment unit); also, the volumes of several
 * reservoirs cannot be aggregated into fewer state variables.
 *
 * \par Several sub-Blocks per stage
 *
 * In the simplest case, an SDDPBlock has \f$ T \f$ sub-Blocks, where the
 * \f$t\f$-th one is a StochasticBlock associated with stage \f$t\f$. However,
 * it may be interesting to have more than one sub-Block associated with each
 * stage in some circumstances. The solution process implemented by
 * SDDPSolver, for instance, involves the solution of a number of subproblems,
 * each of them associated with a stage \f$ t \f$, a particular sub-scenario
 * \f$ S_t^i \f$, and some initial state. There are at least two clear
 * situations under which the presence of multiple sub-Blocks for each stage
 * can be beneficial to SDDPSolver.
 *
 * -# To change the sub-Block as little as possible.
 *
 *    Whenever a subproblem associated with a stage \f$ t \f$ must be solved,
 *    the data of the Block associated with that subproblem must be updated
 *    according to some sub-scenario \f$ S_t^i \f$ and some initial state. In
 *    the case in which the SDDPBlock has \f$ T \f$ sub-Blocks, this means
 *    that its \f$t\f$-th sub-Block must be updated each time a particular
 *    scenario and state is considered. In order to allow reoptimization,
 *    SMS++ is designed to deal with changes in the data of a Block by means
 *    of its Modification mechanism. However, one would expect, in general,
 *    that the less a Block is modified, the faster it can be reoptimized. In
 *    the ideal case, an SDDPBlock would have as many sub-Blocks for each
 *    stage as there are scenarios. In this case, each sub-Block would be
 *    associated with a particular scenario, and the data of each of these
 *    sub-Blocks that depend on the scenarios would be updated only once, in
 *    the beginning. Of course, the sub-Block must still be modified each time
 *    before it is solved, because it also depends on the initial state. But
 *    in general, most of the data in a Block that needs to be updated is
 *    dependent on the scenarios. Therefore, having one sub-Block associated
 *    with each scenario would imply that only the initial state of the
 *    sub-Block must be updated, and the reoptimization could be expected to
 *    be faster (besides avoiding the modification of the scenario of a
 *    sub-Block).
 *
 * -# To allow parallelization in a shared-memory multiprocessing system.
 *
 *    At each iteration, and for each stage, SDDPSolver must solve a set of
 *    subproblems, each one associated with some scenario and some initial
 *    state. If SDDPBlock has only one sub-Block for each stage, then the
 *    process of solving all those subproblems is inevitably sequential, as
 *    solving a particular subproblem requires changing the data of that
 *    sub-Block. The presence of multiple sub-Blocks per stage, however, makes
 *    it possible to parallelize this procedure. Suppose, for instance, that
 *    at each iteration, for each stage, SDDPSolver must solve N subproblems,
 *    each one associated with some scenario and some initial state. Suppose
 *    also that SDDPBlock has B sub-Blocks per stage and a process running
 *    ParallelSDDPSolver::compute() has M threads available. In this case, it
 *    would be possible to allocate min( B , N , M ) subproblems to the
 *    available threads.
 *
 * The decision about the number of Blocks per stage must be well thought out,
 * as it depends, in particular, on the memory resources available. Besides
 * the memory required to store multiple sub-Blocks per stage, one has also to
 * take into account the memory required by the Solver attached to the inner
 * Block of each sub-Block of this SDDPBlock.
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
   delete( block );
  v_Block.clear();
  delete f_scenario_generator;
  f_scenario_generator = nullptr;
 }

/*--------------------------------------------------------------------------*/
 /// loads SDDPBlock out of an istream - not implemented yet

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "SDDPBlock::load: method not implemented yet." ) );
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
  * - The "NumPolyhedralFunctionsPerSubBlock" dimension, containing the number
  *   of PolyhedralFunction that are present in each sub-Block. This dimension
  *   is optional. If it is not provided, then we assume that there is a
  *   single PolyhedralFunction in each sub-Block.
  *
  * - The AbstractPath group containing the description of a vector of
  *   AbstractPath as described in the AbstractPath class. The number of
  *   AbstractPath must be equal to either "NumPolyhedralFunctionsPerSubBlock"
  *   or "NumPolyhedralFunctionsPerSubBlock * TimeHorizon". If the number of
  *   AbstractPath is "NumPolyhedralFunctionsPerSubBlock" then the i-th
  *   PolyhedralFunction of each sub-Block is given by the i-th
  *   AbstractPath. If the number of AbstractPath is
  *   "NumPolyhedralFunctionsPerSubBlock * TimeHorizon" then the i-th
  *   PolyhedralFunction of a sub-Block at stage t is given by the
  *   AbstractPath at position "i + t * NumPolyhedralFunctionsPerSubBlock" for
  *   each t in {0, ..., TimeHorizon - 1} and i in {0, ...,
  *   NumPolyhedralFunctionsPerSubBlock - 1}. An AbstractPath associated with a
  *   stage t is taken with respect to the inner Block of a sub-Block for
  *   stage t of this SDDPBlock.
  *
  * - The "NumSubBlocksPerStage" dimension, containing the number of
  *   sub-Blocks that must be constructed for each stage. This dimension is
  *   optional. If it is not provided, then we assume that there is a single
  *   sub-Block for each stage.
  *
  * - The description of a ScenarioSet, as specified in the comments to
  *   ScenarioSet::deserialize().
  *
  * - The "InitialState" variable, a one-dimensional array of type
  *   netCDF::NcDouble, containing an initial state for the first stage
  *   problem.
  *
  * - The "StateSize" variable, of type netCDF::Uint and being either a scalar
  *   or a one-dimensional array indexed over "TimeHorizon" dimension,
  *   specifying the sizes of the states at each stage. If this variable is a
  *   scalar, then all states are assumed to have the same size given by
  *   "StateSize". If it is an array then, for each t in {0, ...,
  *   TimeHorizon - 1}, StateSize[t] contains the size of the final state at
  *   stage t. The state being a vector, its size is the dimension of the
  *   space in which it lies.
  *
  * - The "AdmissibleState" variable, of type netCDF::NcDouble, containing an
  *   admissible state for each stage. An admissible state for a stage is a
  *   feasible final state for the problem at that stage. If "StateSize" is
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

 void deserialize( const netCDF::NcGroup & group ) override {
  // TimeHorizon

  Index time_horizon;
  ::SMSpp_di_unipi_it::deserialize_dim( group , "TimeHorizon" ,
                                        time_horizon , false );

  // NumSubBlocksPerStage
  // if the dimension is not in the file, the current value (as possibly set
  // beforehand by set_num_sub_blocks_per_stage(), 1 by default) is kept

  ::SMSpp_di_unipi_it::deserialize_dim
   ( group , "NumSubBlocksPerStage" , num_sub_blocks_per_stage );

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
    throw( std::invalid_argument(
     "SDDPBlock::deserialize: The number of AbstractPath to "
     "PolyhedralFunction must be either equal to 1 or equal to "
     "the time horizon." ) );
   else
    throw( std::invalid_argument(
     "SDDPBlock::deserialize: The number of AbstractPath to "
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
      throw( std::invalid_argument(
       "SDDPBlock::deserialize: PolyhedralFunction for stage "
       + std::to_string( t ) + " was not found." ) );
     v_polyhedral_functions.push_back( polyhedral_function );
    }
   }
  }

  // Scenarios
  //
  // Two mutually-exclusive paths are supported:
  //
  // (inline) The scenarios are stored inline at the SDDPBlock group level
  //   (variables "Scenarios", "NumberScenarios", "ScenarioSize" plus the
  //   structural metadata SubScenarioSize / NumberRandomDataGroups /
  //   SizeRandomDataGroups). ScenarioSet::deserialize() handles the lot.
  //
  // (generator) A ScenarioGenerator is provided via a "ScenarioGenerator" sub-
  //   group (with the standard factory schema, i.e. "type" attribute);
  //   the actual scenarios are produced by the generator at runtime, via
  //   init_*_pool() / load_scenarios_from_generator(). In this case the
  //   variables "Scenarios", "NumberScenarios" and "ScenarioSize" are NOT
  //   present at the SDDPBlock level (they would have nothing meaningful
  //   to encode, since the pool is yet to be chosen), but the structural
  //   metadata (SubScenarioSize / NumberRandomDataGroups /
  //   SizeRandomDataGroups) remains at the SDDPBlock level, exactly where
  //   they would be with the inline scenarios. An optional sub-group
  //   "ScenarioGeneratorConfig" can carry a Configuration to be passed to
  //   the generator's set_config().
  //
  // If both forms are present in the same netCDF we throw, to avoid
  // silently picking one over the other. If neither is present, we throw
  // as well (the SDDPBlock has no source of scenarios).

  const auto gen_group = group.getGroup( "ScenarioGenerator" );
  const auto scenarios_var = group.getVar( "Scenarios" );
  const bool has_generator = ! gen_group.isNull();
  const bool has_inline = ! scenarios_var.isNull();

  if( has_generator && has_inline )
   throw( std::logic_error( "SDDPBlock::deserialize: both the "
                            "'Scenarios' variable and the "
                            "'ScenarioGenerator' sub-group are present in "
                            "the netCDF; only one of the two is allowed." ) );

  if( ( ! has_generator ) && ( ! has_inline ) )
   throw( std::logic_error( "SDDPBlock::deserialize: neither the "
                            "'Scenarios' variable nor a 'ScenarioGenerator' "
                            "sub-group are present in the netCDF; no source "
                            "of scenarios available." ) );

  if( has_generator ) {

   // (generator) build the generator via factory, set ourselves as its
   // partner Block, optionally apply ScenarioGeneratorConfig, and
   // populate the structural metadata of scenario_set (the actual data
   // will be filled later by the attached Solver, via
   // load_scenarios_from_generator()).

   f_scenario_generator =
    ::SMSpp_di_unipi_it::ScenarioGenerator::new_ScenarioGenerator(
                                                              gen_group );
   if( ! f_scenario_generator )
    throw( std::logic_error( "SDDPBlock::deserialize: failed to "
                             "deserialize the 'ScenarioGenerator' "
                             "sub-group." ) );

   f_scenario_generator->set_Block( this );

   // optional ScenarioGeneratorConfig: any standard Configuration that
   // the specific :ScenarioGenerator can interpret in its set_config()
   const auto cfg_group = group.getGroup( "ScenarioGeneratorConfig" );
   if( ! cfg_group.isNull() ) {
    auto cfg = ::SMSpp_di_unipi_it::Configuration::new_Configuration(
                                                              cfg_group );
    f_scenario_generator->set_config( cfg );
    delete cfg;
    }

   // read SDDPBlock-level structural metadata. TimeHorizon was already
   // read at the top of this method; the per-stage sub_scenario_size
   // vector is derived differently depending on the generator kind:
   //
   //  - multi-stage generator: derive from the generator itself by
   //    walking a View with descend() and reading get_scenario_size()
   //    at each stage. The 'SubScenarioSize' netCDF attribute is
   //    ignored if present (it would be meaningful only on the
   //    SDDPBlock side if it were authoritative, which it is not here:
   //    the generator owns the per-stage sizes). Requires the generator
   //    to be walkable right after deserialize() (lazy-init contract).
   //
   //  - single-stage generator: read 'SubScenarioSize' from the netCDF
   //    group (optional: if absent, all sub-scenarios share the same
   //    size derived from the generator's scenario_size and time
   //    horizon, as with the inline scenarios).
   std::vector< Index > sub_scenario_size;
   if( auto mgen = dynamic_cast<
       ::SMSpp_di_unipi_it::MultiStageScenarioGenerator * >(
                                                f_scenario_generator ) ) {
    const auto T = mgen->get_stage_number();
    if( T != time_horizon )
     throw( std::logic_error( "SDDPBlock::deserialize: the multi-stage "
                              "generator's stage number (" +
                              std::to_string( T ) + ") does not match "
                              "SDDPBlock's TimeHorizon (" +
                              std::to_string( time_horizon ) + ")." ) );
    sub_scenario_size.reserve( time_horizon );
    auto view = mgen->root_view();
    sub_scenario_size.push_back( view->get_scenario_size() );
    for( Index t = 1 ; t < time_horizon ; ++t ) {
     if( ! view->descend() )
      throw( std::logic_error( "SDDPBlock::deserialize: multi-stage "
                               "generator failed to advance to stage " +
                               std::to_string( t ) + "." ) );
     sub_scenario_size.push_back( view->get_scenario_size() );
     }
    // the View is a position of its own, so the generator is left in the
    // canonical "at first stage" state downstream consumers expect
    }
   else if( ! ::SMSpp_di_unipi_it::deserialize( group , "SubScenarioSize" ,
                                                time_horizon ,
                                                sub_scenario_size ,
                                                true , false ) ) {
    const auto sz = f_scenario_generator->get_scenario_size();
    if( sz % time_horizon != 0 )
     throw( std::logic_error( "SDDPBlock::deserialize: 'SubScenarioSize' "
                              "was not provided with a ScenarioGenerator, "
                              "but the generator's scenario_size (" +
                              std::to_string( sz ) + ") is not a multiple "
                              "of 'TimeHorizon' (" +
                              std::to_string( time_horizon ) + ")." ) );
    sub_scenario_size.assign( time_horizon , sz / time_horizon );
    }

   Index num_random_data_groups = 1;
   std::vector< Index > size_random_data_groups;

   // NumberRandomDataGroups / SizeRandomDataGroups are meaningful only if
   // all sub-scenarios share the same size (as in
   // ScenarioSet::deserialize())
   if( std::adjacent_find( sub_scenario_size.begin() ,
                           sub_scenario_size.end() ,
                           std::not_equal_to<>() ) ==
       sub_scenario_size.end() ) {

    if( ! ::SMSpp_di_unipi_it::deserialize_dim( group ,
                                                "NumberRandomDataGroups" ,
                                                num_random_data_groups ,
                                                true ) )
     num_random_data_groups = 1;

    if( num_random_data_groups > 1 ||
        group.getVar( "SizeRandomDataGroups" ).isNull() == false ) {
     if( ! ::SMSpp_di_unipi_it::deserialize( group ,
                                             "SizeRandomDataGroups" ,
                                             num_random_data_groups ,
                                             size_random_data_groups ,
                                             true , false ) ) {
      size_random_data_groups = { sub_scenario_size.front() };
      if( num_random_data_groups > 1 )
       throw( std::logic_error( "SDDPBlock::deserialize: "
                                "'SizeRandomDataGroups' must be provided "
                                "since 'NumberRandomDataGroups' > 1." ) );
      }
     }
    }

   scenario_set.set_structural_metadata( time_horizon ,
                                         std::move( sub_scenario_size ) ,
                                         num_random_data_groups ,
                                         std::move(
                                          size_random_data_groups ) );
   }
  else {
   // (inline) ScenarioSet handles everything
   scenario_set.deserialize( group );
   }

  // Initial state

  ::SMSpp_di_unipi_it::deserialize( group , "InitialState" ,
                                    initial_state , false );

  // StateSize

  std::vector< Index > state_size;

  ::SMSpp_di_unipi_it::deserialize( group , "StateSize" , time_horizon ,
                                    state_size , false , true );

  bool state_size_is_scalar = ( state_size.size() == 1 );
  if( state_size.size() == 1 )
   state_size.resize( time_horizon , state_size[ 0 ] );
  else if( state_size.size() != time_horizon )
   throw( std::logic_error( "SDDPBlock::deserialize: 'StateSize' must be "
                            "either a scalar or an array with size "
                            "'TimeHorizon'." ) );

  // AdmissibleState

  ::SMSpp_di_unipi_it::deserialize( group , "AdmissibleState" ,
                                    admissible_states , false );

  if( state_size_is_scalar ) {
   if( admissible_states.size() != state_size[ 0 ] &&
       admissible_states.size() != time_horizon * state_size[ 0 ] )
    throw( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                             "array has an invalid size." ) );

   if( admissible_states.size() != time_horizon * state_size[ 0 ] ) {
    std::vector< double > state = admissible_states;
    admissible_states.reserve( time_horizon * state_size[ 0 ] );
    for( Index t = 1 ; t < time_horizon ; ++t )
     admissible_states.insert( admissible_states.cend() ,
                               state.cbegin() , state.cend() );
   }
  }
  else if( admissible_states.size() !=
           std::accumulate( state_size.begin() , state_size.end() ,
                            decltype( state_size )::value_type( 0 ) ) ) {
   throw( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
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

 /// deserialize the random cuts
 /** This function deserializes the random cuts out of the file whose path is
  * given by \p filename. If the path to the file is empty, then no operation
  * is performed. The file must have the following netCDF format:
  *
  * - The netCDF dimension "TimeHorizon" containing the number of stages.
  *
  * - The netCDF dimension "NumberScenarios" containing the number of
  *   scenarios.
  *
  * - The netCDF group "PolyhedralFunction_t_s", for each t in {0, ...,
  *   TimeHorizon - 1} and s in {0, ..., NumberScenarios - 1}, containing the
  *   serialization of the PolyhedralFunction representing the random cuts
  *   associated with stage t and scenario s.
  *
  * @param filename The path to the file containing the netCDF description of
  *        the random cuts. */

 void deserialize_random_cuts( const std::string & filename );

/*--------------------------------------------------------------------------*/

 /// deserialize cuts from the given file and add them to this SDDPBlock
 /** This function deserializes the cuts contained in the file with the given
  * name (path) and adds them to the PolyhedralFunction of every sub-Block of
  * the corresponding stage. If \p filename is empty, then no operation is
  * performed. The file must have one of the two formats described in the
  * comments of serialize_cuts(), which is automatically detected by looking
  * at the leading magic bytes: a netCDF file is dispatched to the netCDF
  * reader, anything else to the CSV one. In both cases the cuts there
  * described are added to (rather than replacing) the current ones of the
  * PolyhedralFunction of each sub-Block of the corresponding stage.
  *
  * @param filename The path to the file containing the description of
  *        the cuts. */

 void deserialize_cuts( const std::string & filename );

/*--------------------------------------------------------------------------*/

 /// returns true if the file has the netCDF (classic or HDF5) magic bytes
 /** Helper telling whether the file with the given name (path) is a netCDF
  * one, so that the format of a cuts file can be automatically detected;
  * see deserialize_cuts(). */

 static bool is_netCDF_file( const std::string & filename );

/*--------------------------------------------------------------------------*/

 /// sets the number of sub-Blocks for each stage
 /** This function sets the number of sub-Blocks that must be constructed at
  * each stage. If this function is invoked after the sub-Blocks of this
  * SDDPBlocks have been constructed, it has no effect. In particular, it has
  * no effect if it is invoked after deserialize() is invoked.
  *
  * @param n The number of sub-Blocks that must be constructed at each stage.
  */
 void set_num_sub_blocks_per_stage( Index n ) {
  if( v_Block.empty() )
   num_sub_blocks_per_stage = n;
 }

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR Saving THE DATA OF THE SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Saving the data of the SDDPBlock
 *  @{ */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// serialize an SDDPBlock into a netCDF::NcGroup
 /** Serialize an SDDPBlock into a netCDF::NcGroup with the format
  * explained in the comments of the deserialize() function.
  *
  * @param group The NcGroup in which this SDDPBlock will be serialized. */

 void serialize( netCDF::NcGroup & group ) const override;

/*--------------------------------------------------------------------------*/

 /// serialize the random cuts
 /** This function serializes the random cuts in the file with the given name
  * (path). If \p filename is empty, then no operation is performed. The file
  * will have the following netCDF format:
  *
  * - The dimension "TimeHorizon" containing the number of stages.
  *
  * - The dimension "NumberScenarios" containing the number of scenarios.
  *
  * - The group "PolyhedralFunction_t_s", for each t in {0, ..., TimeHorizon -
  *   1} and s in {0, ..., NumberScenarios - 1}, containing the serialization
  *   of the PolyhedralFunction representing the random cuts associated with
  *   stage t and scenario s. Each individual group is optional. If the group
  *   "PolyhedralFunction_t_s" is not provided, then the PolyhedralFunction
  *   associated with stage t and scenario s will not be loaded (which means
  *   it will have no cuts).
  *
  * @param filename The name of the file in which the random cuts will be
  *        serialized. */

 void serialize_random_cuts( const std::string & filename ) const;

/*--------------------------------------------------------------------------*/

 /// serialize the cuts
 /** This function serializes the cuts of this SDDPBlock, i.e., the rows of
  * the PolyhedralFunction of (the first sub-Block of) each stage, in the
  * file with the given name (path). If \p filename is empty, then no
  * operation is performed. The format is selected by the extension of
  * \p filename: with ".nc4" or ".nc" the file will have the following
  * netCDF format:
  *
  * - The dimension "TimeHorizon" containing the number of stages.
  *
  * - The group "PolyhedralFunction_t", for each t in {0, ..., TimeHorizon -
  *   1}, containing the serialization of the PolyhedralFunction associated
  *   with stage t.
  *
  * With any other extension the file will be a CSV with a header line
  *
  *     Timestep,a_0,...,a_{n-1},b
  *
  * followed by one line per cut of the form "t,A[i][0],...,A[i][n-1],b[i]",
  * where t is the stage the cut belongs to; this is the historical format
  * of this module, kept so that existing consumers keep working unchanged.
  *
  * @param filename The name of the file in which the cuts will be
  *        serialized. */

 void serialize_cuts( const std::string & filename ) const;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the SDDPBlock
    @{ */

 /// returns the time horizon
 /** This function returns the time horizon associated with this SDDPBlock.
  *
  * @return The time horizon.
  */
 virtual Index get_time_horizon() const {
  return( v_Block.size() / num_sub_blocks_per_stage );
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of PolyhedralFunction
 /** This function returns the vector of PolyhedralFunction associated with
  * this SDDPBlock.
  *
  * @return The vector of PolyhedralFunction.
  */
 const std::vector< PolyhedralFunction * > &
 get_polyhedral_functions() const {
  return( v_polyhedral_functions );
 }

/*--------------------------------------------------------------------------*/

 /// returns a PolyhedralFunction
 /** This function returns a pointer to the i-th PolyhedralFunction of a
  * sub-Block of the given \p stage.
  *
  * @param stage A number between 0 and get_time_horizon() - 1.
  *
  * @param i If there are more than one PolyhedralFunction per sub-Block, this
  *        parameter informs the index of the desired PolyhedralFunction at
  *        the given \p stage.
  *
  * @param sub_block_index The index of the sub-Block, which must be an
  *        integer between 0 and get_num_sub_blocks_per_stage() - 1.
  *
  * @return A pointer to the i-th PolyhedralFunction of the given \p stage. */

 PolyhedralFunction * get_polyhedral_function
 ( Index stage , Index i = 0 , Index sub_block_index = 0 ) const {

  if( ! num_polyhedral_per_sub_block )
   // Well, this is a funny SDDPBlock that has no PolyhedralFunction.
   return( nullptr );

  assert( stage < get_time_horizon() );
  assert( sub_block_index < num_sub_blocks_per_stage );
  assert( i < num_polyhedral_per_sub_block );

  const auto index = ( stage * num_sub_blocks_per_stage + sub_block_index )
   * num_polyhedral_per_sub_block + i;

  return( v_polyhedral_functions[ index ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the number of PolyhedralFunction in each sub-Block
 /** This function returns the number of PolyhedralFunction present in each
  * sub-Block.
  *
  * @return The number of PolyhedralFunction in each sub-Block.
  */
 Index get_num_polyhedral_function_per_sub_block() const {
  return( num_polyhedral_per_sub_block );
 }

/*--------------------------------------------------------------------------*/

 /// returns the number of sub-Blocks for each stage
 /** This function returns the number of sub-Blocks for each stage.
  *
  * @return The number of sub-Blocks for each stage.
  */
 Index get_num_sub_blocks_per_stage() const {
  return( num_sub_blocks_per_stage );
 }

/*--------------------------------------------------------------------------*/

 /// returns a sub-Block of this SDDPBlock
 /** This function returns the sub-Block of index \p sub_block_index at the
  * given \p stage of this SDDPBlock. The given \p stage must be an integer
  * between 0 and get_time_horizon() - 1 and the index of the sub-Block must
  * be an integer between 0 and get_num_sub_blocks_per_stage() 0 - 1. If any
  * of them is an invalid index, an exception is thrown.
  *
  * @param stage The stage of the desired sub-Block. It must be a number
  *        between 0 and get_time_horizon() - 1.
  *
  * @param sub_block_index The index of the desired sub-Block at the given \p
  *        stage. It must be a number between 0 and
  *        get_num_sub_blocks_per_stage() - 1.
  *
  * @return The sub-Block of this SDDPBlock associated with the given \p stage
  *         and having index \p sub_block_index.
  */
 virtual StochasticBlock * get_sub_Block( Index stage ,
                                          Index sub_block_index = 0 ) const;

/*--------------------------------------------------------------------------*/

 /// returns the admissible state associated with the given \p stage
 /** This function returns an iterator to the vector containing the admissible
  * state for the given \p stage.
  *
  * @param stage A stage between 0 and get_time_horizon() - 1.
  *
  * @return An iterator to the vector containing the admissible state for the
  *         given \p stage. */

 std::vector< double >::const_iterator get_admissible_state( Index stage ) const {
  assert( stage < get_time_horizon() );
  return( std::next( admissible_states.cbegin() ,
                     admissible_state_begin[ stage ] ) );
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of the admissible state of the given \p stage
 /** This function returns the size of the admissible state associated with
  * the given \p stage.
  *
  * @param stage A stage between 0 and get_time_horizon() - 1.
  *
  * @return The size of the admissible state for the given \p stage. */

 Index get_admissible_state_size( Index stage ) const {
  assert( stage < get_time_horizon() );
  if( stage == get_time_horizon() - 1 )
   return( admissible_states.size() - admissible_state_begin[ stage ] );
  else
   return( admissible_state_begin[ stage + 1 ] -
           admissible_state_begin[ stage ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the set of scenarios
 /** This function returns the set of scenarios. */
 const ScenarioSet & get_scenario_set() const {
  return( scenario_set );
 }

/*--------------------------------------------------------------------------*/

 /// returns the ScenarioGenerator attached to this SDDPBlock (or nullptr)
 /** Returns a non-owning pointer to the ScenarioGenerator attached to this
  * SDDPBlock, or nullptr if no generator has been provided (in which case
  * the scenarios are entirely owned by get_scenario_set()).
  *
  * If a ScenarioGenerator is present, the typical workflow is for the
  * attached Solver to decide which kind of pool it wants (representative
  * for SDDPSolver, random for SDDPGreedySolver), call the appropriate
  * init_*_pool() on the generator with the desired size, and then call
  * prepare_generator_pool() on this SDDPBlock to snapshot the current
  * pool into the internal cache, after which the data-access helpers
  * #size() / #sub_scenario_begin() / #sub_scenario_end() (and therefore
  * #set_scenario() and ScenarioSimulator) will route reads through the
  * cache instead of through the ScenarioSet storage.
  *
  * With a ScenarioGenerator that is not a MultiStageScenarioGenerator, the
  * scenario produced by get_current_scenario() spans all the stages, and
  * the per-stage decomposition is taken from the netCDF SubScenarioSize
  * variable read by deserialize(), as for the scenarios stored in the
  * ScenarioSet; a MultiStageScenarioGenerator gives instead the scenarios
  * of each stage (see prepare_multi_stage_generator_pool()). */

 ScenarioGenerator * get_scenario_generator() const {
  return( f_scenario_generator );
 }

/*--------------------------------------------------------------------------*/

 /// snapshot the generator's current pool into the internal cache
 /** This method walks the pool currently produced by the ScenarioGenerator
  * attached to this SDDPBlock (if any), copying each scenario into the
  * internal #f_generator_pool_cache buffer. From this moment on, the
  * SDDPBlock-side data-access helpers (#size(), #sub_scenario_begin(),
  * #sub_scenario_end()) route reads through this cache instead of
  * through #scenario_set (whose .scenarios storage is left untouched and
  * is in fact empty when the scenarios come from a generator). The number
  * of scenarios is also mirrored into scenario_set via
  * ScenarioSet::set_num_scenarios(), so that a caller that reads
  * get_scenario_set().size() observes the right count.
  *
  * The generator's pool must already have been initialized (typically by
  * the caller via init_random_pool() or init_representative_pool());
  * otherwise an exception is thrown.
  *
  * Each scenario produced by the generator is expected to be a flat
  * vector spanning all stages, of length get_scenario_size() matching the
  * SubScenarioSize structural metadata read from the netCDF group at
  * deserialize time (see deserialize()). The per-scenario probabilities
  * returned by the generator must all be equal to 1 / pool_size (uniform);
  * non-uniform probabilities make the method throw, since the averaging of
  * the cuts of SDDPSolver uses uniform weights.
  *
  * After this method returns, the generator's pool iteration is left at
  * the beginning (reset_pool() is called at the end) so that subsequent
  * external uses of the generator are unaffected.
  *
  * @throws std::logic_error if no generator is attached, or if the
  *         generator's pool is not initialized, or if probabilities are
  *         non-uniform. */

 void prepare_generator_pool();

/*--------------------------------------------------------------------------*/

 /// snapshot a MultiStageScenarioGenerator's per-stage pools into the cache
 /** Analogue of #prepare_generator_pool() for a
  * MultiStageScenarioGenerator. Walks the attached generator and
  * fills #f_multi_stage_pool_cache, a per-stage table of scenario vectors,
  * where each row cache[t][k] is a flat scenario of length
  * scenario_set.get_sub_scenario_size(t).
  *
  * The walk assumes stage independence, i.e., that the pool of a stage is
  * the same whatever the history of the process before it: for each stage
  * t the method navigates to one realization of the history and enumerates
  * next_scenario() until exhausted. A generator whose pools depend on the
  * history (a tree of scenarios) is not handled.
  *
  * As in #prepare_generator_pool(), the per-stage probabilities must be
  * uniform (otherwise weighted-average cuts in the backward pass would
  * be required, which SDDPSolver does not compute). The pool sizes per stage
  * may differ.
  * After this method returns, the generator's pool iteration is reset
  * (reset_pool()) so that subsequent external uses are unaffected.
  *
  * @throws std::logic_error if no MultiStageScenarioGenerator is
  *         attached, or if the generator's pool is not initialized, or
  *         if probabilities are non-uniform, or if the generator reports
  *         a stage_number() incompatible with the SDDPBlock's
  *         get_time_horizon(). */

 void prepare_multi_stage_generator_pool();

/*--------------------------------------------------------------------------*/

 /// returns the number of scenarios of this SDDPBlock at stage \p stage
 /** Dispatches between three possible scenario sources, in order of
  * precedence:
  *
  *  1. The MultiStageScenarioGenerator cache (#f_multi_stage_pool_cache,
  *     populated by #prepare_multi_stage_generator_pool()):
  *     returns the per-stage pool size cache[\p stage].size().
  *
  *  2. The single-stage ScenarioGenerator cache (#f_generator_pool_cache,
  *     populated by #prepare_generator_pool()): returns
  *     cache.size(), the same value for every stage (the scenario spans
  *     the full horizon).
  *
  *  3. The ScenarioSet storage: returns scenario_set.size(), the
  *     same value for every stage.
  *
  * The parameter \p stage is only consulted in case 1.
  *
  * @param stage A stage in [0, get_time_horizon()); defaults to 0, which
  *        is the right value in cases 2 and 3 (uniform across stages). */

 Index size( Index stage = 0 ) const {
  if( ! f_multi_stage_pool_cache.empty() ) {
   assert( stage < f_multi_stage_pool_cache.size() );
   return( f_multi_stage_pool_cache[ stage ].size() );
   }
  if( f_scenario_generator )
   return( f_generator_pool_cache.size() );
  return( scenario_set.size() );
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of each random data group
 /** Forwards to ScenarioSet::get_size_random_data_groups(), as this is
  * structural metadata that lives in #scenario_set in all paths
  * (inline scenarios, single-stage generator, and multi-stage generator). */

 const std::vector< Index > & get_size_random_data_groups() const {
  return( scenario_set.get_size_random_data_groups() );
 }

/*--------------------------------------------------------------------------*/

 /// returns an iterator to the first element of the sub-scenario (i, t)
 /** Dispatches between three possible scenario sources, in order of
  * precedence (see #size()):
  *
  *  1. The MultiStageScenarioGenerator cache: returns
  *     cache[stage][scenario_id].cbegin(). Each row of cache[stage] is
  *     a flat vector of length scenario_set.get_sub_scenario_size(stage).
  *
  *  2. The single-stage ScenarioGenerator cache: returns an iterator
  *     into cache[scenario_id], advanced to the offset of stage \p
  *     stage as dictated by the structural metadata in #scenario_set.
  *
  *  3. The ScenarioSet storage: delegates to
  *     ScenarioSet::sub_scenario_begin().
  *
  * @param scenario_id The index of a scenario, which must be in
  *        [0, size(stage)).
  *
  * @param stage A stage in [0, get_time_horizon()). */

 std::vector< double >::const_iterator
 sub_scenario_begin( Index scenario_id , Index stage ) const {
  if( ! f_multi_stage_pool_cache.empty() ) {
   assert( stage < f_multi_stage_pool_cache.size() );
   if( scenario_id >= f_multi_stage_pool_cache[ stage ].size() )
    throw( std::invalid_argument
           ( "SDDPBlock::sub_scenario_begin: invalid scenario index "
             + std::to_string( scenario_id ) + " at stage "
             + std::to_string( stage ) + "." ) );
   return( f_multi_stage_pool_cache[ stage ][ scenario_id ].cbegin() );
   }
  if( f_scenario_generator ) {
   if( scenario_id >= f_generator_pool_cache.size() )
    throw( std::invalid_argument
           ( "SDDPBlock::sub_scenario_begin: invalid scenario index "
             + std::to_string( scenario_id ) + "." ) );
   assert( stage < scenario_set.get_time_horizon() );
   return( std::next( f_generator_pool_cache[ scenario_id ].cbegin() ,
                      scenario_set.sub_scenario_begin_offset( stage ) ) );
   }
  return( scenario_set.sub_scenario_begin( scenario_id , stage ) );
 }

/*--------------------------------------------------------------------------*/

 /// returns an iterator past the last element of sub-scenario (i,t)

 std::vector< double >::const_iterator
 sub_scenario_end( Index scenario_id , Index stage ) const {
  if( ! f_multi_stage_pool_cache.empty() ) {
   assert( stage < f_multi_stage_pool_cache.size() );
   if( scenario_id >= f_multi_stage_pool_cache[ stage ].size() )
    throw( std::invalid_argument
           ( "SDDPBlock::sub_scenario_end: invalid scenario index "
             + std::to_string( scenario_id ) + " at stage "
             + std::to_string( stage ) + "." ) );
   return( f_multi_stage_pool_cache[ stage ][ scenario_id ].cend() );
   }
  if( f_scenario_generator ) {
   if( scenario_id >= f_generator_pool_cache.size() )
    throw( std::invalid_argument
           ( "SDDPBlock::sub_scenario_end: invalid scenario index "
             + std::to_string( scenario_id ) + "." ) );
   assert( stage < scenario_set.get_time_horizon() );
   return( std::next( f_generator_pool_cache[ scenario_id ].cbegin() ,
                      scenario_set.sub_scenario_begin_offset( stage + 1 ) ) );
   }
  return( scenario_set.sub_scenario_end( scenario_id , stage ) );
 }

/*--------------------------------------------------------------------------*/

 /// returns the initial state for the first stage problem
 /** This function returns the initial state for the first stage problem. */
 const std::vector< double > & get_initial_state() const {
  return( initial_state );
 }

/*--------------------------------------------------------------------------*/

 /// returns a random cut
 /** This function returns the PolyhedralFunction representing the random cut
  * associated with the given \p stage and the scenario whose index is \p
  * scenario_index. The \p stage argument must be between 0 and
  * get_time_horizon() - 1 while \p scenario_index must be between 0 and
  * get_scenario_set().size() - 1. */

 PolyhedralFunction & get_random_cut( Index stage , Index scenario_index ) {
  if( stage >= random_cuts.size() )
   throw( std::invalid_argument( "SDDPBlock::get_random_cut: no random cut "
                                 "for stage " + std::to_string( stage ) ) );

  if( scenario_index >= random_cuts[ stage ].size() )
   throw( std::invalid_argument
          ( "SDDPBlock::get_random_cut: no random cut for scenario index " +
            std::to_string( scenario_index ) + "." ) );

  return( random_cuts[ stage ][ scenario_index ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the sense of the Objective of the SDDPBlock
 /** This function returns the sense of the Objective of the SDDPBlock, which
  * is defined to be the sense of the Objective of its first inner Block. If
  * this SDDPBlock has no inner Block, this function returns
  * Objective::eUndef.
  *
  * @return the sense of the Objective of the first inner Block of this
  *         SDDPBlock if there is one. Otherwise, it returns
  *         Objective::eUndef. */

 int get_objective_sense() const override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a SDDPBlockSolution representing the current solution
 /** Returns a SDDPBlockSolution representing the current solution status of
  * this SDDPBlock, i.e., the Solution of the inner Block of the
  * BendersBFunction at each stage (the first sub-Block of each stage if
  * there are more than one); see SDDPBlockSolution for details. This is
  * meant to be used after a simulation (see SDDPGreedySolver) has solved
  * the stage sub-problems for some scenario and the solutions have been
  * retrieved (get_var_solution()), so that one SDDPBlockSolution describes
  * the complete trajectory of one scenario.
  *
  * What the Solution contains is "configured" by \p solc, which can be:
  *
  * - a SimpleConfiguration< int >, whose f_value is bit-wise coded as
  *   follows: bit 0 tells that the Solution of the inner Block of each
  *   stage is saved; the remaining bits (f_value >> 1), if nonzero, are
  *   passed as a SimpleConfiguration< int > to get_Solution() of the inner
  *   Block (say, a UCBlock, see UCBlock::get_Solution() for the meaning of
  *   the bits);
  *
  * - a SimpleConfiguration< std::pair< int , Configuration * > >, in which
  *   case f_value.first is used as the bit 0 above and f_value.second is
  *   passed (cloned) to get_Solution() of the inner Block of each stage.
  *
  * If \p solc is nullptr, the f_solution_Configuration of the BlockConfig
  * is used if present, otherwise everything (bit 0 == 1, no inner
  * Configuration) is saved. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods describing the behavior of an SDDPBlock
 * @{ */

 /// add cuts to a sub-Block at the given stage
 /** This function adds cuts to the sub-Block with index \p sub_block_index at
  * the given \p stage. The parameters must satisfy the following
  * requirements:
  *
  * 1. \p A must be a matrix with as many columns as there are cuts to be
  *    added and the number of rows must be equal to the number of Variable
  *    defined in the BendersBlock associated with stage \p stage.
  *
  * 2. \p b must be a vector whose size is equal to the number of
  *    rows of \p A. The cuts are given by Ax + b.
  *
  * 3. \p stage must be an integer between 0 and get_time_horizon() - 1.
  *
  * 4. \p sub_block_index must be an integer between 0 and
  *    get_num_sub_blocks_per_stage() - 1.
  *
  * @param A A matrix containing the coefficients of the cuts to be added.
  *
  * @param b A vector containing the constants of the cuts to be added.
  *
  * @param stage The stage whose cuts should be updated.
  *
  * @param sub_block_index The index of the sub-Block whose cuts will be
  *        updated.
  *
  * @param remove_current_cuts If true, all cuts currently part of sub-Block
  *        with index \p sub_block_index at the given \p stage are removed
  *        before the given cuts are added. If false, the current cuts are
  *        kept. */

 void add_cuts( PolyhedralFunction::MultiVector && A ,
                PolyhedralFunction::RealVector && b , Index stage ,
                Index sub_block_index ,
                Index number_cuts_to_keep = Inf< Index >() );

/*--------------------------------------------------------------------------*/
 /// add cuts to all sub-Blocks at the given stage
 /** This function adds cuts to all subs-Blocks at the given \p stage. The
  * parameters must satisfy the following requirements:
  *
  * 1. \p A must be a matrix with as many columns as there are cuts to be
  *    added and the number of rows must be equal to the number of Variable
  *    defined in the BendersBlock associated with stage \p stage.
  *
  * 2. \p b must be a vector whose size is equal to the number of
  *    rows of \p A. The cuts are given by Ax + b.
  *
  * 3. \p stage must be an integer between 0 and get_time_horizon() - 1.
  *
  * @param A A matrix containing the coefficients of the cuts to be added.
  *
  * @param b A vector containing the constants of the cuts to be added.
  *
  * @param stage The stage whose cuts should be updated.
  *
  * @param remove_current_cuts If true, all cuts currently part of sub-Block
  *        with index \p sub_block_index at the given \p stage are removed
  *        before the given cuts are added. If false, the current cuts are
  *        kept. */

 void add_cuts( PolyhedralFunction::MultiVector && A ,
                PolyhedralFunction::RealVector && b , Index stage ,
                Index number_cuts_to_keep = Inf< Index >() ) {
  for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i ) {
   auto A_ = A;
   auto b_ = b;
   add_cuts( std::move( A_ ) , std::move( b_ ) , stage , i ,
             number_cuts_to_keep );
  }
 }

/*--------------------------------------------------------------------------*/
  /// removes all cuts from each PolyhedralFunction
  /** This function removes all cuts from each PolyhedralFunction, except
   * the ones associated with the last stage: those are data of the problem
   * (the final cuts), only set once at the beginning of compute(). */

 void remove_cuts() {
  for( Index stage = 0 ; stage < get_time_horizon() - 1 ; ++stage )
   for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i )
    for( Index sub_block_index = 0 ;
         sub_block_index < num_sub_blocks_per_stage ; ++sub_block_index )
     get_polyhedral_function( stage , i , sub_block_index )->delete_rows();
 }

/*--------------------------------------------------------------------------*/
 /// store the given random cut
 /** This function store the random cut given by \p coefficients and \p alpha,
  * which must be associated with the given \p stage and with the scenario
  * whose index is \p scenario_index.
  *
  * @param coefficients The coefficients of the cut.
  *
  * @param alpha The constant of the cut.
  *
  * @param stage The stage (a number between 0 and get_time_horizon() - 1)
  *        associated with the given cut.
  *
  * @param scenario_index The index (a number between 0 and
  *        get_scenario_set().size() - 1) of the scenario associated with the
  *        given cut. */

 void store_random_cut( std::vector< double > && coefficients , double alpha ,
                        Index stage , Index scenario_index );

/*--------------------------------------------------------------------------*/
 /// returns the number of cuts currently present at the given \p stage
 /** This function returns the number of cuts currently present in the i-th
  * PolyhedralFunction of the sub-Block with index \p sub_block_index at the
  * given \p stage. The parameter \p i is the index of the PolyhedralFunction
  * to which the cuts should be added (its default value is 0).
  *
  * @param stage The stage whose cuts should be updated.
  *
  * @param sub_block_index The index of a sub-Block at the given \p stage.
  *
  * @param i The index of the PolyhedralFunction in the indicated sub-Block.
  *
  * @return The number of cuts currently present in the i-th
  *         PolyhedralFunction of the sub-Block with index \p sub_block_index
  *         at the given \p stage. */

 Index get_number_cuts( Index stage , Index sub_block_index , Index i = 0 )
  const {

  if( stage >= get_time_horizon() )
   throw( std::invalid_argument( "SDDPBlock::get_num__cuts: invalid stage "
                                 "index: " + std::to_string( stage ) ) );

  return( get_polyhedral_function( stage , i , sub_block_index )->get_nrows() );
 }

/*--------------------------------------------------------------------------*/

 /// returns the current future cost of the given sub-Block at the given stage
 double get_future_cost( Index stage , Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem
  * associated with the sub-Block with index \p sub_block_index at the given
  * \p stage. The size of the \p values array parameter must be equal to the
  * number N of state Variable of the problem at the given \p stage, so that
  * the value of the i-th state Variable will be values( i ), for each i in
  * {0, ..., N-1}.
  *
  * @param values The Eigen::ArrayXd containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set. This must be an
  *              integer between 0 and get_time_horizon() - 1.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_state( const Eigen::ArrayXd & values , Index stage ,
                 Index sub_block_index );

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of all sub-Blocks at the given stage
 /** This function sets the values of the state Variable of all sub-Blocks at
  * the given \p stage. The size of the \p values array parameter must be
  * equal to the number N of state Variable of the problem at the given \p
  * stage, so that the value of the i-th state Variable will be values( i ),
  * for each i in {0, ..., N-1}.
  *
  * @param values The Eigen::ArrayXd containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set. This must be an
  *              integer between 0 and get_time_horizon() - 1. */

 void set_state( const Eigen::ArrayXd & values , Index stage ) {
  for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i )
   set_state( values , stage , i );
 }

/*--------------------------------------------------------------------------*/

 /// returns the values of the state Variable of the problem at the given stage
 /** This function returns the current values of the state Variable of the
  * problem associated with the sub-Block with index \p sub_block_index at the
  * given \p stage.
  *
  * @param stage The stage whose state Variable values are desired.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1.
  *
  * @return The current values of the state Variable of the problem at the
  *         given \p stage. */

 std::vector< double > get_state( Index stage ,
                                  Index sub_block_index = 0 ) const;

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem
  * associated with the sub-Block with index \p sub_block_index at the given
  * \p stage. The size of the \p values array parameter must be equal to the
  * number N of state Variable of the problem at the given \p stage, so that
  * the value of the i-th state Variable will be values[ i ], for each i in
  * {0, ..., N-1}.
  *
  * @param values The vector containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_state( const std::vector< double > & values , Index stage ,
                 Index sub_block_index );

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of all sub-Blocks at the given stage
 /** This function sets the values of the state Variable of all sub-Blocks at
  * the given \p stage. The size of the \p values array parameter must be
  * equal to the number N of state Variable of the problem at the given \p
  * stage, so that the value of the i-th state Variable will be values[ i ],
  * for each i in {0, ..., N-1}.
  *
  * @param values The vector containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set. */

 void set_state( const std::vector< double > & values , Index stage ) {
  for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i )
   set_state( values , stage , i );
 }

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem
  * associated with the sub-Block with index \p sub_block_index at the given
  * \p stage, according to the admissible state of this SDDPBlock.
  *
  * @param stage The stage whose state must be set.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_admissible_state( Index stage , Index sub_block_index = 0 );

/*--------------------------------------------------------------------------*/

 /// updates the sub-Block at the given stage for the given scenario
 /** This function updates the sub-Block whose index is \p sub_block_index at
  * the given \p stage for the given \p scenario.
  *
  * @param scenario_id The id of the scenario that must be set.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_scenario( Index scenario_id , Index stage ,
                    Index sub_block_index = 0 ) {
  auto begin = sub_scenario_begin( scenario_id , stage );

  try {
   get_sub_Block( stage , sub_block_index )->set_data( begin );
  }
  catch( const std::exception & e ) {
   std::cout << "SDDPBlock::set_scenario: exception while setting scenario "
             << scenario_id << " of stage " << stage << ".\n"
             << e.what() << std::endl;
   std::exit( EXIT_FAILURE );
  }
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/
 /// Pointers to the PolyhedralFunction of each sub-Block
 /** This vector stores the pointers to the PolyhedralFunction of each
  * sub-Block of this SDDPBlock. The pointer to the i-th PolyhedralFunction of
  * the j-th sub-Block of stage t is located at position
  *
  *   ( t * num_sub_blocks_per_stage + j ) * num_polyhedral_per_sub_block + i
  */
 std::vector< PolyhedralFunction * > v_polyhedral_functions;

 /// Number of PolyhedralFunctions for each sub-Block
 Index num_polyhedral_per_sub_block = 1;

 /// Number of sub-Blocks for each stage
 Index num_sub_blocks_per_stage = 1;

 /// Simulator for the forward step of the SDDP method
 std::shared_ptr< ScenarioSimulator > simulator_forward;

 /// Simulator for the backward step of the SDDP method
 std::shared_ptr< ScenarioSimulator > simulator_backward;

 /// The set of scenarios
 ScenarioSet scenario_set;

 /// The (optional) ScenarioGenerator owned by this SDDPBlock
 /** Owned pointer (deleted in the destructor). Non-null if and only if the
  * SDDPBlock was deserialized with a "ScenarioGenerator" sub-group; in that
  * case, scenario_set retains only the structural metadata and the actual
  * scenario data is snapshotted into #f_generator_pool_cache by
  * prepare_generator_pool(). Null when the scenarios are stored inline,
  * where they live entirely inside scenario_set. */

 ScenarioGenerator * f_scenario_generator = nullptr;

 /// Per-scenario snapshot of the ScenarioGenerator's current pool
 /** Filled by prepare_generator_pool() from the ScenarioGenerator's pool
  * after the attached Solver has called init_*_pool() on the generator.
  * Each row is a flat scenario of length scenario_set.get_scenario_size(),
  * spanning all stages (the per-stage decomposition is taken from
  * #scenario_set's structural metadata). Empty when the scenarios are
  * stored inline, scenario_set.scenarios being their home instead.
  *
  * Mutually exclusive with #f_multi_stage_pool_cache: at most one of the
  * two is non-empty after a successful prepare_* call.
  *
  * The snapshot is necessary because the ScenarioGenerator API is
  * forward-iterator-only (next_scenario() / reset_pool() / get_current_
  * scenario()) and SDDP needs O(1) random access to any scenario at any
  * time during the backward/forward sweeps. */

 std::vector< std::vector< double > > f_generator_pool_cache;

 /// Per-stage snapshot of a MultiStageScenarioGenerator's pools
 /** Filled by prepare_multi_stage_generator_pool() when the attached
  * generator is a MultiStageScenarioGenerator. The outer dimension is
  * time_horizon; the middle dimension is the per-stage pool size (can
  * differ across stages); the inner dimension is the size of a single
  * sub-scenario for that stage, matching
  * scenario_set.get_sub_scenario_size(t).
  *
  * Mutually exclusive with #f_generator_pool_cache.
  *
  * The walking logic in prepare_multi_stage_generator_pool() assumes
  * stage independence; see the comment there. */

 std::vector< std::vector< std::vector< double > > > f_multi_stage_pool_cache;

 /// The start index of each admissible state
 /** For each t in {0, ..., TimeHorizon - 1}, admissible_state_begin[ t ] is
  * the index in vector #admissible_states at which the admissible state for
  * stage t begins. */

 std::vector< Index > admissible_state_begin;

 /// A vector containing the concatenation of states for each time instant
 std::vector< double > admissible_states;

 /// An initial state for the first stage problem
 std::vector< double > initial_state;

 /// Random cuts for each stage and each scenario
 /** This boost::multi_array stores the random cuts for all stages and all
  * scenarios. A random cut is a cut associated with a particular scenario. */
 boost::multi_array< PolyhedralFunction , 2 > random_cuts;

 /// It indicates whether the random cuts have been initialized
 bool f_random_cuts_initialized = false;

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
  * @param i The index of the sub-Block to be deserialized.
  *
  * @return A pointer to the Block that was deserialized.
  */
 Block * deserialize_sub_Block( const netCDF::NcGroup & group , Index i ) {
  std::string sub_group_name = "StochasticBlock_" + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );

  if( sub_group.isNull() ) {
   sub_group = group.getGroup( "StochasticBlock" );
   if( sub_group.isNull() )
    throw( std::logic_error( "SDDPBlock::deserialize: neither group '" +
                             sub_group_name + "' nor 'StochasticBlock' "
                             "was found." ) );
   sub_group_name = "StochasticBlock";
  }

  auto type = sub_group.getAtt( "type" );
  if( type.isNull() )
   throw( std::logic_error( "SDDPBlock::deserialize: attribute 'type' of '" +
                            sub_group_name + "' must be present." ) );

  std::string type_name;
  type.getValues( type_name );

  if( type_name != "StochasticBlock" )
   throw( std::logic_error( "SDDPBlock::deserialize: attribute 'type' of '" +
                            sub_group_name + "' must contain "
                            "'StochasticBlock'." ) );

  auto sub_Block = new_Block( sub_group , this );

  if( ! sub_Block )
   throw( std::logic_error( "SDDPBlock::deserialize: sub-group '" +
                            sub_group_name + "' is incomplete." ) );

  if( sub_group_name != "StochasticBlock" ) {

   // If StochasticBlock_i does not have sub-group "Block" then
   // "StochasticBlock" must have one.
   if( sub_group.getGroup( "Block" ).isNull() ) {

    auto StochasticBlock_group = group.getGroup( "StochasticBlock" );
    if( StochasticBlock_group.isNull() )
     throw( std::logic_error( "SDDPBlock::deserialize: sub-group 'Block' was "
                              "not provided neither in '" + sub_group_name +
                              "' nor in 'StochasticBlock'" ) );

    auto Block_group = StochasticBlock_group.getGroup( "Block" );
    if( Block_group.isNull() )
     throw( std::logic_error( "SDDPBlock::deserialize: sub-group 'Block' was "
                              "not provided neither in '" + sub_group_name +
                              "' nor in 'StochasticBlock'" ) );

    auto inner_block = new_Block( Block_group, this );
    if( ! inner_block )
     throw( std::logic_error( "SDDPBlock::deserialize: the 'Block' sub-group "
                              "of the 'StochasticBlock' group has an invalid "
                              "or incomplete description." ) );

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
         static_cast< StochasticBlock * >( sub_Block )->get_inner_block() );

      static_cast< StochasticBlock * >( sub_Block )->
       set_data_mappings( std::move( data_mappings ) );
     }
    }
   }
  }

  return( sub_Block );
 }

 /*--------------------------------------------------------------------------*/

 /// initializes the structure that stores the random cuts
 void initialize_random_cuts();

/*--------------------------------------------------------------------------*/

};   // end( class SDDPBlock )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SDDPBlockSolution --------------------------*/
/*--------------------------------------------------------------------------*/
/// a Solution of a SDDPBlock describing the trajectory of one scenario
/** The SDDPBlockSolution class, derived from Solution, represents a solution
 * of a SDDPBlock, i.e., the Solution of the inner Block of the
 * BendersBFunction at each stage (the first sub-Block of the stage if there
 * are more than one). Since a simulation (see SDDPGreedySolver) solves the
 * stage sub-problems in sequence for one given scenario, one
 * SDDPBlockSolution describes the complete trajectory of one scenario. */

class SDDPBlockSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend SDDPBlock;  ///< make SDDPBlock friend

/*------------- CONSTRUCTING AND DESTRUCTING SDDPBlockSolution -------------*/

 /// constructor

 explicit SDDPBlockSolution( void ) : Solution() ,
  f_inner_Config( nullptr ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a SDDPBlockSolution from a netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// destructor

 ~SDDPBlockSolution() override {
  for( auto si : v_stage_solutions )
   delete( si );
  delete( f_inner_Config );
  }

/*--------- METHODS DESCRIBING THE BEHAVIOR OF A SDDPBlockSolution ----------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a SDDPBlockSolution into a netCDF::NcGroup
 /** Serialize a SDDPBlockSolution into a netCDF::NcGroup. The format is the
  * following:
  *
  * - The dimension "TimeHorizon" containing the number of stages. The
  *   dimension is mandatory if "StageSolution_0" is there and optional
  *   otherwise.
  *
  * - If "TimeHorizon" is defined, the groups "StageSolution_T" for T = 0,
  *   ..., TimeHorizon - 1, each containing the Solution object of the inner
  *   Block of the BendersBFunction at stage T. The groups are optional but
  *   either they are all there or none is, hence one can just check the
  *   existence of StageSolution_0: if it exists then all other ones, and
  *   the dimension "TimeHorizon", must exist. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 SDDPBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 SDDPBlockSolution * clone( bool empty = false ) const override;

/*--------------------------------------------------------------------------*/
 /// set the inner Config
 /** Sets the f_inner_Config field, which is used in read() to Config-ure the
  * Solution of the inner Block of each stage; ownership of \p cfg is taken
  * by the SDDPBlockSolution. */

 void set_inner_Config( Configuration * cfg ) {
  delete( f_inner_Config );
  f_inner_Config = cfg;
  }

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "SDDPBlockSolution [" << this << "]: "
         << v_stage_solutions.size() << " stage solutions" << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 Configuration * f_inner_Config;
 ///< the Configuration for the Solution of the inner Block of each stage

 std::vector< Solution * > v_stage_solutions;
 ///< the Solution of the inner Block of each stage

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( SDDPBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* SDDPBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SDDPBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/

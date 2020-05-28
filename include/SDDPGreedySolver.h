/*--------------------------------------------------------------------------*/
/*---------------------- File SDDPGreedySolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SDDPGreedySolver class, implementing the Solver
 * interface, for multistage programming problems defined by the
 * SDDPBlock. The SDDPGreedySolver implements a greedy strategy to solve an
 * SDDPBlock for a fixed scenario as defined below.
 *
 * \version 0.1
 *
 * \date 28 - 05 - 2020
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

#ifndef __SDDPGreedySolver
#define __SDDPGreedySolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "Solver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SDDPGreedySolver_CLASSES Classes in SDDPGreedySolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS SDDPGreedySolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a greedy solver for multistage programming problems
/**
 * The SDDPGreedySolver class derives from Solver and implements a sequential,
 * greedy strategy to solve an SDDPBlock for a fixed scenario. Recall that an
 * SDDPBlock represents an optimization problem of the form
 *
 * \f[
 *   \min_{x_0 \in \mathcal{X}^{n_0}} f_0(x_0) +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_1 \in \mathcal{X}^{n_1}} f_1(x_1) +
 *   \mathbb{E} \left \lbrack \dots +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{X}^{n_{T-1}}} f_{T-1}(x_{T-1})
 *   \right\rbrack \right\rbrack\right\rbrack, \qquad (1)
 * \f]
 *
 * where T is the time horizon, \f$\mathcal{X}^{n_t} \equiv
 * \mathcal{X}^{n_t}(x_{t-1}, \xi_t) \subseteq \mathbb{R}^{n_t}\f$ for each
 * \f$t \in \{0, \dots, T-1\}\f$, and \f$ \xi = \{ \xi_t \}_{t \in \{1, \dots,
 * T-1\}} \f$ is a stochastic process. See SDDPBlock for details. The
 * SDDPGreedySolver considers the problem (1) for a single realization of the
 * stochastic process, i.e., a deterministic problem of the form
 *
 * \f[
 *   \min_{x_0 \in \mathcal{\tilde{X}}^{n_0}} f_0(x_0) +
 *   \left \lbrack
 *   \min_{x_1 \in \mathcal{\tilde{X}}^{n_1}} f_1(x_1) +
 *   \left \lbrack \dots +
 *   \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{\tilde{X}}^{n_{T-1}}} f_{T-1}(x_{T-1})
 *   \right\rbrack \right\rbrack\right\rbrack, \qquad (2)
 * \f]
 *
 * with \f$\mathcal{\tilde{X}}^{n_t} \equiv \mathcal{\tilde{X}}^{n_t}(x_{t-1},
 * \tilde{\xi}_t)\f$ where \f$ \tilde{\xi}_t = \{ \tilde{\xi}_t \}_{t \in \{1,
 * \dots, T-1\}} \f$ is a realization of the stochastic process \f$ \xi
 * \f$. The SDDPGreedySolver is a heuristic as it does not look for an optimal
 * solution to problem (2). The method it employs is very simple: it solves
 * the subproblem at each stage in sequence, from the first to the last one,
 * using the solution found for a stage to define the problem at the next
 * stage. First, for a given \f$ (x_{-1}, \tilde{\xi}_0)\f$, it solves
 * the problem
 *
 * @f{align}
 *   \min       & \ \ f_0(x_0) + \mathcal{P}_{1}(x_0) \qquad (3) \\
 *   {\rm s.t.} & \ \ x_0 \in \mathcal{\tilde{X}}^{n_0}(x_{-1},
 *                            \tilde{\xi}_0)
 * @f}
 *
 * where \f$ \mathcal{P}_{1} \f$ is a polyhedral function that approximates
 * the cost-to-go function. Let \f$ x_0^* \f$ be a solution obtained to
 * problem (3). Next, the following problem is solved
 *
 * @f{align}
 *   \min       & \ \ f_1(x_1) + \mathcal{P}_{2}(x_1)\\
 *   {\rm s.t.} & \ \ x_1 \in \mathcal{\tilde{X}}^{n_1}(x^*_{0},
 *                  \tilde{\xi}_1)
 * @f}
 *
 * and a solution \f$ x_1^* \f$ is obtained. This process continues until the
 * subproblem at stage \f$ T-1 \f$ is solved and a solution \f$ x_{T-1}^* \f$
 * is found for it. In general, for each \f$ t \in \{0, \dots, T-1\}\f$, the
 * subproblem solved at stage \f$ t \f$ is the following
 *
 * @f{align}
 *   \min       & \ \ f_t(x_t) + \mathcal{P}_{t+1}(x_t)\\
 *   {\rm s.t.} & \ \ x_t \in \mathcal{\tilde{X}}^{n_t}(x^*_{t-1},
 *                  \tilde{\xi}_t)
 * @f}
 *
 * where \f$ x^*_{t-1} \f$ is a solution to the subproblem at stage \f$ t-1
 * \f$ if \f$ t > 0 \f$, \f$ x^*_{-1} \equiv x_{-1} \f$, and \f$
 * \mathcal{P}_{t+1} \f$ denotes a polyhedral function that approximates the
 * cost-to-go function at stage \f$ t \f$.
 *
 *     Notice that the SDDPGreedySolver does not solve neither the problem
 *     enconded by SDDPBlock nor the deterministic (single-scenario)
 *     multistage problem defined in (2). It may not even find a feasible
 *     solution to problem (2) even if one exists.
 */

class SDDPGreedySolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 /// "import" Index from SDDPBlock
 using Index = SDDPBlock::Index;

 /// public enum for the possible return values of compute()
 /** Public enum "extending" Solver::sol_type with more detailed values
  * specific to SDDPGreedySolver. Among the values defined in
  * Solver::sol_type, only kOK is not considered in SDDPGreedySolver. A few
  * comments about the specific interpretation of these values is in order. As
  * explained above, the SDDPGreedySolver solves each subproblem
  * sequentially. The following values are returned by compute() depending on
  * the value returned by compute() when solving each subproblem.
  *
  * - kError is returned when an unrecoverable error occurs while solving a
  *   subproblem. The stage of the subproblem at which the error occurred can
  *   then be retrieved by the method get_fault_stage().
  *
  * - kUnbounded is returned when some subproblem is unbounded. The stage of
  *   the unbounded subproblem can then be retrieved by the method
  *   get_fault_stage(). It does not mean, however, that the subproblems at
  *   later stages are all feasible.
  *
  * - kInfeasible is returned when the subproblem at the first stage is
  *   infeasible. This is the only case in which the SDDPGreedySolver
  *   guarantees that the deterministic (single-scenario) multistage problem
  *   is infeasible. In this case, instead of returning kInfeasible we could
  *   return kSubproblemInfeasible (see below), but we decided to return
  *   kInfeasible as this value should be more broadly understood.
  *
  * - kStopTime is returned when all subproblems were "solved", but the
  *   solution of some subproblem terminated with a kStopTime status and all
  *   subproblems at previous stages terminated with either a kOK or a
  *   kLowPrecision status. The stage at which this event occurred can then be
  *   retrieved by the method get_fault_stage(). In other words,
  *   get_fault_stage() will return t such that the solution of subproblem at
  *   stage t terminated with a kStopTime status and the solution of each
  *   subproblem at stage in {0, ..., t-1} terminated with either a kOK or a
  *   kLowPrecision status.
  *
  * - kStopIter is returned when all subproblems were "solved", but the
  *   solution of some subproblem terminated with a kStopIter status and all
  *   subproblems at previous stages terminated with either a kOK or a
  *   kLowPrecision status. The stage at which this event occurred can then be
  *   retrieved by the method get_fault_stage(). In other words,
  *   get_fault_stage() will return t such that the solution of subproblem at
  *   stage t terminated with a kStopIter status and the solution of each
  *   subproblem at stage in {0, ..., t-1} terminated with either a kOK or a
  *   kLowPrecision status.
  *
  * - kLowPrecision is returned when every subproblem is "solved" and
  *   terminated with either a kOK or kLowPrecision status. In this case, the
  *   solution found is feasible for the deterministic (single-scenario)
  *   multistage problem but there is no guarantee that it is optimal.
  */

 enum sddp_greedy_sol_type {
 kSubproblemInfeasible ,  ///< some subproblem may be infeasible
 /**< Means that a subproblem at some stage, let say t, different than the
  * first one turned out to be infeasible. Since the feasible region of a
  * subproblem may depend on the solution of the subproblem at the previous
  * stage, it does not mean that the deterministic (single-scenario)
  * multistage problem is infeasible, as there could be another solution for
  * the subproblem at the previous stage that would make the subproblem at
  * stage t feasible. The stage t of the infeasible problem can be retrieved
  * by the method get_fault_stage().
  */

 };  // end( sddp_greedy_sol_type )

/**@} ----------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING SDDPGreedySolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SDDPGreedySolver
 *  @{ */

 /// constructor
 SDDPGreedySolver( void ) { }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~SDDPGreedySolver() { }

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the model encoded by the current Block
 *  @{ */

 /// (try to) solve the model encoded in the SDDPBlock for a single scenario
 /**
  */

 virtual int compute( bool changedvars = true ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 * @{ */

 virtual void get_var_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 /** If a call to compute() returns kError, kStopTime, or kStopIter, this
  * method returns the stage at which the associated event has
  * occurred. Otherwise, this method returns Inf<Index>().
  *
  * @return the stage at which a fault has occurred.
  */
 Index get_fault_stage() const {
  return fault_stage;
 }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE SDDPGreedySolver -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the state of the SDDPGreedySolver
 *  @{ */

 /// returns the time horizon of the problem associated with the SDDPBlock
 /** This function returns the time horizon of the problem associated with the
  * SDDPBlock with which this SDDPGreedySolver is attached.
  *
  * @return The time horizon of the problem associated with the SDDPBlock.
  */
 Index get_time_horizon( void ) const;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 /// The id of the scenario that should be considered
 Index scenario_id;

 /// The stage at which some special event has happened
 Index fault_stage;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int solve( Index stage , bool write_solution = false );

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

 /// returns a pointer to the BendersBFunction associated with the given stage
 /** This method returns a pointer to the BendersBFunction associated with the
  * given \p stage.
  *
  * @param[in] stage An Index in the interval [0, T-1], where T is the time
  * horizon.
  *
  * @return A pointer to the BendersBFunction associated with the given stage.
  */
 BendersBFunction * get_benders_function( Index stage ) const;

/*--------------------------------------------------------------------------*/

 void set_scenario( void );

/*--------------------------------------------------------------------------*/

 void set_initial_state( void );

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Indicates whether the scenario has already been set
 bool scenario_is_set = false;

 /// Indicates whether the initial state has already been set
 bool initial_state_is_set = false;

};   // end( class SDDPGreedySolver )

/** @} end( group( SDDPGreedySolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SDDPGreedySolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File SDDPGreedySolver.h ------------------------*/
/*--------------------------------------------------------------------------*/

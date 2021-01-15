/*--------------------------------------------------------------------------*/
/*---------------------- File SDDPGreedySolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file

 * Header file for the SDDPGreedySolver class, implementing the Solver
 * interface, for multistage programming problems defined by the
 * SDDPBlock. The SDDPGreedySolver implements a greedy strategy to try to
 * solve a deterministic (single-scenario) multistage problem encoded by an
 * SDDPBlock as defined below.
 *
 * \version 0.1
 *
 * \date 15 - 01 - 2021
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

#include "SDDPBlock.h"
#include "Solver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class BendersBFunction;      // forward declaration of BendersBFunction

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
 *     encoded by SDDPBlock nor the deterministic (single-scenario) multistage
 *     problem defined in (2). It may not even find a feasible solution to
 *     problem (2) even if one exists.
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
 kSubproblemInfeasible = kInfeasible + 1 ,
 ///< some subproblem may be infeasible
 /**< It means that a subproblem at some stage, let say t, different than the
  * first one turned out to be infeasible. Since the feasible region of a
  * subproblem may depend on the solution of the subproblem at the previous
  * stage, it does not mean that the deterministic (single-scenario)
  * multistage problem is infeasible, as there could be another solution for
  * the subproblem at the previous stage that would make the subproblem at
  * stage t feasible. The stage t of the infeasible problem can be retrieved
  * by the method get_fault_stage().
  */

 kSolutionNotFound ,
 ///< the solution to some subproblem has not been found
 /**< It means that a solution to a subproblem has not been found for whatever
  * reason. The stage at which the solution could not be found can be
  * retrieved by the method get_fault_stage().
  */

 };  // end( sddp_greedy_sol_type )

/*--------------------------------------------------------------------------*/

 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters of
  * "int" type that the SDDPGreedySolver has, besides those defined in
  * Solver. The value intLastAlgPar is provided so that the list can be easily
  * further extended by derived classes. */

 enum int_par_type_SDDP_Greedy_S {

  intScenarioId = int_par_type_S::intLastAlgPar ,
  ///< The id of the scenario that must be considered
  /**< This is the id of the scenario that must be considered when trying to
   * solve the deterministic (single-scenario) multistage problem. It must be
   * a valid id for a scenario handled by the SDDPBlock. */

  intLastAlgPar
  ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set of
   * int algorithmic parameters. */

 };  // end( int_par_type_SDDP_Greedy_S )

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

 /// set a given integer (int) numerical parameter
 /** Set a given integer (int) numerical parameter. Besides
  * considering the integer parameters defined in #int_par_type_S,
  * this function also accepts the following parameters:
  *
  * - #intScenarioId
  *
  * Please refer to the #int_par_type_SDDP_Greedy_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const int value ) override {
  switch( par ) {
   case( intScenarioId ): set_scenario_id( value ); return;
   case( intLogVerb ): log_verbosity = value; return;
  }
  Solver::set_par( par , value );
 }

/**@} ----------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the SDDPGreedySolver
 *  @{ */

 /// get the number of int parameters
 /** Get the number of int parameters.
  *
  * @return The number of int parameters.
  */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of an int parameter
 /** Get the default value of the int parameter with given index. Please see
  * the #int_par_type_SDDP_Greedy_S and #int_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 int get_dflt_int_par( const idx_type par ) const override {
  switch( par ) {
   case( intScenarioId ): return 0;
   case( intLogVerb ): return 0;
  }
  return Solver::get_dflt_int_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get a specific integer (int) numerical parameter
 /** Get a specific integer (int) numerical parameter. Please see the
  * #int_par_type_SDDP_Greedy_S and #int_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 int get_int_par( const idx_type par ) const override {
  switch( par ) {
   case( intScenarioId ): return scenario_id;
   case( intLogVerb ): return log_verbosity;
  }
  return( Solver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the index of the int parameter with given string \p name
 /** This method takes a string, which is assumed to be the name of an int
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it. The method is given a void
  * implementation (throwing exception), rather than being pure virtual, so
  * that derived classes not having any int parameter do not have to bother
  * with implementing it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type int_par_str2idx( const std::string & name ) const override {
  if( name == "intScenarioId" ) return intScenarioId;
  return Solver::int_par_str2idx( name );
 }

/*--------------------------------------------------------------------------*/
 /// returns the string name of the int parameter with given index
 /** This method takes an int parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & int_par_idx2str( const idx_type idx ) const override {

  static const std::vector<std::string> parameter_names = { "intScenarioId" };

  if( idx >= int_par_type_S::intLastAlgPar && idx < intLastAlgPar )
   return parameter_names[ idx - int_par_type_S::intLastAlgPar ];

  return Solver::int_par_idx2str( idx );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the model encoded by the current Block
 *  @{ */

 /// (try to) solve the model encoded in the SDDPBlock for a single scenario
 /** This method tries to solve the deterministic (single-scenario) multistage
  * problem defined in (2) above. The problem is determined by the scenario
  * whose id is returned by the get_scenario_id() method and whose initial
  * state is defined in the SDDPBlock. This method does not really try to
  * solve the problem (2), but employs a procedure that may find a feasible
  * solution for that problem and can be interpreted as a simulation.
  *
  * Beginning at the first stage, this method tries to solve the subproblem
  * associated with each stage, sequentially, until the last one. The
  * subproblem at stage t is encoded by the inner Block of the
  * BendersBFunction associated with the stage t. The subproblem at stage t >
  * 0 may depend on the variables of the subproblem at stage t-1. After
  * solving the subproblem at stage t-1, the solution found to this subproblem
  * is used to update the next subproblem according to this dependency (which
  * is characterized by the BendersBFunction).
  *
  * At any given stage, the subproblem may be successfully solved or not. If a
  * solution to a subproblem is not found (for instance, if the subproblem
  * turns out to be infeasible, unbounded, or an error occurs while solving
  * it), then this method stops with the corresponding status as described in
  * #sddp_greedy_sol_type. In this case, no solution for the deterministic
  * (single-scenario) multistage problem can be provided.
  *
  * Notice that a feasible solution may not be found for the deterministic
  * (single-scenario) multistage problem (2) even if one exists.
  *
  * @return Please refer to #sddp_greedy_sol_type for a description of each
  *         value that this method may return.
  */

 int compute( bool changedvars = true ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR CHANGING THE DATA OF THE SDDPGreedySolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the SDDPGreedySolver
 *  @{ */

 /// sets the scenario that should be considered
 /** This method defines which scenario should be considered when trying to
  * solve the deterministic single-scenario problem. The \p scenario_id
  * parameter must be the id of a scenario handled by the SDDPBlock.
  *
  * @param scenario_id The id of the scenario
  */
 void set_scenario_id( Index scenario_id ) {
  if( this->scenario_id == scenario_id )
   return;
  this->scenario_id = scenario_id;
  scenario_is_set = false;
  status_compute = Solver::kUnEval;
 }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 * @{ */

 bool has_var_solution( void ) override {
  return ( status_compute == Solver::kLowPrecision ) ||
   ( status_compute == Solver::kStopIter ) ||
   ( status_compute == Solver::kStopTime );
 }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override {
  if( ( get_objective_sense() == Objective::eMax ) && has_var_solution() )
   return solution_value;
  return( - std::numeric_limits<OFValue>::infinity() );
 }

/*--------------------------------------------------------------------------*/

 OFValue get_ub( void ) override {
  if( ( get_objective_sense() == Objective::eMin ) && has_var_solution() )
   return solution_value;
  return( std::numeric_limits<OFValue>::infinity() );
 }

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

 /// returns the time horizon of the SDDPBlock attached to this SDDPGreedySolver
 /** Returns the time horizon of the problem represented by the SDDPBlock
  * attached to this SDDPGreedySolver.
  *
  * @return The time horizon of the problem represented by the SDDPBlock
  *         attached to this SDDPGreedySolver.
  */
 Index get_time_horizon( void ) const;

/*--------------------------------------------------------------------------*/

 /// returns the id of the scenario being currently considered
 /** This method returns the id of the scenario being currently considered
  * when trying to solve the deterministic (single-scenario) multistage
  * problem.
  *
  * @return The id of the scenario to be considered
  */
 Index get_scenario_id( void ) const {
  return scenario_id;
 }

/*--------------------------------------------------------------------------*/

 /// returns the status of the most recent call to compute()
 /** Returns the status of the most recent call to compute().
  *
  * @return The status of the most recent call to compute().
  */
 Index get_status( void ) const {
  return status_compute;
 }

/*--------------------------------------------------------------------------*/

 /// sets the scenario to be considered
 /** This method updates the sub-Blocks of the SDDPBlock with the data
  * provided by the scenario whose id is given by the method
  * get_scenario_id().
  */
 void set_scenario( void );

/*--------------------------------------------------------------------------*/

 /// sets the initial state
 /** This function sets the state of the subproblem at the first stage
  * according to the initial state present in the SDDPBlock.
  */
 void set_initial_state( void );

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 /// The id of the scenario that should be considered
 Index scenario_id = 0;

 /// The stage at which some special event has happened
 Index fault_stage = Inf<Index>();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// solves the subproblem associated with the given stage
 /** This method solves the subproblem associated with the given \p stage. If
  * \p write_solution is true then the solution found for the subproblem (if
  * any) is written into its Block. To solve the subproblem, the method
  * compute() of the Solver attached to its Block is invoked and the status
  * returned by that method is returned here.
  *
  * @param stage The stage associated with the subproblem to be solved.
  *
  * @param write_solution Indicates whether the solution found for the
  *        subproblem (if any) must be written into its Block.
  *
  * @return The status returned by compute() when solving the subproblem.
  */
 int solve( Index stage , bool write_solution = false );

/*--------------------------------------------------------------------------*/

 /// returns the Solver attached to the subproblem at the given stage
 /** This method returns a pointer to the Solver attached to the subproblem
  * associated with the given \p stage.
  *
  * @param stage The stage associated with the subproblem whose Solver is
  *        desired.
  *
  * @return A pointer to the Solver attached to the subproblem associated with
  *         the given \p stage.
  */
 Solver * get_sub_solver( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// returns a pointer to the BendersBFunction associated with the given stage
 /** This method returns a pointer to the BendersBFunction associated with the
  * given \p stage.
  *
  * @param[in] stage An Index in the interval [0, T-1], where T is the time
  *            horizon.
  *
  * @return A pointer to the BendersBFunction associated with the given stage.
  */
 BendersBFunction * get_benders_function( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// returns the solution associated with the problem at the given stage
 /** This function returns the solution of the problem associated with the
  * given \p stage, which is part of the state variables of the next stage.
  *
  * @param stage The stage whose solution is required.
  *
  * @return The vector containing the solution of the problem at the given
  *         stage.
  */
 std::vector<double> get_solution( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// sets the state variables of the subproblem at the given stage
 /** This function sets the state variables associated with the subproblem at
  * the given \p stage.
  *
  * @param The vector containing the state of the problem at the given stage.
  *
  * @param stage The stage whose state must be set.
  */
 void set_state( const std::vector<double> & state , Index stage ) const;

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/

 double get_sub_solution_value( Index stage ) const {
  assert( stage < get_time_horizon() );
  auto sub_solver = get_sub_solver( stage );
  if( sub_solver->is_var_feasible() )
   return sub_solver->get_var_value();
  else if( get_objective_sense( stage ) == Objective::eMin )
   return sub_solver->get_ub();
  else
   return sub_solver->get_lb();
 }

/*--------------------------------------------------------------------------*/

 int get_objective_sense( Index stage = 0 ) const {
  assert( stage < get_time_horizon() );
  auto benders_function = get_benders_function( stage );
  auto inner_block = benders_function->get_inner_block();
  assert( inner_block );
  return inner_block->get_objective_sense();
 }

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// The status returned by compute()
 int status_compute = Solver::kUnEval;

 /// Indicates whether the scenario has already been set
 bool scenario_is_set = false;

 /// Indicates whether the initial state has already been set
 bool initial_state_is_set = false;

 /// The value of the solution (if any).
 double solution_value = 0.0;

 /// It indicates the level of verbosity of the log
 int log_verbosity = 0;

};   // end( class SDDPGreedySolver )

/** @} end( group( SDDPGreedySolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SDDPGreedySolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File SDDPGreedySolver.h ------------------------*/
/*--------------------------------------------------------------------------*/

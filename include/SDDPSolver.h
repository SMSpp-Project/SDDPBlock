/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SDDPSolver class, implementing the Solver interface,
 * for multistage linear stochastic programming problems defined by the
 * SDDPBlock. This is a wrapper for the SDDP method implemented by the
 * STochastic OPTimization library (StOpt):
 *
 * https://gitlab.com/stochastic-control/StOpt
 *
 * \version 0.1
 *
 * \date 16 - 12 - 2020
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

#ifndef __SDDPSolver
#define __SDDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <boost/bimap.hpp>
#include <Eigen/Dense>
#include "ScenarioSimulator.h"
#include "Solver.h"
#include "StOpt/sddp/OptimizerSDDPBase.h"
#include "StOpt/sddp/SDDPFinalCut.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

class SDDPBlock;
class StochasticBlock;

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SDDPSolver_CLASSES Classes in SDDPSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SDDPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// an SDDP solver for multistage linear stochastic programming problems
/**
 * The SDDPSolver class derives from Solver and implements the stochastic dual
 * dynamic programming (SDDP) method for multistage linear stochastic
 * problems. In fact, this works as a wrapper for the SDDP solver implemented
 * by the STochastic OPTimization (StOpt) library, which is publicly available
 * at https://gitlab.com/stochastic-control/StOpt. This SDDPSolver can be
 * attached to an SDDPBlock whose subproblems are linear, i.e., one that has
 * the following form:
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
 * We shall distinguish two types of random variables: the convex and
 * the non-convex random variables.
 *
 * - The *convex random variables* are those that can only appear in
 *   the right-hand side of the constraints, i.e., they can only be
 *   part of the vectors \f$ b_t \f$ for \f$ t \in \{1, \dots, T-1\}
 *   \f$.
 *
 * - The *non-convex random variables* are those that can only appear
 *   in the left-hand side of the constraints or in the objective
 *   function, i.e., they can only be part of the matrices \f$ A_t \f$
 *   or the vectors \f$ c_t \f$ for \f$ t \in \{1, \dots, T-1\} \f$.
 *
 * Among the convex random variables, we further consider two
 * types. The *time-related random variables* are those that depend on
 * some random data of the previous stage. The *time-independent
 * random variables* are those that do not depend on random variables
 * of previous stages. Notice that we do not allow non-convex random
 * variables to be time-related. We also denote the \f$t\f$-th random
 * variable of the stochastic process \f$ \xi \f$ as
 *
 * \f[  \xi_t = ( \omega_t , \nu_t ) \f]
 *
 * where \f$ \omega_t \f$ is a convex random vector and \f$ \nu_t \f$
 * is a non-convex random vector, for \f$ t \in \{1, \dots, T-1\}
 * \f$. The random vector \f$ \omega_t \f$, for \f$ t \in \{1, \dots,
 * T-1\} \f$, is divided into time-related and time-independent random
 * vectors as follows:
 *
 * \f[
 *    \omega_t = ( \omega_t^{\text{dep}} , \omega_t^{\text{ind}} ).
 * \f]
 *
 * This means that the first components of the random vector \f$
 * \omega_t \f$ are time-related random variables and the last
 * components are time-independent random variables. We also define
 * \f$ \omega_{-1}^{\text{dep}} \f$ in order for the first stage
 * problem to (artificially) present the same dependence as the other
 * subproblems do. For each \f$ t \in \{0, \dots, T-1\}\f$, we call
 *
 * \f[
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t}\\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{V}_{t+1}(x_t, \omega_t^{\text{dep}})
 * \f]
 *
 * the problem associated with stage \f$ t \f$, where
 *
 * \f[
 *    \mathcal{V}_{t+1}(x_t, \omega_t^{\text{dep}}) =
 *      \mathbb{E}
 *        \left\lbrack
 *          V_{t+1}(x_t, \xi_{t+1}) \mid \omega_t^{\text{dep}}
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
 *    \mathcal{V}_{t+1}(x_t, \omega_t^{\text{dep}})
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
 *    \mathcal{P}_{t+1}(x_t)
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
 */

class SDDPSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic
  * parameters of "int" type that the SDDP solver has, besides those
  * defined in Solver. The value intLastAlgPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum int_par_type_SDDP_S {

  intNStepConv = int_par_type_S::intLastAlgPar ,
  ///< Frequency in which the convergence is checked
  /**< The method stops when either the maximum number of iterations
   * is reached (an iteration performs one backward pass and one (or
   * two) forward pass(es); see parameter #intMaxIter for the maximum
   * number of iterations) or if the method converges according to
   * the criterion explained in the following. The convergence
   * criterion is checked every intNStepConv iterations. Thus, the
   * method stops if, at some iteration \f$ k \f$ that is a multiple
   * of intNStepConv, the following condition holds:
   *
   *   \f[
   *     \left | \frac{ \bar{z}_{k} - \underline{z}_{k} }
   *                  { \bar{z}_{k} } \right |         \le \epsilon.
   *   \f]
   *
   * where \f$ \epsilon \f$ is the value defined by the #dblAccuracy
   * parameter. \f$ \underline{z}_{k} \f$ is the lower bound obtained
   * at iteration \f$ k \f$, from the backward pass (such a lower
   * bound is computed at every iteration and requires no extra
   * effort). The upper bound \f$ \bar{z}_{k} \f$, on the other hand,
   * is only computed at iterations that are multiple of the value
   * given by the intNStepConv parameter. At such an iteration, an
   * extra forward pass is performed, in which the number of
   * simulations considered is given by the value of the parameter
   * #intNbSimulCheckForSimu. \f$ \bar{z}_{k} \f$ is the upper bound
   * computed by this forward pass. The default value for
   * intNStepConv is 1. */

  intPrintTime ,
  ///< Indicates whether computational time should be displayed
  /**< This parameter determines whether the computational time spent
   * at each step of the backward and forward passes should be
   * displayed. If the value of this parameter is zero, the
   * computation time is not displayed. Otherwise, the computational
   * time is displayed at every step. The default value for
   * intPrintTime is 1. */

  intNbSimulCheckForSimu ,
  ///< Number of simulations considered when checking convergence
  /**< This parameter determines the number of simulations that must
   * be considered during the forward pass that computes the upper
   * bound used for checking convergence. See the comments for the
   * parameter #intNStepConv for more details. The default value for
   * intNbSimulCheckForSimu is 1. */

  intLogVerbosity ,
  ///< It indicates the verbosity of the log
  /**< This parameter indicates the verbosity of the log. If it is less than
   * or equal to zero, no log is output. The higher this value, the more
   * detailed is the log. */

  intLastAlgPar
  ///< First allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes
   * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different types of algorithmic
  * parameters of "double" type that the SDDPSolver has, besides those
  * defined in Solver. The value dblLastAlgPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum dbl_par_type_SDDP_S {

  dblAccuracy = dbl_par_type_S::dblLastAlgPar ,
  ///< relative accuracy for declaring a solution optimal
  /**< The algorithmic parameter for setting the *relative* accuracy
   * required to the solution of the SDDPBlock. Please see the
   * comments of the #intNStepConv parameter for a detailed
   * explanation of its meaning. */

  dblLastAlgPar
  ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend
   * the set of double algorithmic parameters. */

 };  // end( dbl_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the string algorithmic parameters
 /** Public enum describing the different types of algorithmic
  * parameters of "string" type that the SDDPSolver has, besides those
  * defined in Solver. The value strLastAlgPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum str_par_type_SDDP_S {

  strRegressorsFilename = str_par_type_S::strLastAlgPar ,
  ///< name of the file that will store regressors
  /**< Name of the file in which the regressors will be stored. The
   * default value is "regressors.sddp". */

  strCutsFilename ,
  ///< name of the file that will store the cuts
  /**< Name of the file in which the cuts will be stored. The
   * default value is "cuts.sddp". */

  strVisitedStatesFilename ,
  ///< name of the file that will store the visited states
  /**< Name of the file in which the visited states will be stored.
   * The default value is "visited_states.sddp". */

  strLastAlgPar
  ///< first allowed new string parameter for derived classes
  /**< Convenience value for easily allow derived classes
   * to extend the set of string algorithmic parameters. */

 };  // end( str_par_type_SDDP_S )

/**@} ----------------------------------------------------------------------*/
/*---------------- CONSTRUCTING AND DESTRUCTING SDDPSolver -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SDDPSolver
 *  @{ */

 /// constructor
 SDDPSolver( void ) {

  // Default values of the parameters

  // int

  maximum_number_iterations = get_dflt_int_par( intMaxIter );
  convergence_frequency = get_dflt_int_par( intNStepConv );
  print_cpu_time = get_dflt_int_par( intPrintTime );
  number_simulations_for_convergence =
   get_dflt_int_par( intNbSimulCheckForSimu );

  // double

  accuracy = get_dflt_dbl_par( dblAccuracy );

  // string

  regressors_filename = get_dflt_str_par( strRegressorsFilename );
  cuts_filename = get_dflt_str_par( strCutsFilename );
  visited_states_filename = get_dflt_str_par( strVisitedStatesFilename );

  // vector

  // number_meshes = get_dflt_str_par( vecMeshForReg ); // TODO

  // SDDPOptimizer

  sddp_optimizer = std::make_shared<SDDPOptimizer>( this );
 }

/*--------------------------------------------------------------------------*/

 void set_Block( Block * block ) override {
  if( f_Block == block )  // registering to the same Block
   return;                // cowardly and silently return

  Solver::set_Block( block );

  if( ! block )
   return;

  if( auto sddp_block = dynamic_cast< SDDPBlock * >( block ) ) {
   auto scenario_set = sddp_block->get_scenario_set();
   std::static_pointer_cast< SDDPOptimizer >( sddp_optimizer )->
    set_scenarios( scenario_set );
  }
  else
   throw( std::invalid_argument( "SDDPSolver::set_Block: An SDDPSolver can "
                                 "only be attached to an SDDPBlock." ) );
 }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~SDDPSolver() { }

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
  * - #intMaxIter
  *
  * - #intNStepConv
  *
  * - #intPrintTime
  *
  * - #intNbSimulCheckForSimu
  *
  * - #intLogVerbosity
  *
  * Please refer to the #int_par_type_SDDP_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const int value ) override {
  switch( par ) {
  case( intMaxIter ): maximum_number_iterations = value; return;
  case( intNStepConv ): convergence_frequency = value; return;
  case( intPrintTime ): print_cpu_time = value; return;
  case( intNbSimulCheckForSimu ):
   number_simulations_for_convergence = value; return;
  case( intLogVerbosity ):
   log_verbosity = value; return;
  }
  Solver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set a given float (double) numerical parameter
 /** Set a given float (double) numerical parameter. Besides
  * considering the integer parameters defined in #dbl_par_type_S,
  * this function also accepts the following parameters:
  *
  * - #dblAccuracy
  *
  * - #dblLastAlgPar
  *
  * Please refer to the #dbl_par_type_SDDP_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const double value ) override {
  if( par == dblAccuracy ) {
   accuracy = value;
   return;
  }
  Solver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set a given string parameter
 /** Set a given string parameter. Besides considering the integer
  * parameters defined in #str_par_type_S, this function also accepts
  * the following parameters:
  *
  * - #strRegressorsFilename
  *
  * - #strCutsFilename
  *
  * - #strVisitedStatesFilename
  *
  * - #strLastAlgPar
  *
  * Please refer to the #str_par_type_SDDP_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const std::string & value ) override {
  switch( par ) {
  case( strRegressorsFilename ): regressors_filename = value; return;
  case( strCutsFilename ): cuts_filename = value; return;
  case( strVisitedStatesFilename ): visited_states_filename = value; return;
  }
  Solver::set_par( par , value );
 }

/**@} ----------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the SDDPSolver
 *  @{ */

 /// get the number of int parameters
 /** Get the number of int parameters.
  *
  * @return The number of int parameters.
  */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastAlgPar ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of double parameters
 /** Get the number of double parameters.
  *
  * @return The number of double parameters.
  */

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastAlgPar ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of string parameters
 /** Get the number of string parameters.
  *
  * @return The number of string parameters.
  */

 idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of an int parameter
 /** Get the default value of the int parameter with given index.
  * Please see the #int_par_type_SDDP_S and #int_par_type_S
  * enumerations for a detailed explanation of the possible
  * parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 int get_dflt_int_par( const idx_type par ) const override {
  switch( par ) {
  case( intNStepConv ): return 1;
  case( intPrintTime ): return 1;
  case( intNbSimulCheckForSimu ): return 1;
  case( intLogVerbosity ): return 0;
  }
  return Solver::get_dflt_int_par( par );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the default value of a double parameter
 /** Get the default value of the double parameter with given index.
  * Please see the #dbl_par_type_SDDP_S and #dbl_par_type_S
  * enumerations for a detailed explanation of the possible
  * parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 double get_dflt_dbl_par( const idx_type par ) const override {
  if( par == dblAccuracy ) return 1.0e-8;
  return Solver::get_dflt_dbl_par( par );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the default value of a string parameter
 /** Get the default value of the string parameter with given index.
  * Please see the #str_par_type_SDDP_S and #str_par_type_S
  * enumerations for a detailed explanation of the possible
  * parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::string & get_dflt_str_par( const idx_type par ) const override {

  static const std::vector<std::string> default_values =
   { "regressors.sddp" , "cuts.sddp" , "visited_states.sddp" };

  if( par >= str_par_type_S::strLastAlgPar && par < strLastAlgPar )
   return default_values[ par - str_par_type_S::strLastAlgPar ];

  return Solver::get_dflt_str_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get a specific integer (int) numerical parameter
 /** Get a specific integer (int) numerical parameter. Please see the
  * #int_par_type_SDDP_S and #int_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 int get_int_par( const idx_type par ) const override {
  switch( par ) {
  case( intMaxIter ): return maximum_number_iterations;
  case( intNStepConv ): return convergence_frequency;
  case( intPrintTime ): return print_cpu_time;
  case( intNbSimulCheckForSimu ): return number_simulations_for_convergence;
  case( intLogVerbosity ): return log_verbosity;
  }
  return( Solver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a specific float (double) numerical parameter
 /** Get a specific float (double) numerical parameter. Please see the
  * #dbl_par_type_SDDP_S and #dbl_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 double get_dbl_par( const idx_type par ) const override {
  if( par == dblAccuracy )
   return accuracy;
  return( get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a specific string numerical parameter
 /** Get a specific string numerical parameter. Please see the
  * #str_par_type_SDDP_S and #str_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::string & get_str_par( const idx_type par ) const override {
  switch( par ) {
  case( strRegressorsFilename ): return regressors_filename;
  case( strCutsFilename ): return cuts_filename;
  case( strVisitedStatesFilename ): return visited_states_filename;
  }
  return Solver::get_str_par( par );
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
  if( name == "intNStepConv" ) return intNStepConv;
  if( name == "intPrintTime" ) return intPrintTime;
  if( name == "intNbSimulCheckForSimu" ) return intNbSimulCheckForSimu;
  if( name == "intLogVerbosity" ) return intLogVerbosity;
  return Solver::int_par_str2idx( name );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the index of the double parameter with given string name
 /** This method takes a string, which is assumed to be the name of a double
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  if( name == "dblAccuracy" ) return dblAccuracy;
  return Solver::dbl_par_str2idx( name );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the index of the string parameter with given string name
 /** This method takes a string, which is assumed to be the name of a string
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type str_par_str2idx( const std::string & name ) const override {
  if( name == "strRegressorsFilename" ) return strRegressorsFilename;
  if( name == "strCutsFilename" ) return strCutsFilename;
  if( name == "strVisitedStatesFilename" ) return strVisitedStatesFilename;
  return Solver::str_par_str2idx( name );
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

  static const std::vector<std::string> parameter_names =
   { "intNStepConv", "intPrintTime", "intNbSimulCheckForSimu" ,
     "intLogVerbosity" };

  if( idx >= int_par_type_S::intLastAlgPar && idx < intLastAlgPar )
   return parameter_names[ idx - int_par_type_S::intLastAlgPar ];

  return Solver::int_par_idx2str( idx );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the string name of the double parameter with given index
 /** This method takes a double parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & dbl_par_idx2str( const idx_type idx ) const override {
  static const std::string dblAccuracy_name = "dblAccuracy";
  if( idx == dblAccuracy ) return dblAccuracy_name;
  return Solver::dbl_par_idx2str( idx );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the string name of the string parameter with given index
 /** This method takes a string parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & str_par_idx2str( const idx_type idx ) const override {

  static const std::vector<std::string> parameter_names =
   { "strRegressorsFilename", "strCutsFilename", "strVisitedStatesFilename" };

  if( idx >= str_par_type_S::strLastAlgPar && idx < strLastAlgPar )
   return parameter_names[ idx - str_par_type_S::strLastAlgPar ];

  return Solver::str_par_idx2str( idx );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the model encoded by the current Block
 *  @{ */

 /// (try to) solve the model encoded in the SDDPBlock
 /**
  */

 int compute( bool changedvars = true ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 * @{ */

 void get_var_solution( Configuration *solc = nullptr ) override {
  // TODO
 }

/*--------------------------------------------------------------------------*/

 double get_var_value( void ) override {
  return backward_value;
 }

/*--------------------------------------------------------------------------*/

 double get_lb( void ) override {
  return backward_value;
 }

/*--------------------------------------------------------------------------*/

 double get_ub( void ) override {
  return backward_value;
 }

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE SDDPSolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the state of the SDDPSolver
 *  @{ */

/*--------------------------------------------------------------------------*/

 /// returns the number of iterations performed in the last call to compute()
 /** This function returns the number of iterations performed by the solver
  * during the last call to compute().
  *
  * @return The number of iterations performed by the solver during the last
  *         call to compute().
  */
 int get_number_iterations_performed() const {
  return number_iterations_performed;
 }

/*--------------------------------------------------------------------------*/

 /// returns the time horizon of the problem associated with the SDDPBlock
 /** This function returns the time horizon of the problem associated with the
  * SDDPBlock with which this SDDPSolver is attached.
  *
  * @return The time horizon of the problem associated with the SDDPBlock.
  */
 SDDPBlock::Index get_time_horizon() const;

/*--------------------------------------------------------------------------*/

/// returns the solution associated with the problem at the given stage
/** This function returns the solution of the problem associated with the
 * given \p stage, which is part of the state variables of the next stage.
 *
 * @param stage The stage whose solution is required.
 *
 * @return The array containing the solution of the problem at the given
 *         stage.
 */

 template< class T = Eigen::ArrayXd >
 T get_solution( SDDPBlock::Index stage ) const;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/// returns a pointer to the BendersBFunction associated with the given \p stage
/** This function returns a pointer to the BendersBFunction associated with
 * the given \p stage, which must be an integer between 0 and
 * get_time_horizon() - 1.
 *
 * @param stage An index between 0 and get_time_horizon() - 1.
 */
 BendersBFunction * get_benders_function( SDDPBlock::Index stage ) const;

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 /// Initial state
 /** The initial state at the beginning of the simulation. */
 Eigen::ArrayXd initial_state;

 /// Number of meshes in each direction
 /** This array stores the number of meshes in each direction. The
  * i-th component of this array contains the number of meshes
  * (number of steps) at direction i. */
 Eigen::ArrayXi number_meshes;

 /// The cuts used at the last time step
 /** The cuts used at the last time step: when the final value
  * function is zero, the last cut is given by an all zero array of
  * size nbstate + 1. */
 StOpt::SDDPFinalCut final_cut;

 /// Number of iterations performed by the method
 /** Number of iterations performed by the method at the last call of
  * compute(). */
 int number_iterations_performed;

 /// Accuracy achieved by the method
 /** Accuracy achieved by the method at the last call to compute(),
  * which is given by
  *
  * | backwardValue - forwardValue | / forwardValue
  *
  * where backwardValue is the value of the last backward pass and
  * forwardValueForConv is the value obtained in the forward pass
  * when checking for convergence.
  */
 double accuracy_achieved;

 /// It indicates the level of verbosity of the log
 int log_verbosity = 0;

 // PARAMETERS

 /// Name of the file in which regressors will be stored
 std::string regressors_filename;

 /// Name of the file in which the cuts will be stored
 std::string cuts_filename;

 /// Name of the file in which the visited states will be stored
 std::string visited_states_filename;

 /// Maximum number of iterations that the method should perform
 int maximum_number_iterations;

 /** Frequency in which the convergence check is performed. The
  * convergence is checked every "convergence_frequency" steps of the
  * method. See the comments about the intNStepConv parameter for
  * more details. */
 int convergence_frequency;

 /** Indicates whether the CPU time spent at each backward and
  * forward steps should be printed. */
 bool print_cpu_time;

 /** The number of simulations (number of forward passes called) when
  * we have to check the convergence by comparing the outcome given
  * by the forward pass and the one given by the backward pass. */
 int number_simulations_for_convergence;

 /** Relative accuracy for declaring a solution optimal. See the
  * comments about the dblAccuracy parameter for more details. */
 double accuracy;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/// add cuts to the sub-problem at the given stage
/** This function adds cuts to the sub-problem at the given \p stage.
 *
 * @param cuts An Eigen::ArrayXXd containing the cuts to be added. It must be
 *        a matrix with as many columns as there are cuts to be added and the
 *        number of rows must be equal to one plus the number of Variable
 *        defined in the BendersBlock associated with stage \p stage.
 *
 * @param stage The stage at which cuts should be updated, which must be an
 *        integer between 0 and get_time_horizon() - 1.
 *
 * @param range The indices of the cuts in \p cuts that should be added.
 */

 void add_cuts( const Eigen::ArrayXXd & cuts , SDDPBlock::Index stage ,
                Block::Range range =
                std::make_pair( 0 , Inf<SDDPBlock::Index>() )  ) const;

/*--------------------------------------------------------------------------*/

 void set_state( const Eigen::ArrayXd & state ,
                 SDDPBlock::Index stage ) const;

/*--------------------------------------------------------------------------*/

 void set_scenario( SDDPBlock::Index scenario_id ,
                    SDDPBlock::Index stage ) const;

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification();

/*--------------------------------------------------------------------------*/

/// solves the subproblem associated with the given stage
/** This function solves the subproblem associated with the given \p stage,
 * which must be an integer between 0 and get_time_horizon() - 1.
 *
 * @param stage The stage whose associated subproblem must be solved.
 */
 double solve( SDDPBlock::Index stage );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the value of the last backward pass
 double backward_value;

 /// the value obtained during the forward pass when checking for convergence
 double forward_value;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE CLASSES -------------------------------*/
/*--------------------------------------------------------------------------*/

 class SDDPOptimizer : public StOpt::OptimizerSDDPBase {

 public:

  /// constructor taking a pointer to an SDDPSolver
  /** Constructs an SDDPOptimizer associated with the given SDDPSolver.
   *
   * @param solver A pointer to an SDDPSolver. This parameter is
   * optional and its default value is nullptr.
   */
  SDDPOptimizer( SDDPSolver * solver = nullptr ) {
   sddp_solver = solver;
  }

/*--------------------------------------------------------------------------*/

  Eigen::ArrayXd oneStepBackward
  ( const StOpt::SDDPCutOptBase & p_linCut,
    const std::tuple< std::shared_ptr<Eigen::ArrayXd>, int, int > & p_aState,
    const Eigen::ArrayXd & p_particle, const int & p_isample) const override;

/*--------------------------------------------------------------------------*/

  double oneStepForward
  ( const Eigen::ArrayXd &p_aParticle, Eigen::ArrayXd &p_state,
    Eigen::ArrayXd &p_stateToStore,
    const StOpt::SDDPCutOptBase &p_linCut,
    const int &p_isimu ) const override;

/*--------------------------------------------------------------------------*/

  /// updates this SDDPOptimizer for a new stage
  /** This function updates this SDDPOptimizer for a new stage. The
   * \p date and \p date_next parameters have different meanings in
   * the backward and forward steps:
   *
   * - In a backward step, the optimization problem that must be
   *   solved is that associated with the stage \p date_next. The
   *   parameter \p date indicates, therefore, the previous stage.
   *
   * - In a forward step, the optimization problem that must be
   *   solved is that associated with the stage \p date. The
   *   parameter \p date_next indicates, therefore, the next stage.
   *
   *   We assume that the given arguments correspond to consecutive
   *   stages, i.e., \p date_next == \p date + 1. Moreover, we
   *   assume that -1 <= \p date < T, where T is the time
   *   horizon. Also, the following conditions should hold:
   *
   * -# -1 <= \p date <= T - 2 for any backward step;
   *
   * -#  0 <= \p date <= T - 1 for any forward step.
   *
   * @param date A stage.
   *
   * @param date_next Another stage.
   */
  void updateDates( const double & date, const double & date_next ) override;

/*--------------------------------------------------------------------------*/

  /// returns an initial state for the problem at the given stage
  /** This function must return an initial state for the
   * optimization problem associated with the given date. If the
   * given date is t, then the optimization problem associated with
   * time t has variables x_t and depends on the state (x_{t-1},
   * w_{t-1}^{dep}). Thus, this function must return an array
   * containing values for x_{t-1} and w_{t-1}^{dep}. Notice that
   * the subvector w_{t-1}^{dep} is only present in the state if
   * some random variables of the problem associated with time t
   * depend on the random variables w_{t-1}^{dep} of the problem at
   * stage t-1.
   *
   * @param stage The stage for which an initial state must be
   * provided.
   *
   * @return An initial state for the optimization problem
   * associated with the given stage.
   */
  Eigen::ArrayXd oneAdmissibleState( const double & stage ) override;

/*--------------------------------------------------------------------------*/

  /// return the size of the state vector
  /** This function returns the size of the state vector. It assumes
   * that the states of all stages have the same size. If the state
   * at a time t is given by (x_t, w_t^{dep}), then this function
   * should return the size of x_t plus the size of w_t^{dep}.
   *
   * @return The size of the state vector.
   */
  int getStateSize() const override {
   return sddp_solver->initial_state.size();
  }

/*--------------------------------------------------------------------------*/

  /// returns the simulator for the backward pass
  /** This function returns the simulator that is used during the
   * backward pass.
   *
   * @return The simulator associated with the backward pass.
   */
  std::shared_ptr< StOpt::SimulatorSDDPBase >
  getSimulatorBackward() const override {
   return simulator_backward;
  }

/*--------------------------------------------------------------------------*/

  /// returns the simulator for the forward pass
  /** This function returns the simulator that is used during the
   * forward pass.
   *
   * @return The simulator associated with the forward pass.
   */
  std::shared_ptr< StOpt::SimulatorSDDPBase >
  getSimulatorForward() const override {
   return simulator_forward;
  }

/*--------------------------------------------------------------------------*/

  /// set the SDDPSolver with which this SDDPOptimizer will be associated
  /** This method is used to set the (pointer to the) SDDPSolver
   * with which this SDDPOptimizer will be associated.
   *
   * @param solver A pointer to an SDDPSolver.
   */
  void set_solver( SDDPSolver * solver ) {
   sddp_solver = solver;
  }

/*--------------------------------------------------------------------------*/

  void set_scenarios( const ScenarioSet & scenario_set ) {
   if( simulator_forward )
    simulator_forward->set_scenarios( scenario_set );
   else {
    simulator_forward = std::make_shared< ScenarioSimulator >( scenario_set ,
                                                               false );
    simulator_forward->set_number_simulations
     ( std::min( 3u , scenario_set.size() ) );
   }
   if( simulator_backward )
    simulator_backward->set_scenarios( scenario_set );
   else {
    simulator_backward = std::make_shared< ScenarioSimulator >( scenario_set ,
                                                                true );
    simulator_backward->set_number_simulations( scenario_set.size() );
   }
  }

/*--------------------------------------------------------------------------*/

  /// returns the (pointer to the) Block associated with the given stage
  /** This method takes a double as parameter, representing a stage, and
   * returns a pointer to the Block associated with this stage.
   *
   * @param[in] stage A number in the interval [0, T-1], where T is the time
   *            horizon. Although the parameter is of type double, its value
   *            must be actually an integer.
   *
   * @return A pointer to the Block associated with the given stage.
   */

  StochasticBlock * get_block( const double & stage ) const;

/*--------------------------------------------------------------------------*/

 private:

  double date;
  double date_next;
  std::shared_ptr< ScenarioSimulator > simulator_forward;
  std::shared_ptr< ScenarioSimulator > simulator_backward;

  SDDPSolver * sddp_solver;
 };   // end( class SDDPOptimizer )

/*--------------------------------------------------------------------------*/

 friend SDDPOptimizer;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 std::shared_ptr< StOpt::OptimizerSDDPBase > sddp_optimizer;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

};   // end( class SDDPSolver )

/** @} end( group( SDDPSolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SDDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SDDPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/

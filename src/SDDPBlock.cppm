/*--------------------------------------------------------------------------*/
/*-------------------------- File SDDPBlock.cppm ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.sddpblock, which exports every name
 * that the headers of SDDPBlock declare at namespace scope, together with the
 * modules it depends on; a translation unit can then write
 *
 *     import smspp.sddpblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
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
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "ParallelSDDPSolver.h"
#include "SDDPGreedySolver.h"
#include "SDDPSimulationSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.sddpblock;

export import smspp.stochasticblock;

export namespace SMSpp_di_unipi_it {

 // ParallelSDDPSolver.h
 using SMSpp_di_unipi_it::ParallelSDDPSolver;

 // SDDPBlock.h
 using SMSpp_di_unipi_it::SDDPBlock;
 using SMSpp_di_unipi_it::SDDPBlockSolution;

 // SDDPGreedySolver.h
 using SMSpp_di_unipi_it::SDDPGreedySolver;

 // SDDPSimulationSolver.h
 using SMSpp_di_unipi_it::SDDPSimulationSolver;

 // SDDPSolver.h
 using SMSpp_di_unipi_it::SDDPSolver;
 using SMSpp_di_unipi_it::SDDPSolverState;

 // ScenarioSet.h
 using SMSpp_di_unipi_it::ScenarioSet;

 // ScenarioSimulator.h
 using SMSpp_di_unipi_it::ScenarioSimulator;
}

/*--------------------------------------------------------------------------*/
/*------------------------ End File SDDPBlock.cppm -------------------------*/
/*--------------------------------------------------------------------------*/

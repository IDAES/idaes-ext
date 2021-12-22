/*###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
###############################################################################*/

// AMPL solver interface for PETSc
// Author: John Eslick

#ifndef PETSC_H
#define PETSC_H

#include <sys/types.h> // fix an issue with ssize_t on CentOS 7 build
#include<asl.h>
#undef filename
#include<petscsnes.h>
#include<petscts.h>

#define MSG_BUF_SIZE 2000

/* Define some default settings, changable through commandline options */
#define DEFAULT_LINEAR_PACK MATSOLVERMUMPS //be sure to build petsc with mumps
#define DEFAULT_SNES SNESNEWTONLS // newton line search
#define DEFAULT_PC PCLU //LU decomposition direct solve
#define DEFAULT_KSP KSPPREONLY //default preconditioner is direct solve
#define DEFAULT_SNES_ATOL 1e-9
#define DEFAULT_SNES_RTOL 1e-9
#define DEFAULT_KSP_ATOL 1e-10
#define DEFAULT_KSP_RTOL 1e-15
#define DEFAULT_SNES_MAX_IT 2000
#define DEFAULT_SNES_MAX_FUNC 50000
#define DEFAULT_TS_SNES_MAX_FAIL 500

typedef enum{  //keep these under 50 and shouldn't confilict with PETSc codes
   P_EXIT_NORMAL = 0, //Finished okay (solved is another matter)
   P_EXIT_INTEGER = 1, //Exited due to integer variables
   P_EXIT_DOF = 2, //Exited on DOF != 0
   P_EXIT_INEQ = 3, //Exited on inequalities
   P_EXIT_NO_NL_FILE_ERROR = 5, //Exited due to file not specified
   P_EXIT_NL_FILE_ERROR = 6, //Exited due to nl file read error
   P_EXIT_DOF_DAE = 7, //DOF wrong for DAE problem
   P_EXIT_VAR_DAE_MIS = 8, //number of derivatives mismatch
   P_EXIT_MULTIPLE_TIME = 9 //more than one time variable
}P_EXIT_CODES;

typedef struct{ // Sturcture for solver options
  PetscMPIInt    mpi_size; //Number of processors (should be 1 for now)
  PetscBool      show_cl; //show the command line, and transformed CL
  PetscBool      ignore_scaling;
  char           stub[PETSC_MAX_PATH_LEN]; //File name (with or without ext)
  char           typ_file[PETSC_MAX_PATH_LEN]; //output file with DAE types
  fint           stublen; // Stub string length
  PetscBool      got_stub;  // file stub was specified with -s
  PetscBool      show_con;  // Option to show initial constraint values
  PetscBool      show_init; // show initial values for x vec
  PetscBool      show_jac;  // show jacobian at intial value
  PetscBool      show_scale_factors;  // show jacobian at intial value
  PetscBool      ampl_opt;  // -AMPL specified I catch it but ignore
  PetscBool      use_bounds; // give solver variable bounds
  PetscBool      scale_var; // scale the variables based on jacobian at init
  PetscBool      scale_eq;  // scale the equations based on jacobian at init
  PetscBool      dae_solve; //use dae solver (requires suffix information)
}Solver_options;

typedef struct{ // solver context
  ASL *asl; // ASL context
  Solver_options opt; // command-line options
  SufDesc *dae_suffix_var; // DAE suffixes on variables
  SufDesc *dae_link_var; // DAE link derivatives to vars
  SufDesc *scaling_factor_var; //user scale factors for vars
  SufDesc *scaling_factor_con; //user scale factors for constraints
  int dae_map_t; // ASL index of time variable (-1 for none)
  int *dae_map_x; // PETSc index in x vec -> ASL index
  int *dae_map_xdot; // PETSc index in xdot vec -> ASL index
  int *dae_map_back; // ASL var index -> PETSc index (in x or xdot)
  int *dae_link; // ASL index -> ASL index of linked derivative or differential var
  int n_ineq; // Number of inequality constraints
  int n_var_diff; //DAE number of differential vars
  int n_var_deriv; //DAE number of derivative vars
  int n_var_state; //DAE number of state vars
  int n_var_alg;  //DAE number of algebraic variables
  int explicit_time; //DAE includes time variable? 1=yes 0=no
  int dof; // Degrees of freedom
  FILE *nl; // nl-file pointer
}Solver_ctx;

// Initialize a solver context
void sol_ctx_init(Solver_ctx *ctx);

// Transform arguments from AMPL-style and double dash to form expected by PETSc
char **transform_args(int argc, char** argv, int *size);

// Read scaling suffixes and set variable and constraint scaling as appropriate
int ScaleVarsUser(Solver_ctx *sol_ctx);
int ScaleEqsUser(Solver_ctx *sol_ctx);

// SNES Function and Jacobian callbacks
PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
PetscErrorCode FormFunction(SNES,Vec,Vec,void*);

// TS Function and Jacobian callbacks
PetscErrorCode FormDAEFunction(TS, PetscReal, Vec, Vec, Vec, void*);
PetscErrorCode FormDAEJacobian(TS, PetscReal, Vec, Vec, PetscReal, Mat, Mat, void*);

// Get DAE suffixes and map variables relationships
void get_dae_info(Solver_ctx *sol_ctx);
void dae_var_map(Solver_ctx *sol_ctx);

// get solver status for sol file
int get_snes_sol_message(char *msg, SNESConvergedReason term_reason, ASL *asl);

// Extra diagnostic printing functions.
void print_commandline(const char* msg, int argc, char **argv);
void print_x_asl(ASL *asl);
void print_jac_asl(ASL *asl);
void print_var_scale_factors_asl(ASL *asl);
void print_con_scale_factors_asl(ASL *asl);
void print_init_diagnostic(Solver_ctx *sol_ctx);

#endif

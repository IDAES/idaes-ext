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

// TODO<jce> TS Solver return code
// TODO<jce> Add TAO optimization solvers

#include"help.h"
#include"petsc.h"
#include<stdio.h>
#include<math.h>

int main(int argc, char **argv){
  PetscErrorCode ierr;        // Error code from PETSc functions
  PetscInt       its;         // Number of solver iterations
  ASL            *asl;        // ASL context
  Solver_ctx     sol_ctx;     // solver context
  int            err;         // Error code  from ASL fulctions
  int            i=0;         // Loop counters
  int            argc_new;    // new number of arguments reformated for PETSc
  char           **argv_new;  // argv transformed to PETSc's format
  static SufDecl suftab[] = { // suffixes to read in
    //doc for this at https://ampl.com/netlib/ampl/solvers/README.suf
    {"dae_suffix", NULL, ASL_Sufkind_var}, //var kinds for DAE solver
    {"dae_link", NULL, ASL_Sufkind_var}, //link derivatives to vars
    {"scaling_factor", NULL, ASL_Sufkind_var|ASL_Sufkind_real}, //var scale factors
    {"scaling_factor", NULL, ASL_Sufkind_con|ASL_Sufkind_real} //constraint scale factors
  };
  Vec            x,r,xl,xu; // solution, residual, bound vectors
  Mat            J;            // Jacobian matrix
  TS             ts;           // DAE solver
  SNES           snes;         // nonlinear solver context
  KSP            ksp;          // linear solver context
  PC             pc;           // linear preconditioner context
  SNESLineSearch linesearch;   // line search context
  SNESConvergedReason cr; // reason for convergence (or lack of convergence)
  char           msg[MSG_BUF_SIZE]; // just a string buffer
  PetscScalar    *xx, *xxl, *xxu; // for accessing x, xlb, and xub vectors
  real t; //time for DAE solution
  TSConvergedReason tscr;
  real *x_asl;
  PetscViewer pv;

  // Set some initial values in sol_ctx
  sol_ctx_init(&sol_ctx);
  sol_ctx.opt.stub[0] = '\0';
  // Change the AMPL style args into what PETSc would expect
  argv_new = transform_args(argc, argv, &argc_new);
  // Initialize up PETSc stuff
  PetscInitialize(&argc_new, &argv_new, (char*)0, help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &sol_ctx.opt.mpi_size);CHKERRQ(ierr);
  // Get added options
  PetscOptionsGetString(NULL, NULL, "-s", sol_ctx.opt.stub, PETSC_MAX_PATH_LEN-1, &sol_ctx.opt.got_stub);
  PetscOptionsHasName(NULL, NULL, "-show_jac", &sol_ctx.opt.show_jac);
  PetscOptionsHasName(NULL, NULL, "-show_initial", &sol_ctx.opt.show_init);
  PetscOptionsHasName(NULL, NULL, "-show_scale_factors", &sol_ctx.opt.show_scale_factors);
  PetscOptionsHasName(NULL, NULL, "-use_bounds", &sol_ctx.opt.use_bounds);
  PetscOptionsHasName(NULL, NULL, "-AMPL", &sol_ctx.opt.ampl_opt); // I don't use this
  PetscOptionsHasName(NULL, NULL, "-dae_solve", &sol_ctx.opt.dae_solve);
  PetscOptionsHasName(NULL, NULL, "-show_cl", &sol_ctx.opt.show_cl);
  PetscOptionsHasName(NULL, NULL, "-ignore_scaling", &sol_ctx.opt.ignore_scaling);

  // If show_cl otion, show the original and transformed command line
  if(sol_ctx.opt.show_cl){
    PetscPrintf(PETSC_COMM_SELF, "-----------------------------------------------------------------\n");
    print_commandline("Original Exec:\n  ", argc, argv);
    print_commandline("Transformed Exec:\n  ", argc_new, argv_new);
    PetscPrintf(PETSC_COMM_SELF, "-----------------------------------------------------------------\n");
  }
  // Make sure that a file was specified, and get the string length
  if(argc < 2) exit(P_EXIT_NL_FILE_ERROR);
  if(!sol_ctx.opt.got_stub) strcpy(sol_ctx.opt.stub, argv[1]); //assume first arg is file if no -s
  if(sol_ctx.opt.stub[0]=='-') exit(P_EXIT_NL_FILE_ERROR); // is an option name not filename
  if(strlen(sol_ctx.opt.stub)==0) exit(P_EXIT_NL_FILE_ERROR);
  sol_ctx.opt.stublen = strlen(sol_ctx.opt.stub);
  // Create ASL context and read nl file
  sol_ctx.asl = ASL_alloc(ASL_read_fg); // asl context
  asl = sol_ctx.asl;
  // set the suffix data structure
  suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
  // get file pointer and basic problem info from nl file
  sol_ctx.nl = jac0dim(sol_ctx.opt.stub, sol_ctx.opt.stublen);
  if(sol_ctx.nl==NULL){
      PetscPrintf(PETSC_COMM_SELF, "Could not read nl file %s\n", sol_ctx.opt.stub);
      exit(P_EXIT_NL_FILE_ERROR);
  }
  // Allocated space for some ASL items and test calculations
  X0 = (real*)Malloc(n_var*sizeof(real));  /* Initial X values */
  LUv = (real*)Malloc(n_var*sizeof(real)); /* Variable lower bounds */
  Uvx = (real*)Malloc(n_var*sizeof(real)); /* Variable upper bounds */
  LUrhs = (real*)Malloc(n_con*sizeof(real)); /* Lower constraint right side */
  Urhsx = (real*)Malloc(n_con*sizeof(real)); /* Upper constraint right side */

  // Read nl file and make function for jacobian
  err = fg_read(sol_ctx.nl, 0);
  PetscPrintf(PETSC_COMM_SELF, "Called fg_read, err: %d (0 is good)\n", err);

  // count inequalities
  for(i=0; i<n_con; ++i) if(LUrhs[i] - Urhsx[i] > 1e-9) {
    ++sol_ctx.n_ineq;
    PetscPrintf(PETSC_COMM_SELF, "%d, ineq (%f < body < %f)", i, LUrhs[i], Urhsx[i]);
  }
  // count degrees of freedom (n_var and n_con are macros from asl.h)
  if(sol_ctx.opt.dae_solve){
    get_dae_info(&sol_ctx);
    dae_var_map(&sol_ctx);
    sol_ctx.dof = n_var - n_con + sol_ctx.n_ineq -
      sol_ctx.n_var_deriv - sol_ctx.explicit_time;
  }
  else{
    sol_ctx.dof = n_var - n_con + sol_ctx.n_ineq;
  }
  // Print basic problem information
  PetscPrintf(PETSC_COMM_SELF, "---------------------------------------------------\n");
  PetscPrintf(PETSC_COMM_SELF, "DAE: %d\n", sol_ctx.opt.dae_solve);
  PetscPrintf(PETSC_COMM_SELF, "Reading nl file: %s\n", sol_ctx.opt.stub);
  PetscPrintf(PETSC_COMM_SELF, "Number of constraints: %d\n", n_con);
  PetscPrintf(PETSC_COMM_SELF, "Number of nonlinear constraints: %d\n", nlc);
  PetscPrintf(PETSC_COMM_SELF, "Number of linear constraints: %d\n", n_con-nlc);
  PetscPrintf(PETSC_COMM_SELF, "Number of inequalities: %d\n", sol_ctx.n_ineq);
  PetscPrintf(PETSC_COMM_SELF, "Number of variables: %d\n", n_var);
  PetscPrintf(PETSC_COMM_SELF, "Number of integers: %d\n", niv);
  PetscPrintf(PETSC_COMM_SELF, "Number of binary: %d\n", nbv);
  PetscPrintf(PETSC_COMM_SELF, "Number of objectives: %d (Ignoring)\n", n_obj);
  PetscPrintf(PETSC_COMM_SELF, "Number of non-zeros in Jacobian: %d \n", nzc);
  // If DAES, get DAE var types and map vars between ASL and PETSc
  if(sol_ctx.opt.dae_solve){
    PetscPrintf(PETSC_COMM_SELF, "Explicit time variable: %d\n", sol_ctx.explicit_time);
    PetscPrintf(PETSC_COMM_SELF, "Number of derivatives: %d\n", sol_ctx.n_var_deriv);
    PetscPrintf(PETSC_COMM_SELF, "Number of differential vars: %d\n", sol_ctx.n_var_diff);
    PetscPrintf(PETSC_COMM_SELF, "Number of algebraic vars: %d\n", sol_ctx.n_var_alg);
    PetscPrintf(PETSC_COMM_SELF, "Number of state vars: %d\n", sol_ctx.n_var_state);
  }
  PetscPrintf(PETSC_COMM_SELF, "Number of degrees of freedom: %d\n", sol_ctx.dof);
  PetscPrintf(PETSC_COMM_SELF, "---------------------------------------------------\n");

  // There are some restrictions (at least for now) to check
  if(nbv + niv > 0){ // no integer vars (nbv and niv are ASL macros)
    PetscPrintf(PETSC_COMM_SELF, "ERROR: Contains integer or binary variables.");
    ASL_free(&(sol_ctx.asl));
    exit(P_EXIT_INTEGER);
  }
  if(sol_ctx.dof != 0){ //dof must == 0 for nonlinear solve
    PetscPrintf(PETSC_COMM_SELF, "ERROR: Degrees of freedom not equal to 0\n");
    ASL_free(&(sol_ctx.asl));
    exit(P_EXIT_DOF);
  }
  if(sol_ctx.n_ineq > 0){ // no inequalities for nonlinear sys or DAE
    PetscPrintf(PETSC_COMM_SELF, "ERROR: contains inequalities");
    ASL_free(&(sol_ctx.asl));
    exit(P_EXIT_INEQ);
  }
  if(sol_ctx.explicit_time > 1){
    PetscPrintf(PETSC_COMM_SELF, "ERROR: DAE: Multiple time variables (allowed 1 at most)");
    ASL_free(&(sol_ctx.asl));
    exit(P_EXIT_MULTIPLE_TIME);
  }

  // Equation/variable scaling
  sol_ctx.dae_suffix_var = suf_get("dae_suffix", ASL_Sufkind_var);
  sol_ctx.dae_link_var = suf_get("dae_link", ASL_Sufkind_var);
  sol_ctx.scaling_factor_var = suf_get("scaling_factor", ASL_Sufkind_var);
  sol_ctx.scaling_factor_con = suf_get("scaling_factor", ASL_Sufkind_con);

  // Equation/variable scaling
  if(!sol_ctx.opt.ignore_scaling){
    ScaleVarsUser(&sol_ctx);
    ScaleEqsUser(&sol_ctx);
  }

  if(sol_ctx.opt.dae_solve){  //This block sets up DAE solve and solves
    ierr = TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr); //create an x vector
    ierr = VecSetSizes(x, PETSC_DECIDE, sol_ctx.n_var_state); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr); //command line options for vec
    ierr = VecDuplicate(x, &r);CHKERRQ(ierr);  // duplicate x for resuiduals
    ierr = VecDuplicate(x,&xl);CHKERRQ(ierr); // duplicate x for lower bounds
    ierr = VecDuplicate(x,&xu);CHKERRQ(ierr); // duplicate x for upper bounds
    /* Make x vec set initial guess from nl file, also get lb and ub */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(xl,&xxl);CHKERRQ(ierr);
    ierr = VecGetArray(xu,&xxu);CHKERRQ(ierr);
    for(i=0;i<n_var;++i){
      if(sol_ctx.dae_suffix_var->u.i[i]!=2 && sol_ctx.dae_suffix_var->u.i[i]!=3){
        xx[sol_ctx.dae_map_back[i]] = X0[i];
        xxl[sol_ctx.dae_map_back[i]] = LUv[i]; //lower bound
        xxu[sol_ctx.dae_map_back[i]] = Uvx[i]; //upper bound
      }
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(xl,&xxl);CHKERRQ(ierr);
    ierr = VecRestoreArray(xu,&xxu);CHKERRQ(ierr);
    // Print the variable types for reading trajectories
    strcpy(sol_ctx.opt.typ_file, sol_ctx.opt.stub);
    if(sol_ctx.opt.stublen > 3){
      if((sol_ctx.opt.stub[sol_ctx.opt.stublen - 1] == 'l') &&
         (sol_ctx.opt.stub[sol_ctx.opt.stublen - 2] == 'n') &&
         (sol_ctx.opt.stub[sol_ctx.opt.stublen - 3] == '.')){
        strcpy(sol_ctx.opt.typ_file + sol_ctx.opt.stublen - 3, ".typ");
      }
      else{ //just stub no ".nl" extension
        strcpy(sol_ctx.opt.typ_file + sol_ctx.opt.stublen, ".typ");
      }
    }
    else{// not long enough for there to be an extension
      strcpy(sol_ctx.opt.typ_file + sol_ctx.opt.stublen, ".typ");
    }

    err = PetscViewerASCIIOpen(PETSC_COMM_WORLD, sol_ctx.opt.typ_file, &pv);
    for(i=0;i<n_var;++i){
      PetscViewerASCIIPrintf(pv, "%d\n", sol_ctx.dae_suffix_var->u.i[i]);
    }
    err = PetscViewerDestroy(&pv);

    print_init_diagnostic(&sol_ctx); //print initial diagnostic information
    ierr = VecRestoreArray(x, &xx);CHKERRQ(ierr);
    // Make Jacobian matrix (by default sparse AIJ)
    ierr = MatCreate(PETSC_COMM_WORLD,&J); CHKERRQ(ierr);
    ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, n_con, sol_ctx.n_var_state); CHKERRQ(ierr);
    ierr = MatSetFromOptions(J); CHKERRQ(ierr); //command line options override defaults
    ierr = MatSetUp(J); CHKERRQ(ierr);
    /* finish up jacobian */
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /* Set residual and jacobian functions */
    ierr = TSSetIFunction(ts, r, FormDAEFunction, &sol_ctx);CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts, J, J, FormDAEJacobian, &sol_ctx);CHKERRQ(ierr);
    /* First set a bunch of default options, then read CL to override */
    ierr = TSSetProblemType(ts, TS_NONLINEAR);
    ierr = TSSetEquationType(ts, TS_EQ_IMPLICIT);
    ierr = TSGetSNES(ts, &snes);CHKERRQ(ierr);
    ierr = TSSetMaxSNESFailures(ts, DEFAULT_TS_SNES_MAX_FAIL);CHKERRQ(ierr);
    ierr = SNESGetKSP(snes, &ksp);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(snes, &linesearch);CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHBT);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,DEFAULT_PC);CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverType(pc, DEFAULT_LINEAR_PACK);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,DEFAULT_KSP);CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, 1);
    ierr = TSSetMaxTime(ts, 10);
    ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
    // Set up solver from CL options
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
    if(sol_ctx.opt.use_bounds) ierr = SNESVISetVariableBounds(snes, xl, xu);CHKERRQ(ierr);
    // Solve
    ierr = TSSolve(ts, x);
    ierr = TSGetTime(ts, &t);
    ierr = TSGetConvergedReason(ts, &tscr); CHKERRQ(ierr);
    /* Get the results */
    x_asl = (real*)malloc((n_var)*sizeof(real));
    ierr = VecGetArray(x, &xx);CHKERRQ(ierr);
    for(i=0;i<sol_ctx.n_var_state;++i) x_asl[sol_ctx.dae_map_x[i]] = xx[i];
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
    if(sol_ctx.explicit_time) x_asl[sol_ctx.dae_map_t] = t;
    /* write the AMPL solution file */
    get_ts_sol_message(msg, tscr, sol_ctx.asl);
    PetscPrintf(PETSC_COMM_SELF, "TSConvergedReason = %s\n", msg);
    write_sol(msg, x_asl, NULL, NULL); // write ASL sol file
    ierr = TSDestroy(&ts);
  } //end ts solve
  else{ // nonlinear solver setup and solve
    print_init_diagnostic(&sol_ctx); // Print requested diagnostic info
    /*Create nonlinear solver context*/
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);
    /*Create vectors for solution and nonlinear function*/
    ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr); //create an x vector
    ierr = VecSetSizes(x,PETSC_DECIDE,n_var);CHKERRQ(ierr); //vecs are n_vars long
    ierr = VecSetFromOptions(x);CHKERRQ(ierr); //command line options for vec
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);  // duplicate x for resuiduals
    ierr = VecDuplicate(x,&xl);CHKERRQ(ierr); // duplicate x for lower bounds
    ierr = VecDuplicate(x,&xu);CHKERRQ(ierr); // duplicate x for upper bounds
    //Make Jacobian matrix (by default sparse AIJ)
    ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
    ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n_con,n_var);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr); //command line options override defaults
    ierr = MatSetUp(J);CHKERRQ(ierr);
    /* finish up jacobian */
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /* Set residual and jacobian functions */
    ierr = SNESSetFunction(snes, r, FormFunction, &sol_ctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, J, J, FormJacobian, &sol_ctx);CHKERRQ(ierr);
    /* Default solver setup override from CL later*/
    ierr = SNESSetType(snes, DEFAULT_SNES);CHKERRQ(ierr);
    ierr = SNESGetKSP(snes, &ksp);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(snes, &linesearch);CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHBT);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,DEFAULT_PC);CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverType(pc, DEFAULT_LINEAR_PACK);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,DEFAULT_KSP);CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, DEFAULT_SNES_ATOL, DEFAULT_SNES_RTOL,
          0, DEFAULT_SNES_MAX_IT, DEFAULT_SNES_MAX_FUNC);CHKERRQ(ierr);
    ierr = KSPSetTolerances(
          ksp,DEFAULT_KSP_ATOL,DEFAULT_KSP_RTOL,PETSC_DEFAULT,200);CHKERRQ(ierr);
    /* Read command line options for solver */
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
    /* Make x vec set initial guess from nl file, also get lb and ub */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(xl,&xxl);CHKERRQ(ierr);
    ierr = VecGetArray(xu,&xxu);CHKERRQ(ierr);
    for(i=0;i<n_var;++i){
        xx[i] = X0[i]; //initial values
        xxl[i] = LUv[i]; //lower bound
        xxu[i] = Uvx[i]; //upper bound
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(xl,&xxl);CHKERRQ(ierr);
    ierr = VecRestoreArray(xu,&xxu);CHKERRQ(ierr);
    /* Set upper and lower bound most solver don't like but a few can use
       if you include bounds and solver can't use will cause failure */
    if(sol_ctx.opt.use_bounds) ierr = SNESVISetVariableBounds(snes, xl, xu);CHKERRQ(ierr);
    /* Solve it */
    ierr = SNESSolve(snes, NULL, x);CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snes, &cr); CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes, &its);CHKERRQ(ierr);
    /* Get the results */
    ierr = VecGetArray(x, &xx);CHKERRQ(ierr);
    /* write the AMPL solution file */
    get_snes_sol_message(msg, cr, sol_ctx.asl);
    PetscPrintf(PETSC_COMM_SELF, "SNESConvergedReason = %s, in %d iterations\n", msg, its);
    write_sol(msg, (real*)xx, NULL, NULL);
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
    ierr = VecDestroy(&xl);CHKERRQ(ierr);
    ierr = VecDestroy(&xu);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr); // ksp and pc are part of this
  } //end snes solve

  /* Should free stuff, but program ending anyway, so what's the point? */
  ASL_free(&asl);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = PetscFinalize(); //All done with PETSc
  return P_EXIT_NORMAL;
}

char **transform_args(int argc, char** argv, int *size){
  /* change the format of the command line arguments from ASL to PETSc*/
  const int max_args=200, max_arg_len=300;
  char argv2[max_args][max_arg_len]; // buffer for reformatted args
  char **argv3;
  int argc2 = argc; //reformatted number of args
  int i=0, j=0, k=0, h=0; // another counter

  for(i=0; i<argc; ++i){
    if(memcmp(argv[i], "-v", sizeof(char)*3) == 0){
      memcpy(argv2[i], "-version", sizeof(char)*9);
      ++k;
      continue;
    }
    for(j=0;argv[i][j]!='\0';++j){
      if (j==0){
        if (argv[i][j] == '-' && argv[i][j+1] == '-'){ //convert leading -- to -
          continue;
        }
      }
      if(argv[i][j]==' '||argv[i][j]=='\t'||argv[i][j]=='\n'){} //strip white spc
      else if(argv[i][j]=='='){ // split on '=''
        h=0;
        ++k;
      }
      else{
        argv2[k][h] = argv[i][j]; //keep anything else
        ++h;
      }
    }
    ++k; //on to next arg
    h=0;
  }
  argc2 = k;
  argv3 = (char**)malloc(sizeof(char*)*argc2);
  for(i=0; i<argc2; ++i){
    argv3[i] = (char*)malloc(sizeof(char)*(strlen(argv2[i])+2));
    memcpy(argv3[i], argv2[i], sizeof(char)*(strlen(argv2[i])+2));
  }
  *size = argc2;
  return argv3;
}

void sol_ctx_init(Solver_ctx *ctx){
  //Initialize some values in the solver context struct
  ctx->dae_map_t=-1;  //ASL index of time variable (-1 for none)
  ctx->n_var_diff=0; //DAE number of differential vars
  ctx->n_var_deriv=0; //DAE number of derivative vars
  ctx->n_var_state=0; //DAE number of state vars
  ctx->n_var_alg=0;  //DAE number of algebraic variables
  ctx->n_ineq=0;  // Number of inequality constraints
  ctx->explicit_time=0; //DAE includes time variable? 1=yes 0=no
  ctx->dof=0; //degrees of freedom
}

int ScaleEqsUser(Solver_ctx *sol_ctx){
    ASL *asl = sol_ctx->asl;
    int i = 0, err=0;
    real s;

    if (sol_ctx->scaling_factor_con->u.r == NULL) return 0; //no scaling factors provided
    for(i=0;i<n_con;++i){ //n_con is asl vodoo incase you wonder where it came from
      s = sol_ctx->scaling_factor_con->u.r[i];
      if(s != 0.0) conscale(i, s, &err);
    }
    return err;
}

int ScaleVarsUser(Solver_ctx *sol_ctx){
    //Use scaling factors set in the scaling_factor suffix, for DAEs ignore
    //scaling on the derivatives and use scaling from the differential vars
    //instead variables and there derivatives should be scaled the same
    int i=0, j=0, err=0;
    real s;
    ASL *asl = sol_ctx->asl;
    if (sol_ctx->scaling_factor_var->u.r == NULL) return 0; //no scaling factors
    if (sol_ctx->opt.dae_solve){ //dae variable scaling
      for(i=0;i<n_var;++i){ //n_var is asl vodoo
        if (sol_ctx->dae_suffix_var->u.i[i] == 2) {
          j = sol_ctx->dae_link[i]; //use differential var scale
          s = sol_ctx->scaling_factor_var->u.r[j];
          if(s != 0.0) varscale(i, 1.0/s, &err);
        }
        else if (sol_ctx->dae_suffix_var->u.i[i] == 3) {
          continue; //for now ignore time scaling
        }
        else {
          s = sol_ctx->scaling_factor_var->u.r[i];
          if(s != 0.0) varscale(i, 1.0/s, &err);
        }
      }
    }
    else{ // no dae so use scale factors given
      for(i=0;i<n_var;++i){ //n_var is asl vodoo
        s = sol_ctx->scaling_factor_var->u.r[i];
        if(s != 0.0) varscale(i, 1.0/s, &err); // invert scale to match Ipopt
      }
    }
    return err;
}

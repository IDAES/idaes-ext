#include"petsc.h"

void print_commandline(const char* msg, int argc, char **argv){
  /* print command line arguments */
  int i=0;
  PetscPrintf(PETSC_COMM_SELF, msg);
  for(i=0; i<argc; ++i) PetscPrintf(PETSC_COMM_SELF, "%s ", argv[i]);
  PetscPrintf(PETSC_COMM_SELF, "\n");
}

void print_init_diagnostic(Solver_ctx *sol_ctx){
  if(sol_ctx->opt.show_jac){
    print_jac_asl(sol_ctx->asl);
  }
  if(sol_ctx->opt.show_init){
    print_x_asl(sol_ctx->asl);
  }
  if(sol_ctx->opt.show_scale_factors){
    print_var_scale_factors_asl(sol_ctx->asl);
    print_con_scale_factors_asl(sol_ctx->asl);
  }
}

void print_x_asl(ASL *asl){
  int i=0;
  PetscPrintf(PETSC_COMM_SELF, "Initial Values (scaled)\n");
  for (i=0;i<n_var;++i){
     PetscPrintf(PETSC_COMM_SELF, "v%d: %e <= %e <= %e\n", i, LUv[i], X0[i], Uvx[i]);
  }
}

void print_var_scale_factors_asl(ASL *asl){
  int i;
  PetscPrintf(PETSC_COMM_SELF, "Variable Scale Factors:\n");
  if(asl->i.vscale!=NULL){
    for(i=0;i<n_var;++i) PetscPrintf(PETSC_COMM_SELF, "  v%d: %f\n", i, asl->i.vscale[i]);
  }
  else{
    PetscPrintf(PETSC_COMM_SELF, "  None\n");
  }
}

void print_con_scale_factors_asl(ASL *asl){
  int i;
  PetscPrintf(PETSC_COMM_SELF, "Constraint Scale Factors:\n");
  if(asl->i.cscale!=NULL){
    for(i=0;i<n_con;++i) PetscPrintf(PETSC_COMM_SELF, "  c%d: %f\n", i, asl->i.cscale[i]);
  }
  else{
    PetscPrintf(PETSC_COMM_SELF, "  None\n");
  }
}

void print_jac_asl(ASL *asl){
  /* Print sparse Jacobian  highlight elements over u or under l*/
  cgrad          *cg;  /* sparse jacobian elements*/
  real           *Jac; /* ASL test Jacobian */
  int            err;  /* Error code  from ASL fulctions */
  int            i=0;

  Jac = (real *)Malloc(nzc*sizeof(real)); /* Jacobian space */
  jacval(X0, Jac, &err); /*calculate jacobian */
  PetscPrintf(PETSC_COMM_SELF, "Computed Jacobian, err = %d\n", err);
  PetscPrintf(PETSC_COMM_SELF, "Computed Jacobian values (scaled):\n");
  for(i=n_conjac[0];i<n_conjac[1]; ++i){ /*i is constraint index */
    cg = Cgrad[i];
    PetscPrintf(PETSC_COMM_SELF, "c%d", i);
    while(cg!=NULL){
      PetscPrintf(PETSC_COMM_SELF, " v%d(%e)", cg->varno, Jac[cg->goff]);
      cg=cg->next;
    }
    PetscPrintf(PETSC_COMM_SELF, "\n");
  }
  free(Jac);
}


/* Help message prepended to the regular PETSc help */

static char help[] ="\
--------------------------------------------------------------------------\n\
Use ASL and PETSc to solve a problem (non-linear equations or DAE) defined \n\
in an AMPL nl file. Optimization solver support will be added soon.\n\n\
   Added Options: \n\
     [filename]: File name with or without extension (when -s is not specified) \n\
     -s <stub>: File name with or without extension\n\
     -show_scale_factors: Show the calculated or user specified scale factors\n\
     -show_jac: Show non-zero jacobian values at initial point\n\
     -show_intial: Show the guess intial values \n\
     -show_cl: Show the command line input and transformation from AMPL format\n\
     -perturb_test <factor>: Test using nl file where initial value is solution \n\
        This just multiplies the inital value by factor and resolves, then it \n\
        compares the new solution to the intial value and report differnces > \n\
        1e-6. It's a fairly crude test but it should ensure the Jacobian and \n\
        function are calculating right.\n\
     -dae_solve: Run DAE solver, must provide appropriate suffixes\n\
     -jac_explicit_diag: Create explicit Jacobian entries in diagonal\n\
    ";

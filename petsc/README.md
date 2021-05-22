# README

This is a first pass at creating an AMPL executable for the PETSc solver library, allowing PETSc to be used with Pyomo.  This version supports the TS (time-stepping) and SNES (nonlinear) solvers; however in the future TAO (optimization) solvers and additional features may also be supported.

## Usage

The petsc executable solver will run like any other AMPL executable and can take the standard PETSc options.  The options can be set from within Pyomo.  For use with Pyomo it is easiest if the petsc executable is in your executable path.

<TODO> Will fill out later, for now see examples and PETSc docs.  PETSc command line options work with the PETSc solvers. The are some additional command line arguments (see ```petsc -help```).  Command line arguments can be passed through Pyomo's solver interface.

## Testing

The included tl.nl file can be used to test the petsc solver. The problem is an old version of the IDAES MEA model, but that's not important. The initial values of variables in the file are the solution.  To test the solver the initial values can be perturbed before solving the problem the new solution can be compared to the old initial values to ensure the problem solved.

To test the petsc executable run:

petsc -s t1 -snes_monitor -perturb_test 1.1

This multiplies the initial values by 1.1 and resolves.  The results shows any differences between the new and old solutions greater than 1e-6, and the values of those variables.  Depending on the value of the variables, difference of more than 1e-6 are not necessarily bad.

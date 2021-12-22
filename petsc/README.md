# README

This is a first pass at creating an AMPL executable for the PETSc solver library,
allowing PETSc to be used with Pyomo via the AMPL Solver Library (ALS).  This
version supports the TS (time-stepping) and SNES (nonlinear) solvers; however in
the future TAO (optimization) solvers and additional features may also be
supported.

## Building

The solver wrapper build uses the PETSc build system.  On Linux, compile PETSc
with the desired features (https://petsc.org/release/install/). The ASL can be
obtained from http://www.netlib.org/ampl/solvers/. Before building the wrapper
set the following environment variables:

* ASL_INC - directory with ASL header files
* ASL_LIB - ASL library file
* PETSC_DIR - see PETSc docs
* PETSC_ARCH - see PETSc docs

## Usage

The SNES solvers can be used with Pyomo to solve sets of nonlinear equations
as a standard AMPL solver. Tools to use the time stepping solver with Pyomo.DAE
are under development as part of the IDAES project, and example can be found
there (https://github.com/IDAES/idaes-pse/pull/552).

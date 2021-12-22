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

// AMPL solver interface for PETSc, help for additional cl args
// Author: John Eslick

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
     -dae_solve: Run DAE solver, must provide appropriate suffixes\n\
     -ignore_scaling: ignore scaling suffixes\n\
    ";

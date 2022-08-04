/*-------------------------------------------------------------------------------+
| The Institute for the Design of Advanced Energy Systems Integrated Platform    |
| Framework (IDAES IP) was produced under the DOE Institute for the              |
| Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021      |
| by the software owners: The Regents of the University of California, through   |
| Lawrence Berkeley National Laboratory,  National Technology & Engineering      |
| Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University |
| Research Corporation, et al.  All rights reserved.                             |
|                                                                                |
| Please see the files COPYRIGHT.md and LICENSE.md for full copyright and        |
| license information.                                                           |
+-------------------------------------------------------------------------------*/

#include<vector>
#include <cstddef>


#ifndef _INCLUDE_SOLVER_H_
#define _INCLUDE_SOLVER_H_

/*--------------------------------------------------------------------------------
  Function pointer types for generic solver implimentations
--------------------------------------------------------------------------------*/

typedef void (*fgh1_ptr)(double, std::vector<double>*, void*); // 1 var, gad, hes
typedef double (*f1_ptr)(double, void*); // 1 var

/*-------------------------------------------------------------------------------
  General purpose solvers for to f(x) = 0
-------------------------------------------------------------------------------*/

// Bracketing solvers bisection and false-position, there are used when the
// solution is in a know range, but newton solvers are not dependable.  These
// can also be used to get a good initial guess.
int bracket(
  f1_ptr f,
  double xa,
  double xb,
  double *sol,
  int max_it,
  double ftol,
  double xtol,
  void *data = NULL
);

// Halley is a second order newton method (needs 2nd derivatives)
int halley(
  fgh1_ptr f,
  double x0,
  double *sol,
  std::vector<double> *fg,
  int max_it,
  double ftol,
  void *data = NULL
);

// 1D Newton with optional line search
int newton_ls(
  fgh1_ptr f,
  f1_ptr fls,
  double x0,
  double *sol,
  std::vector<double> *fg,
  int max_it,
  double ftol,
  void *data = NULL,
  bool lsearch=1,
  double c=0.001,
  double t=0.5
);

#endif

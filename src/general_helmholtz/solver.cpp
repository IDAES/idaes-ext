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

/*--------------------------------------------------------------------------------
 Specific water functions from:

 Generic 1D solver implimentations:
     1) False-Position,
     2) Halley (2nd order Newton), and
     3) Newton with line search

 Flase-position can be used when you know the range where the solution is, but
 need a closer starting point to get the Newton method to converge. At this
 point, Halley's method is used for the Newton solver in most places since we
 are already calcuating second derivatives.  The Newton lineseach could be
 useful when starting from less good initial points, and may speed things up,
 but more testing is required.

 Author: John Eslick
 File: solver.cpp
--------------------------------------------------------------------------------*/

#include "solver.h"
#include<cmath>


/*------------------------------------------------------------------------------
 False-position type bracketing solver.  This is mostly used to confine the
 search for a good starting point to a specific area, where a newton method
 given a starting point too far way may not converge or converge to the wrong
 solution. It may be a bit slow, but this can also be used for solve
 a problem competely.
------------------------------------------------------------------------------*/
int bracket(
  f1_ptr f,
  double xa,
  double xb,
  double *sol,
  int max_it,
  double ftol,
  double xtol,
  void *data
){
  double fa, fb, fc, a=xa, b=xb, c;
  int it;
  bool prev_a=0, prev_b=0;

  fa = f(a, data);
  fb = f(b, data);

  if (fa*fb > 0){ //no solution (or multiple solutions)
    if (fabs(fa) < fabs(fb)){
      *sol = xa;
      return -1; // no solution xa is closest
    }
    *sol = xb;
    return -2; // no solution xb is closest
  }

  for(it=0; it<max_it; ++it){
    c = b - fb*(b - a)/(fb - fa);
    fc = f(c, data);
    if(fc*fa >= 0){
      a = c;
      fa = fc;
      if (prev_a){ //Illinois mod
        fb *= 0.5;
      }
      prev_a = 1;
      prev_b = 0;
    }
    else{
      b = c;
      fb = fc;
      if (prev_b){ //Illinois mod
        fa *= 0.5;
      }
      prev_a = 0;
      prev_b = 1;
    }
    if (fabs(fa) < ftol && !prev_b) { // Sometimes one side will converge to solution
      *sol = a;
      return it;
    }
    if (fabs(fb) < ftol && !prev_a) { // Sometimes one side will converge to solution
      *sol = b;
      return it;
    }
    *sol = (a + b)/2.0;
    if(fabs(b - a) < xtol) {
      return it;
    }
  }
  return -3; // max iterations sol is up to date here
}

/*------------------------------------------------------------------------------
  Halley's method
------------------------------------------------------------------------------*/
int halley(
  fgh1_ptr f,
  double x0,
  double *sol,
  std::vector<double> *fg,
  int max_it,
  double ftol,
  void *data
){
  int it = 0;
  double fun_prev = 0, x = x0;
  f(x0, fg, data);
  while(fabs((*fg)[0]) > ftol && it < max_it && (*fg)[0] != fun_prev){
    x -= (*fg)[0] * (*fg)[1] / ((*fg)[1]*(*fg)[1] - 0.5 * (*fg)[0] * (*fg)[2]);
    fun_prev = (*fg)[0];
    f(x, fg, data);
    ++it;
  }
  *sol = x;
  if(it == max_it){
    return -3;
  }
  return it;
}

/*------------------------------------------------------------------------------
  1D Newton's method with backtracking line search.
------------------------------------------------------------------------------*/
int newton_ls(
  fgh1_ptr f,
  f1_ptr fls,
  double x0,
  double *sol,
  std::vector<double> *fg,
  int max_it,
  double ftol,
  void *data,
  bool lsearch,
  double c,
  double t
){
  int it=0, it2=0;
  const int max_ls_it = 20;
  double x = x0, alpha, p, funx, fun;
  f(x0, fg, data);

  while(fabs((*fg)[0]) > ftol && it < max_it){
    alpha = -(*fg)[0]/(*fg)[1];
    if (lsearch){
      p = (alpha < 0) ? -1 : 1;
      funx = (*fg)[0];
      fun = fls(x + alpha, data);
      it2 = 0;  // line search iterations
      //  The while loop cuts step if
      //    1) improvment is not enough,
      //    2) function is NaN, or
      //    3) function is inf.
      //  The c parameter determines what is considered enough improvment to
      //  accept the step and t determines how much to cut the step.
      while(
        (fabs(funx) - fabs(fun) < c*p*(*fg)[1] || std::isnan(fun) || std::isinf(fun)) &&
        it2 < max_ls_it
      ){
        alpha *= t;
        fun = fls(x + alpha, data);
        ++it2;
      }
    }
    x += alpha;
    f(x, fg, data);
    ++it;
  }
  *sol = x;
  if(it == max_it){
    return -3;
  }
  return it;
}

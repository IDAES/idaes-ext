#ifndef _INCLUDE_HELMHOLTZ_SOLVE_H_
#define _INCLUDE_HELMHOLTZ_SOLVE_H_

#include "helmholtz_config.h"


// These are function pointer types for generic solver implimentations here,
// one for functions of one variable and one for functions of two variables,
// The pointer arguments are for the gradient and hessian.

typedef s_real (*f_ptr1)(s_real, s_real*, s_real*);
typedef s_real (*f_ptr2)(s_real, s_real, s_real*, s_real*);
typedef s_real (*f_ptr3)(s_real, s_real, s_real, s_real*, s_real*);
typedef s_real (*f_ptr2n)(s_real, s_real); //two vars no derivatives

// For a lot of solves I want to use binary functions where on of the two
// arguments is fixed.  To avoid implimneting the solvers three times for
// f(x), f(a, x), f(x, a), I'll use a functor to wrap up functions for solves
class FuncWrapper {
public:
  FuncWrapper(char apos, s_real a, s_real c);
  s_real operator() (s_real, s_real *grad, s_real *hes);
  s_real operator() (s_real, s_real, s_real *grad, s_real *hes);
  s_real operator() (s_real, s_real);
  s_real operator() (s_real);
  int grad_pos;
  int hes_pos;
  void set_f1(f_ptr1 f1);
  void set_f2(f_ptr2 f2);
  void set_f3(f_ptr3 f3);
  void set_f2n(f_ptr2n f2);
private:
  s_real _a;
  s_real _c;
  unsigned char _apos;
  f_ptr1 _f1;
  f_ptr2 _f2;
  f_ptr3 _f3;
  f_ptr2n _f2n;
};

/*------------------------------------------------------------------------------
  Geranl purpose 1 and 2 dimensional solvers. These solvers find the roots
  f(x) - c = 0
------------------------------------------------------------------------------*/

int bracket(FuncWrapper *f, s_real xa, s_real xb,
            s_real *sol, int max_it, s_real ftol, s_real xtol);
int halley(FuncWrapper *f, s_real x0, s_real *sol, s_real *grad, s_real *hes,
           int max_it, s_real ftol);
int newton_1d(FuncWrapper *f, s_real x0, s_real *sol, s_real *grad, s_real *hes,
              int max_it, s_real ftol);
int newton_2d(FuncWrapper *f0, FuncWrapper *f1, s_real x00, s_real x10,
              s_real *grad0, s_real *hes0, s_real *grad1, s_real *hes1,
              s_real *sol0, s_real *sol1, int max_it, s_real ftol);

/*------------------------------------------------------------------------------
  Function to calculate delta (rho/rho_c) from P and tau (T_c/T)
------------------------------------------------------------------------------*/

s_real delta_p_tau_rf(s_real pr, s_real tau, s_real a, s_real b);
s_real delta_vap(s_real p, s_real tau, s_real *grad=NULL, s_real *hes=NULL);
s_real delta_liq(s_real p, s_real tau, s_real *grad=NULL, s_real *hes=NULL);

/*------------------------------------------------------------------------------
  Functions for saturation pressure and density as a function to tau (T_c/T).
  These use the method of Akasaka (2008), which is good at least up to 647.093 K
  In this case the last bit of the density curves from 647.090K to 647.096K was
  filled in with cubics that match density, first and second derivatives at the
  transition point and critical density at the critical point.
------------------------------------------------------------------------------*/
s_real sat_delta_liq(s_real tau); // saturated liquid density at tau
s_real sat_delta_vap(s_real tau); // saturated vapor density at tau
s_real sat_p(s_real tau);         // saturation pressure at tau
int sat(s_real tau, s_real *delta_l_sol, s_real *delta_v_sol); //sat solver

s_real sat_tau_with_derivs(s_real pr, s_real *grad, s_real *hes);
s_real tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes);
s_real mem_tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes);
s_real tau_from_sp_with_derivs(s_real st, s_real pr, s_real *grad, s_real *hes);
s_real tau_from_up_with_derivs(s_real ut, s_real pr, s_real *grad, s_real *hes);

s_real p_from_htau_with_derivs(s_real ht, s_real tau, s_real *grad, s_real *hes);
s_real p_from_stau_with_derivs(s_real st, s_real tau, s_real *grad, s_real *hes);
s_real p_from_utau_with_derivs(s_real ut, s_real tau, s_real *grad, s_real *hes);



#endif

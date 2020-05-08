#ifndef _INCLUDE_HELMHOLTZ_SOLVE_H_
#define _INCLUDE_HELMHOLTZ_SOLVE_H_

#include "helmholtz_config.h"


/*------------------------------------------------------------------------------
  Function to calculate delta (rho/rho_c) from P and tau (T_c/T)
------------------------------------------------------------------------------*/

// Solve for delta using Halley's method delta0 in an inital guess, tol is
// the absolute residual tolerance, nit provides a means to return number of
// iterations for testing there is more than 1 solution maybe even more than one
// physically meaningfull solution, so delta0 is very important
s_real delta_p_tau(s_real p, s_real tau, s_real delta_0, s_real tol=1e-10,
                   int *nit=NULL, s_real *grad=NULL, s_real *hes=NULL);

s_real delta_p_tau_rf(s_real pr, s_real tau, s_real a, s_real b, bool bisect);

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

s_real sat_tau_with_derivs(s_real pr, s_real *grad, s_real *hes, int *nit=NULL);
s_real tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes);
s_real tau_from_sp_with_derivs(s_real st, s_real pr, s_real *grad, s_real *hes);
s_real tau_from_up_with_derivs(s_real st, s_real pr, s_real *grad, s_real *hes);

#endif

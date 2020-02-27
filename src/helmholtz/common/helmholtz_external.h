/*------------------------------------------------------------------------------
 Institute for the Design of Advanced Energy Systems Process Systems
 Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
 software owners: The Regents of the University of California, through
 Lawrence Berkeley National Laboratory,  National Technology & Engineering
 Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
 University Research Corporation, et al. All rights reserved.

 Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
 license information, respectively. Both files are also available online
 at the URL "https://github.com/IDAES/idaes".
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 Main IAPWS R6-95(2016) steam calculations. For now, the non-analytic terms are
 not included, because they cause a singularity at the critical point.  It is
 assumed that we will generally not be operating very near the critical point,
 so well behaved calculations are prefered to high accuracy in near the critical
 point.

 For references see iapws95.h

 Author: John Eslick
 File iapws95_external.h
------------------------------------------------------------------------------*/

#include "helmholtz_config.h"
#include "helmholtz_phi.h"
#include "helmholtz_solve.h"

#ifndef _INCLUDE_HELMHOLTZ_EXTERNAL_H_
#define _INCLUDE_HELMHOLTZ_EXTERNAL_H_

/*------------------------------------------------------------------------------
  Basic thermodynamic quantities as functions of delta (rho/rho_c) and
  tau (T_c/T)
------------------------------------------------------------------------------*/
s_real p(s_real delta, s_real tau);
s_real u(s_real delta, s_real tau);
s_real s(s_real delta, s_real tau);
s_real h(s_real delta, s_real tau);
s_real f(s_real delta, s_real tau);
s_real cv(s_real delta, s_real tau);
s_real cp(s_real delta, s_real tau);
s_real g(s_real delta, s_real tau);
s_real w(s_real delta, s_real tau);

/*------------------------------------------------------------------------------
  Basic thermodynamic functions with gradient and Hessian as functions of
  delta (rho/rho_c) and tau (T_c/T)
------------------------------------------------------------------------------*/
s_real p_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real u_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real h_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real s_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real f_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real g_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real cp_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real cv_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);
s_real w_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes);

s_real hvpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes);
s_real hlpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes);
s_real vf_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes);

/*------------------------------------------------------------------------------
  Functions for saturation pressure and reduced density as a function of tau
  (T_c/T) with gradient and Hessian.
------------------------------------------------------------------------------*/
s_real sat_delta_liq_with_derivs(s_real tau, s_real *grad, s_real *hes);
s_real sat_delta_vap_with_derivs(s_real tau, s_real *grad, s_real *hes);
s_real sat_p_with_derivs(s_real tau, s_real *grad, s_real *hes, bool limit=1);

/*------------------------------------------------------------------------------
  Functions for saturation delta (rho/roh_c) as a function of P and tau (T_c/T)
  with gradient and Hessian.  These functions are well tested and just feed the
  right initial guess to delta_p_tau() gauranteeing the correct solution.
------------------------------------------------------------------------------*/
s_real delta_vap(s_real p, s_real tau, s_real *grad=NULL, s_real *hes=NULL, int *nit=NULL);
s_real delta_liq(s_real p, s_real tau, s_real *grad=NULL, s_real *hes=NULL, int *nit=NULL);
#endif

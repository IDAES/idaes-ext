/*------------------------------------------------------------------------------
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 Main IAPWS R6-95(2016) steam calculations. This file contains functions to
 solve for liquid and vapor denstiy from temperature and pressure.  It also
 contains functions for solving for saturation conditions.

 For references see iapws95.h

 Author: John Eslick
 File iapws95_solve.cpp
------------------------------------------------------------------------------*/

#include<stdio.h>
#include<cmath>
#include<iostream>

#include"iapws95_memo.h"
#include"iapws95_config.h"
#include"iapws95_deriv_parts.h"
#include "iapws95_solve.h"

/*------------------------------------------------------------------------------
In this section you will find function for calculating density from pressure
and teperature. This lets you calculate properties in the more standard way as
a function of temperature and pressure instead of density and pressure.

Unfortunatly this is difficult and a bit messy.
*-----------------------------------------------------------------------------*/
s_real delta_p_tau_rf(s_real pr, s_real tau, s_real a, s_real b, bool bisect){
  /*----------------------------------------------------------------------------
  Bracketing methods, false position and bisection, for finding a better initial
  guess for density when solving for density from temperature and pressure. This
  is only used in particularly diffucult areas.  At this point it probably
  overused, but there are places where it is probably necessary.

  Args:
    pr: pressure (kPa)
    tau: inverse of reduced pressure Tc/T
    bisect: 1 to use bisection (probably slow), 0 to use false position
            bisection isn't really used, but I kept it around for debugging
    a: first density bound (kg/m3)  | really density, why didn't use delta?
    b: second density bound (kg/m3) | it was easier to think in terms of density
  Returns:
    delta: reduced density at pr and tau (more or less approximate)
  ----------------------------------------------------------------------------*/
  s_real c=a, fa, fb, fc;
  int it = 0;
  // If right by critical point guess critical density. (okay this isn't part
  // of a backeting method but it was conveneint).
  if( fabs(T_c/tau - T_c) < 1e-7 && fabs(pr - P_c) < 1e-4) return 1;
  // solve f(delta, tau) = 0; f(delta, tau) = p(delta, tau) - pr
  a /= rho_c; // convert to reduced density
  b /= rho_c; // convert to reduced density
  fa = p(a, tau) - pr; // initial f(a, tau)
  fb = p(b, tau) - pr; // initial f(b, tau)
  while(it < MAX_IT_BRACKET && (a - b)*(a - b) > TOL_BRACKET){
    if(bisect) c = (a + b)/2.0; // bisection
    else c = b - fb*(b - a)/(fb - fa); //regula falsi
    fc = p(c, tau) - pr; // calcualte f(c)
    if(fc*fa >= 0){a = c; fa = fc;}
    else{b = c; fb = fc;}
    ++it;
  }
  return (a+b)/2.0;
}

s_real delta_p_tau_liq_guess(s_real p, s_real tau){
  s_real delta0, a=265, b=435;
  bool use_rf = 0;
  if(tau <= 1.0 || p >= P_c){
    //over critical T or P
    s_real T = T_c/tau;
    if(p >= 201.130822*T - 108125.908560 && p <= 513.615856*T - 308340.302065){
       //Supercritical and above 230 kg/m3 below and 500 kg/m3
      use_rf = 1; a = 210; b = 515;
    }
    else if(p >= -0.030275366391367*T*T + 800.362108600627*T - 468414.629148942){//rho >= 600
      delta0 = 1100.0/rho_c;}
    else if(p >= 0.456576012809866*T*T - 224.822653977863*T - 23457.7703540548){ //rho >= 425
      delta0 = 500/rho_c;}
    else if(p <= -2.26154966031622E-06*T*T + 0.467571780470989*T - 4.3044421839477){//rho <= 1
      delta0 = 0.5/rho_c;}
    else if(p <= -0.001481257848455*T*T + 15.4473970626314*T - 2739.92167514421){//rho <= 25
      delta0 = 10.0/rho_c;}
    else {
      a = 10;
      b = 450;
      use_rf = 1;
    }
  }
  else if(p > 21.5 && tau < T_c/645){
    use_rf = 1;
    a = 322;
    b = 480;
  }
  else delta0=1100/rho_c;
  if(use_rf) delta0 = delta_p_tau_rf(p, tau, a, b, 0); //bracket for better i.g.
  return delta0;
}

s_real delta_p_tau(s_real pr, s_real tau, s_real delta_0, s_real tol, int *nit,
  s_real *grad, s_real *hes){
  /*----------------------------------------------------------------------------
  Halley's method to calculate density from temperature and pressure

  Args:
   pr: pressure (kPa)
   tau: inverse of reduced pressure Tc/T
   delta0: initial guess for delta
   tol: absolute residual tolerance for convergence
   nit: pointer to return number of iterations to, or NULL
   grad: location to return gradient of delta wrt pr and tau or NULL
   hes: location to return hessian (upper triangle) or NULL
  Returns:
   delta: reduced density (accuracy depends on tolarance and function shape)
  ----------------------------------------------------------------------------*/
  s_real delta = delta_0, fun, gradp[2], hesp[3];
  int it = 0; // iteration count
  fun = p_with_derivs(delta, tau, gradp, hesp) - pr;
  while(fabs(fun) > tol && it < MAX_IT_DELTA){
    delta = delta - fun*gradp[0]/(gradp[0]*gradp[0] - 0.5*fun*hesp[0]);
    fun = p_with_derivs(delta, tau, gradp, hesp) - pr;
    ++it;
  }
  if(nit != NULL) *nit = it;
  if(grad != NULL){ // calculate gradient if needed
    grad[0] = 1.0/gradp[0];
    grad[1] = -gradp[1]*grad[0]; //triple product
    if(hes != NULL){ // calculate hession if needed.
      hes[0] = -hesp[0]*grad[0]/gradp[0]/gradp[0];
      hes[1] = -(hesp[1] + hesp[0]*grad[1])/gradp[0]/gradp[0];
      hes[2] = -(grad[0]*(hesp[2] + grad[1]*hesp[1]) + gradp[1]*hes[1]);}}
  return delta;
}

s_real delta_liq(s_real p, s_real tau, s_real *grad, s_real *hes, int *nit){
  /*----------------------------------------------------------------------------
  Get a good liquid or super critical inital density guess then call
  delta_p_tau() to calculate density. In difficult cases an interval with the
  desired soultion is used with a backeting method to produce a better inital
  guess.  There is a area around the critical point and in the super critical
  region where this is hard to solve so there is some pretty extensive code for
  guessing.

  Since solving this takes some time and this function is called a lot with the
  exact same inputs when using the Pyomo wrapper, this function is memoized.

  Args:
   p: pressure (kPa)
   tau: inverse of reduced pressure Tc/T
   grad: location to return gradient of delta wrt pr and tau or NULL
   hes: location to return hessian (upper triangle) of NULL
   nit: pointer to return number of iterations to, or NULL
  Returns:
   delta: reduced density (accuracy depends on tolarance and function shape)
  ----------------------------------------------------------------------------*/
  s_real val = memoize::get_bin(memoize::delta_liq, p, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  s_real delta, delta0, tol=TOL_DELTA_LIQ;
  delta0 = delta_p_tau_liq_guess(p, tau);
  delta = delta_p_tau(p, tau, delta0, tol, nit, grad, hes); //solve
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){
    // This is just to avoid evaluation errors.  Want to be able to calucalte
    // vapor properties even when vapor doesn't exist.  In the IDAES Framework
    // these properties may be calculated and multipled by a zero liquid fraction,
    // so it doesn't mater that they are wrong.
    delta = 3.1;
    zero_derivs2(grad, hes);
  }
  memoize::add_bin(memoize::delta_liq, p, tau, delta, grad, hes); //store
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   //   function
  return delta;
}

s_real delta_vap(s_real p, s_real tau, s_real *grad, s_real *hes, int *nit){
  /*----------------------------------------------------------------------------
  Get a good vapor or super critical inital density guess then call
  delta_p_tau() to calculate density. In the supercritical region this just
  calls the liquid function. In the rest of the vapor region the inital guess
  is pretty easy.

  Since solving this takes some time and this function is called a lot with the
  exact same inputs when using the Pyomo wrapper, this function is memoized.

  Args:
   p: pressure (kPa)
   tau: inverse of reduced pressure Tc/T
   grad: location to return gradient of delta wrt pr and tau or NULL
   hes: location to return hessian (upper triangle) of NULL
   nit: pointer to return number of iterations to, or NULL
  Returns:
   delta: reduced density (accuracy depends on tolarance and function shape)
  ----------------------------------------------------------------------------*/
  s_real val = memoize::get_bin(memoize::DV_FUNC, p, tau, grad, hes);
  if(!std::isnan(val)) return val; // return stored result if available
  s_real delta, delta0 = 0.01;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // If supercritical use the liquid calculation, which includes sc region
  if(tau <= 1.0 || p >= P_c) return delta_liq(p, tau, grad, hes);
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  delta0 = 0.0001; // this guess should work for vapor region
  delta = delta_p_tau(p, tau, delta0, TOL_DELTA_VAP, nit, grad, hes);
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){
    // This is just to avoid evaluation errors.  Want to be able to calucalte
    // vapor properties even when vapor doesn't exist.  In the IDAES Framework
    // these properties may be calculated and multipled by a zero vapor fraction,
    // so it doesn't mater that they are wrong.
    delta = 0.001;
    zero_derivs2(grad, hes);
  }
  memoize::add_bin(memoize::DV_FUNC, p, tau, delta, grad, hes); // store result
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return delta;
}


/*------------------------------------------------------------------------------
In this section you will find functions for calculating the saturation curve.
Staturation pressure and density as a function of temperature.
*-----------------------------------------------------------------------------*/
s_real sat_delta_liq(s_real tau){ //caculate saturated liquid density from tau
  s_real delta_l, delta_v;
  sat(tau, &delta_l, &delta_v);
  return delta_l;
}

s_real sat_delta_vap(s_real tau){ //caculate saturated vapor density from tau
  s_real delta_l, delta_v;
  sat(tau, &delta_l, &delta_v);
  return delta_v;
}

s_real p_sat_iapws97(s_real tau){ //saturation pressure from tau IAPWS-97 eq.
  //the IAPWS-97 isn't as consitent as IAPWS-95, but this provides a good guess
  s_real T = T_c/tau;
  s_real tt = T + n_psat[8]/(T - n_psat[9]);
  s_real A = tt*tt + n_psat[0]*tt + n_psat[1];
  s_real B = n_psat[2]*tt*tt + n_psat[3]*tt + n_psat[4];
  s_real C = n_psat[5]*tt*tt + n_psat[6]*tt + n_psat[7];
  return 1000*pow(2*C/(-B + pow(B*B - 4*A*C, 0.5)), 4);
}

s_real delta_sat_v_approx(s_real tau){ //approximate saturated vapor density
  // This equation is from the original IAPWS-95 paper
  s_real XX = 1 - 1.0/tau;
  s_real delta = exp(-2.03150240*pow(XX,2.0/6.0)
            - 2.68302940*pow(XX,4.0/6.0)
            - 5.38626492*pow(XX,8.0/6.0)
            - 17.2991605*pow(XX,18.0/6.0)
            - 44.7586581*pow(XX,37.0/6.0)
            - 63.9201063*pow(XX,71.0/6.0));
  return delta_p_tau(p_sat_iapws97(tau), tau, delta);
}

s_real delta_sat_l_approx(s_real tau){ //approximate saturated vapor density
  // This equation is from the original IAPWS-95 paper.
  s_real XX = 1 - 1.0/tau;
  s_real delta = 1.001
           + 1.99274064*pow(XX,1.0/3.0)
           + 1.09965342*pow(XX,2.0/3.0)
           - 0.510839303*pow(XX,5.0/3.0)
           - 1.75493479*pow(XX,16.0/3.0)
           - 45.5170352*pow(XX,43.0/3.0)
           - 6.74694450e5*pow(XX,110.0/3.0);
  return delta_p_tau(p_sat_iapws97(tau), tau, delta);
}

inline s_real J(s_real delta, s_real tau){
  // Term from Akasaka method for saturation state
  return delta*(1+delta*phir_delta(delta, tau));
}

inline s_real K(s_real delta, s_real tau){
  // Term from Akasaka method for saturation state
  return delta*phir_delta(delta,tau) + phir(delta, tau) + log(delta);
}

inline s_real J_delta(s_real delta, s_real tau){
  return 1.0 + 2.0*delta*phir_delta(delta, tau) + delta*delta*phir_delta2(delta, tau);
}

inline s_real K_delta(s_real delta, s_real tau){
  return 2.0*phir_delta(delta, tau) + delta*phir_delta2(delta, tau) + 1.0/delta;
}

inline s_real Delta_Aka(s_real delta_l, s_real delta_v, s_real tau){
  return J_delta(delta_v, tau)*K_delta(delta_l, tau) -
         J_delta(delta_l, tau)*K_delta(delta_v, tau);
}

int sat(s_real tau, s_real *delta_l_sol, s_real *delta_v_sol){
    //Get stautated phase densities at tau by Akasaka (2008) method
    s_real delta_l, delta_v, fg, gradl[1], hesl[1], gradv[1], hesv[1];
    int n=0, max_it=MAX_IT_SAT;

    if(tau - 1 < 1e-12){
      delta_l = 1.0;
      delta_v = 1.0;
      max_it=0;
    }
    else{
      // okay so you've decided to solve this thing
      delta_l = delta_sat_l_approx(tau); // guess based on IAPWS-97
      delta_v = delta_sat_v_approx(tau); // guess based on IAPWS-97
    }
    // Since the equilibrium conditions are gl = gv and pl = pv, I am using the
    // the relative differnce in g as a convergence criteria, that is easy to
    // understand.  fg < tol for convergence, fg is calucalted upfront in the
    // off chance that the guess is the solution
    *delta_l_sol = delta_l; // just in case we don't do at least 1 iteration
    *delta_v_sol = delta_v; // just in case we don't do at least 1 iteration
    fg = fabs((g(delta_v, tau) - g(delta_l, tau))/g(delta_l, tau));
    while(n<max_it && fg > TOL_REL_SAT_G){
      ++n; // Count iterations
      //calculations deltas at next step (Akasaka (2008))
      *delta_l_sol = delta_l + SAT_GAMMA/Delta_Aka(delta_l, delta_v, tau)*(
             (K(delta_v, tau) - K(delta_l, tau))*J_delta(delta_v,tau) -
             (J(delta_v, tau) - J(delta_l, tau))*K_delta(delta_v,tau));
      *delta_v_sol = delta_v + SAT_GAMMA/Delta_Aka(delta_l, delta_v, tau)*(
             (K(delta_v, tau) - K(delta_l, tau))*J_delta(delta_l,tau) -
             (J(delta_v, tau) - J(delta_l, tau))*K_delta(delta_l,tau));
      delta_v = *delta_v_sol; //step
      delta_l = *delta_l_sol;
      //calculate convergence criterium
      fg = fabs((g(delta_v, tau) - g(delta_l, tau))/g(delta_l, tau));
    }

    //Calculate grad and hes for and memoize

    gradv[0] = LHM/LGM;
    gradl[0] = gradv[0]*LBV/LBL + (LCV - LCL)/LBL;
    hesv[0] = LdHdt(delta_l, delta_v, tau, gradl[0], gradv[0])/LGM
             - LHM/LGM/LGM*LdGdt(delta_l, delta_v, tau, gradl[0], gradv[0]);
    hesl[0] = hesv[0]*LBV*LFL + gradv[0]*(LBVt + LBVd*gradv[0])*LFL
              + gradv[0]*LBV*(LFLt + LFLd*gradl[0]) + (LFLt + LFLd*gradl[0])*(LCV - LCL)
              + LFL*(LCVt - LCLt + LCVd*gradv[0] - LCLd*gradl[0]);

    memoize::add_un(memoize::DL_SAT_FUNC, tau, delta_l, gradl, hesl);
    memoize::add_un(memoize::DV_SAT_FUNC, tau, delta_v, gradv, hesv);
    return n;
}

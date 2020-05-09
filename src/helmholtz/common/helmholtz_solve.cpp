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
 This file contains functions to solve for liquid and vapor density from
 temperature and pressure.  It also contains functions for solving for
 saturation conditions (pressure, liquid density and vapor density).

 Author: John Eslick
 File helmholtz_solve.cpp
------------------------------------------------------------------------------*/

#include<stdio.h>
#include<cmath>
#include<iostream>

#include"helmholtz_memo.h"
#include"helmholtz_config.h"
#include"helmholtz_deriv_parts.h"
#include"helmholtz_solve.h"

/*------------------------------------------------------------------------------
In this section you will find function for calculating density from pressure
and temperature. This lets you calculate properties in the more standard way as
a function of temperature and pressure instead of density and pressure.

Unfortunatly this is difficult and a bit messy.
*-----------------------------------------------------------------------------*/
s_real delta_p_tau_rf(s_real pr, s_real tau, s_real a, s_real b, bool bisect){
  /*----------------------------------------------------------------------------
  Bracketing methods, false position and bisection, for finding a better initial
  guess for density when solving for density from temperature and pressure. This
  is only used in particularly difficult areas.  At this point it probably
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
  // of a bracketing method but it was conveneint).
  if( fabs(T_c/tau - T_c) < 1e-7 && fabs(pr - P_c) < 1e-4) return 1;
  // solve f(delta, tau) = 0; f(delta, tau) = p(delta, tau) - pr
  fa = p(a, tau) - pr; // initial f(a, tau)
  fb = p(b, tau) - pr; // initial f(b, tau)
  while(it < MAX_IT_BRACKET && (b - a) > TOL_BRACKET){
    if(bisect) c = (a + b)/2.0; // bisection
    else c = b - fb*(b - a)/(fb - fa); //false position
    fc = p(c, tau) - pr; // calcualte f(c)
    if(fc*fa >= 0){a = c; fa = fc;}
    else{b = c; fb = fc;}
    ++it;
  }
  return (a+b)/2.0;
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
   delta: reduced density (accuracy depends on tolerance and function shape)
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
    if(hes != NULL){ // calculate Hessian if needed.
      hes[0] = -hesp[0]*grad[0]/gradp[0]/gradp[0];
      hes[1] = -(hesp[1] + hesp[0]*grad[1])/gradp[0]/gradp[0];
      hes[2] = -(grad[0]*(hesp[2] + grad[1]*hesp[1]) + gradp[1]*hes[1]);}}
  return delta;
}

s_real delta_liq(s_real p, s_real tau, s_real *grad, s_real *hes, int *nit){
  /*----------------------------------------------------------------------------
  Get a good liquid or super critical initial density guess then call
  delta_p_tau() to calculate density. In difficult cases an interval with the
  desired soultion is used with a bracketing method to produce a better initial
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
   delta: reduced density (accuracy depends on tolerance and function shape)
  ----------------------------------------------------------------------------*/
  s_real val = memoize::get_bin(memoize::DL_FUNC, p, tau, grad, hes);
  if(!std::isnan(val)) return val;
  s_real delta;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  delta = delta_p_tau(p, tau, LIQUID_DELTA_GUESS, TOL_DELTA_LIQ, nit, grad, hes); //solve
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){
    // This is just to avoid evaluation errors.  Want to be able to calculate
    // vapor properties even when vapor doesn't exist.  In the IDAES Framework
    // these properties may be calculated and multipled by a zero liquid fraction,
    // so it doesn't matter that they are wrong.
    delta = 3.1;
    zero_derivs2(grad, hes);
  }
  memoize::add_bin(memoize::DL_FUNC, p, tau, delta, grad, hes); //store
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   //   function
  return delta;
}

s_real delta_vap(s_real p, s_real tau, s_real *grad, s_real *hes, int *nit){
  /*----------------------------------------------------------------------------
  Get a good vapor or super critical initial density guess then call
  delta_p_tau() to calculate density. In the supercritical region this just
  calls the liquid function. In the rest of the vapor region the initial guess
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
   delta: reduced density (accuracy depends on tolerance and function shape)
  ----------------------------------------------------------------------------*/
  s_real val = memoize::get_bin(memoize::DV_FUNC, p, tau, grad, hes);
  if(!std::isnan(val)) return val; // return stored result if available
  s_real delta;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  delta = delta_p_tau(p, tau, VAPOR_DELTA_GUESS, TOL_DELTA_VAP, nit, grad, hes);
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){
    // This is just to avoid evaluation errors.  Want to be able to calculate
    // vapor properties even when vapor doesn't exist.  In the IDAES Framework
    // these properties may be calculated and multipled by a zero vapor fraction,
    // so it doesn't matter that they are wrong.
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
Saturation pressure and density as a function of temperature.
*-----------------------------------------------------------------------------*/
s_real sat_delta_liq(s_real tau){ //calculate saturated liquid density from tau
  s_real delta_l, delta_v;
  sat(tau, &delta_l, &delta_v);
  return delta_l;
}

s_real sat_delta_vap(s_real tau){ //calculate saturated vapor density from tau
  s_real delta_l, delta_v;
  sat(tau, &delta_l, &delta_v);
  return delta_v;
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
    //Get saturated phase densities at tau by Akasaka (2008) method
    s_real delta_l, delta_v, fg, gradl[1], hesl[1], gradv[1], hesv[1];
    int n=0, max_it=MAX_IT_SAT;

    if(tau - 1 < 1e-12){
      delta_l = 1.0;
      delta_v = 1.0;
      max_it=0;
    }
    else{
      // okay so you've decided to solve this thing
      delta_l = DELTA_LIQ_SAT_GUESS;
      delta_v = DELTA_VAP_SAT_GUESS;
    }
    // Since the equilibrium conditions are gl = gv and pl = pv, I am using the
    // the relative difference in g as a convergence criteria, that is easy to
    // understand.  fg < tol for convergence, fg is calculated upfront in the
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


/*------------------------------------------------------------------------------
We often use enthaply and pressure as state variables so there are functions for
calculatig Tsat from P and T from H and P below
-------------------------------------------------------------------------------*/
s_real sat_tau_with_derivs(s_real pr, s_real *grad, s_real *hes, int *nit){
  //Before getting into the real calculation, check if outside the allowed range
  //of 240K to T_c, the low end of this range doen't mean anything.  Its too cold
  //to expect liquid, but the calculations hold up there so for numerical reasons
  //I'll allow it.  Above the critical temperature there is only a single phase
  if(pr > P_c){ // above critical T
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return 1.0;
  }
  if(pr < P_t){ // above critical T
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return T_c/T_t;
  }
  //Now check if there is a stored value
  s_real val = memoize::get_un(memoize::TAU_SAT_FUNC, pr, grad, hes);
  if(!std::isnan(val)) return val;
  // grad and/or hes not provided so allocate for memoization
  bool free_grad = 0, free_hes = 0;
  if(grad==NULL){grad = new s_real[1]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[1]; free_hes = 1;}
  s_real tau = 1.5, fun, gradp[1], hesp[1], tol=TOL_SAT_TAU;
  if(P_c - pr < 1e-3){
    pr = P_c - 1e-3;
  }
  int it = 0; // iteration count
  fun = sat_p_with_derivs(tau, gradp, hesp, 0) - pr;
  while(fabs(fun) > tol && it < MAX_IT_SAT_TAU){
    tau = tau - fun*gradp[0]/(gradp[0]*gradp[0] - 0.5*fun*hesp[0]);
    fun = sat_p_with_derivs(tau, gradp, hesp, 0) - pr;
    ++it;
  }
  grad[0] = 1.0/gradp[0];
  hes[0] = -grad[0]*grad[0]*grad[0]*hesp[0];
  memoize::add_un(memoize::TAU_SAT_FUNC, pr, tau, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   // function
  if(nit != NULL) *nit = it;
  return tau;
}

s_real tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes){

    s_real val = memoize::get_bin(memoize::TAU_FUNC, ht, pr, grad, hes);
    if(!std::isnan(val)) return val;

    //Even if you aren't asking for derivatives, calculate them for memo
    bool free_grad = 0, free_hes = 0;
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

    s_real tau_sat, hv=1.0, hl=1.0, fun, tau, gradh[2], hesh[3], tol = 1e-11;
    int it = 0, max_it = 20;

    tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
    if(pr <= P_c && pr >= P_t){
      hv = hvpt_with_derivs(pr, tau_sat, NULL, NULL);
      hl = hlpt_with_derivs(pr, tau_sat, NULL, NULL);
    }

    if( (hl > ht && pr > P_t) || pr > P_c){ // to see if it's liquid check enthalpy and make sure above tripple point
      if (pr >= P_c){
        tau = 1.5;
      }
      else{
        tau = T_c/T_t;
      }
      if(hlpt_with_derivs(pr, tau, gradh, hesh) - ht < 0 && pr < P_c){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Tsat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = tau_sat;
        b = tau;
        fa = hlpt_with_derivs(pr, a, gradh, hesh) - ht;
        fb = hlpt_with_derivs(pr, b, gradh, hesh) - ht;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = hlpt_with_derivs(pr, c, gradh, hesh) - ht;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        tau = (a+b)/2.0;
      }
      fun = hlpt_with_derivs(pr, tau, gradh, hesh) - ht;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = hlpt_with_derivs(pr, tau, gradh, hesh) - ht;
        ++it;
      }
    }
    else if (hv < ht  || pr < P_t){
      tau = T_c/(T_c/tau_sat + 100);
      if(hvpt_with_derivs(pr, tau, gradh, hesh) - ht > 0 && (pr > P_t)){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Tsat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = tau;
        b = tau_sat;
        fa = hvpt_with_derivs(pr, a, gradh, hesh) - ht;
        fb = hvpt_with_derivs(pr, b, gradh, hesh) - ht;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = hvpt_with_derivs(pr, c, gradh, hesh) - ht;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        tau = (a+b)/2.0;
      }
      fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
        ++it;
      }
    }
    else{
      zero_derivs2(grad, hes);
      return tau_sat;
    }
    if(tau < 0.0 || tau > TAU_HIGH){
        std::cerr << "WARNING: External Helmholtz EOS low temperature clip, h= " << ht << " P= " << pr << " T= " << T_c/tau << " Tsat= " << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    else if(tau < TAU_LOW){
        std::cerr << "WARNING: External Helmholtz EOS high temperature clip, h= " << ht << " P= " << pr << " T= " << T_c/tau << " Tsat= " << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    if(grad != NULL){
        grad[0] = 1.0/gradh[1];
        grad[1] = -grad[0]*gradh[0];
        if(hes != NULL){
          hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
          hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
          hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
        }
    }
    memoize::add_bin(memoize::TAU_FUNC, ht, pr, tau, grad, hes);
    // If we alocated grad and hes here, free them
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
    return tau;
}


s_real tau_from_sp_with_derivs(s_real st, s_real pr, s_real *grad, s_real *hes){

    s_real val = memoize::get_bin(memoize::TAU_ENTR_FUNC, st, pr, grad, hes);
    if(!std::isnan(val)) return val;

    //Even if you aren't asking for derivatives, calculate them for memo
    bool free_grad = 0, free_hes = 0;
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
    s_real tau_sat, sv=1.0, sl=1.0, fun, tau, gradh[2], hesh[3], tol = 1e-11;
    int it = 0, max_it = 20;

    tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
    if(pr <= P_c && pr >= P_t){
      sv = svpt_with_derivs(pr, tau_sat, NULL, NULL);
      sl = slpt_with_derivs(pr, tau_sat, NULL, NULL);
    }
    if( (sl > st && pr > P_t) || pr > P_c){ // to see if it's liquid check enthalpy and make sure above tripple point
      if (pr >= P_c){
        tau = 1.5;
      }
      else{
        tau = T_c/T_t;
      }
      if(slpt_with_derivs(pr, tau, gradh, hesh) - st < 0 && pr < P_c){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Tsat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = tau_sat;
        b = tau;
        fa = slpt_with_derivs(pr, a, gradh, hesh) - st;
        fb = slpt_with_derivs(pr, b, gradh, hesh) - st;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = slpt_with_derivs(pr, c, gradh, hesh) - st;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        tau = (a+b)/2.0;
      }
      fun = slpt_with_derivs(pr, tau, gradh, hesh) - st;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = slpt_with_derivs(pr, tau, gradh, hesh) - st;
        ++it;
      }
    }
    else if (sv < st  || pr < P_t){
      tau = T_c/(T_c/tau_sat + 100);
      if(svpt_with_derivs(pr, tau, gradh, hesh) - st > 0 && (pr > P_t)){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Tsat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = tau;
        b = tau_sat;
        fa = svpt_with_derivs(pr, a, gradh, hesh) - st;
        fb = svpt_with_derivs(pr, b, gradh, hesh) - st;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = svpt_with_derivs(pr, c, gradh, hesh) - st;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        tau = (a+b)/2.0;
      }
      fun = svpt_with_derivs(pr, tau, gradh, hesh) - st;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = svpt_with_derivs(pr, tau, gradh, hesh) - st;
        ++it;
      }
    }
    else{
      zero_derivs2(grad, hes);
      return tau_sat;
    }
    if(tau < 0.0 || tau > TAU_HIGH){
        std::cerr << "WARNING: External Helmholtz EOS low temperature clip, s= " << st << " P= " << pr << " T= " << T_c/tau << " Tsat= " << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    else if(tau < TAU_LOW){
        std::cerr << "WARNING: External Helmholtz EOS high temperature clip, s= " << st << " P= " << pr << " T= " << T_c/tau << " Tsat= " << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    if(grad != NULL){
        grad[0] = 1.0/gradh[1];
        grad[1] = -grad[0]*gradh[0];
        if(hes != NULL){
          hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
          hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
          hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
        }
    }
    memoize::add_bin(memoize::TAU_ENTR_FUNC, st, pr, tau, grad, hes);
    // If we alocated grad and hes here, free them
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
    return tau;
}

s_real tau_from_up_with_derivs(s_real ut, s_real pr, s_real *grad, s_real *hes){

    s_real val = memoize::get_bin(memoize::TAU_INTEN_FUNC, ut, pr, grad, hes);
    if(!std::isnan(val)) return val;

    //Even if you aren't asking for derivatives, calculate them for memo
    bool free_grad = 0, free_hes = 0;
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
    s_real tau_sat, uv=1.0, ul=1.0, fun, tau, gradh[2], hesh[3], tol = 1e-11;
    int it = 0, max_it = 20;

    tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
    if(pr <= P_c && pr >= P_t){
      uv = uvpt_with_derivs(pr, tau_sat, NULL, NULL);
      ul = ulpt_with_derivs(pr, tau_sat, NULL, NULL);
    }
    if( (ul > ut && pr > P_t) || pr > P_c){ // to see if it's liquid check enthalpy and make sure above tripple point
      if (pr >= P_c){
        tau = 1.5;
      }
      else{
        tau = T_c/T_t;
      }
      if(ulpt_with_derivs(pr, tau, gradh, hesh) - ut < 0 && pr < P_c){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Tsat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = tau_sat;
        b = tau;
        fa = ulpt_with_derivs(pr, a, gradh, hesh) - ut;
        fb = ulpt_with_derivs(pr, b, gradh, hesh) - ut;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = ulpt_with_derivs(pr, c, gradh, hesh) - ut;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        tau = (a+b)/2.0;
      }
      fun = ulpt_with_derivs(pr, tau, gradh, hesh) - ut;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = ulpt_with_derivs(pr, tau, gradh, hesh) - ut;
        ++it;
      }
    }
    else if (uv < ut  || pr < P_t){
      tau = T_c/(T_c/tau_sat + 100);
      if(uvpt_with_derivs(pr, tau, gradh, hesh) - ut > 0 && (pr > P_t)){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Tsat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = tau;
        b = tau_sat;
        fa = uvpt_with_derivs(pr, a, gradh, hesh) - ut;
        fb = uvpt_with_derivs(pr, b, gradh, hesh) - ut;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = uvpt_with_derivs(pr, c, gradh, hesh) - ut;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        tau = (a+b)/2.0;
      }
      fun = uvpt_with_derivs(pr, tau, gradh, hesh) - ut;
      while(fabs(fun) > tol && it < max_it){
        tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
        fun = uvpt_with_derivs(pr, tau, gradh, hesh) - ut;
        ++it;
      }
    }
    else{
      zero_derivs2(grad, hes);
      return tau_sat;
    }
    if(tau < 0.0 || tau > TAU_HIGH){
        std::cerr << "WARNING: External Helmholtz EOS low temperature clip, u= " << ut << " P= " << pr << " T= " << T_c/tau << " Tsat= " << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    else if(tau < TAU_LOW){
        std::cerr << "WARNING: External Helmholtz EOS high temperature clip, u= " << ut << " P= " << pr << " T= " << T_c/tau << " Tsat= " << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    if(grad != NULL){
        grad[0] = 1.0/gradh[1];
        grad[1] = -grad[0]*gradh[0];
        if(hes != NULL){
          hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
          hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
          hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
        }
    }
    memoize::add_bin(memoize::TAU_INTEN_FUNC, ut, pr, tau, grad, hes);
    // If we alocated grad and hes here, free them
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
    return tau;
}

s_real p_from_htau_with_derivs(s_real ht, s_real tau, s_real *grad, s_real *hes){

    s_real val = memoize::get_bin(memoize::P_ENTH_FUNC, ht, tau, grad, hes);
    if(!std::isnan(val)) return val;

    //Even if you aren't asking for derivatives, calculate them for memo
    bool free_grad = 0, free_hes = 0;
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
    s_real p_sat, hv=1.0, hl=1.0, fun, pr, gradh[2], hesh[3], tol = 1e-11, T=T_c/tau;
    int it = 0, max_it = 20;

    std::cerr << "Hi in p" << std::endl;

    p_sat = sat_p_with_derivs(tau, NULL, NULL);
    std::cerr << "p_sat " << p_sat << std::endl;
    if(T <= T_c && T >= T_t){
      hv = hvpt_with_derivs(p_sat, tau, NULL, NULL);
      hl = hlpt_with_derivs(p_sat, tau, NULL, NULL);
    }
    std::cerr << "hl, ht" << hl << " " << ht << std::endl;
    if( (hl > ht && T > T_t) || T >= T_c ){ // to see if it's liquid check enthalpy and make sure above tripple point
      if (T >= T_c){
        pr = P_c*1.2;
      }
      else{
        pr = (p_sat + P_c)/2.0;
      }
      std::cerr << "liquid P = " << pr << std::endl;

      if(hlpt_with_derivs(pr, tau, gradh, hesh) - ht < 0 && T < T_c){
        // Unfortunatly if the initial guess isn't good you can get on the wrong
        // side of Psat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = p_sat;
        b = pr;
        fa = hlpt_with_derivs(a, tau, gradh, hesh) - ht;
        fb = hlpt_with_derivs(b, tau, gradh, hesh) - ht;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = hlpt_with_derivs(c, tau, gradh, hesh) - ht;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        pr = (a+b)/2.0;
      }
      std::cerr << "after bracket liquid P = " << pr << std::endl;
      fun = hlpt_with_derivs(pr, tau, gradh, hesh) - ht;
      while(fabs(fun) > tol && it < max_it){
        pr = pr - fun*gradh[0]/(gradh[0]*gradh[0] - 0.5*fun*hesh[0]);
        fun = hlpt_with_derivs(pr, tau, gradh, hesh) - ht;
        ++it;
      }
    }
    else if (hv < ht  || T < T_t){
      pr = P_t;
      std::cerr << "vap P = " << std::endl;

      if(hvpt_with_derivs(pr, tau, gradh, hesh) - ht > 0 && (T > T_t)){
        // Unfotunatly if the initial guess isn't good you can get on the wrong
        // side of Psat, then you have trouble. This false position method up
        // front keeps the temperature on the right side while refining the
        // guess.  With the better guess the newton method shouldn't get out of
        // control
        s_real a, b, c, fa, fb, fc;
        a = pr;
        b = p_sat;
        fa = hvpt_with_derivs(a, tau, gradh, hesh) - ht;
        fb = hvpt_with_derivs(a, tau, gradh, hesh) - ht;
        for(it=0;it<5;++it){
          c = b - fb*(b - a)/(fb - fa);
          fc = hvpt_with_derivs(c, tau, gradh, hesh) - ht;
          if(fc*fa >= 0){a = c; fa = fc;}
          else{b = c; fb = fc;}
          if(b - a < 1e-4) {break;}
        }
        pr = (a+b)/2.0;
      }
      fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
      while(fabs(fun) > tol && it < max_it){
        pr = pr - fun*gradh[0]/(gradh[0]*gradh[0] - 0.5*fun*hesh[0]);
        fun = hvpt_with_derivs(pr, tau, gradh, hesh) - ht;
        ++it;
      }
    }
    else{
      zero_derivs2(grad, hes);
      return p_sat;
    }
    if(pr > P_HIGH){
        std::cerr << "WARNING: External Helmholtz EOS high pressure clip, h= " << ht << " P= " << pr << " T= " << T_c/tau << " Psat= " << p_sat << std::endl;
        return 0.0/0.0;
    }
    else if(pr <= P_LOW){
        std::cerr << "WARNING: External Helmholtz EOS low pressure clip, h= " << ht << " P= " << pr << " T= " << T_c/tau << " Psat= " << p_sat << std::endl;
        return 0.0/0.0;
    }
    if(grad != NULL){
        grad[0] = 1.0/gradh[0];
        grad[1] = -grad[0]*gradh[0];
        if(hes != NULL){
          hes[0] = -grad[0]*grad[0]*grad[0]*hesh[0];
          hes[1] = -grad[0]*grad[0]*(hesh[0] + hesh[0]*grad[1]);
          hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
        }
    }
    memoize::add_bin(memoize::P_ENTH_FUNC, ht, tau, pr, grad, hes);
    // If we alocated grad and hes here, free them
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
    return tau;
}

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
  Function wrapper functor
------------------------------------------------------------------------------*/

FuncWrapper::FuncWrapper(char apos, s_real a, s_real c){
  _apos = apos;
  _a = a;
  _c = c;
  grad_pos = apos;  // this is used for binary->unary functions
  hes_pos = 2*apos; // this is used for binary->unary functions
}

void FuncWrapper::set_f1(f_ptr1 f1){
  this->_f1 = f1;
}

void FuncWrapper::set_f2(f_ptr2 f2){
  this->_f2 = f2;
}

void FuncWrapper::set_f2n(f_ptr2_nod f2){
  this->_f2n = f2;
}

s_real FuncWrapper:: operator () (s_real x, s_real *grad, s_real *hes){
  // There is a lack of error checking here, but since this is ment for internal
  // use, I'll be careful to use it right.
  if(this->_f1){
    return (*_f1)(x, grad, hes) - this->_c;
  }
  else if(this->_apos==0){
    return (*_f2)(x, this->_a, grad, hes) - this->_c;
  }
  return (*_f2)(this->_a, x, grad, hes) - this->_c;
}

s_real FuncWrapper:: operator () (s_real x){
  // There is a lack of error checking here, but since this is ment for internal
  // use, I'll be careful to use it right.
  if(this->_f1){
    return (*_f1)(x, NULL, NULL) - this->_c;
  }
  else if(this->_f2n){
    return (*_f2n)(x, this->_a) - this->_c;
  }
  else if(this->_apos==0){
    return (*_f2)(x, this->_a, NULL, NULL) - this->_c;
  }
  return (*_f2)(this->_a, x, NULL, NULL) - this->_c;
}

s_real FuncWrapper:: operator () (s_real x0, s_real x1, s_real *grad, s_real *hes){
  return (*_f2)(x0, x1, grad, hes) - this->_c;
}

s_real FuncWrapper:: operator () (s_real x0, s_real x1){
  if(this->_f2n){
    return (*_f2n)(x0, x1) - this->_c;
  }
  return (*_f2)(x0, x1, NULL, NULL) - this->_c;
}

/*------------------------------------------------------------------------------
 False-position type bracketing solver.  This is mostly used to confine the
 search for a good starting point to a specific area, where a newton method
 given a starting point too far way may not converge or converge to the wrong
 solution. It may be a bit slow, but this can also be used for solve
 a problem competely.
------------------------------------------------------------------------------*/
int bracket(FuncWrapper *f, s_real xa, s_real xb, s_real *sol,
            int max_it, s_real ftol, s_real xtol){
  s_real fa, fb, fc, a=xa, b=xb, c;
  int it;
  bool prev_a=0, prev_b=0;

  fa = (*f)(a);
  fb = (*f)(b);

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
    fc = (*f)(c);
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
int halley(FuncWrapper *f, s_real x0, s_real *sol, s_real *grad, s_real *hes,
           int max_it, s_real ftol){
  int it=0;
  s_real fun = (*f)(x0, grad, hes);
  s_real x = x0;

  while(fabs(fun) > ftol && it < max_it){
   x = x - fun*grad[f->grad_pos]/
           (grad[f->grad_pos]*grad[f->grad_pos] - 0.5*fun*hes[f->hes_pos]);
   fun = (*f)(x, grad, hes);
   //std::cerr << it << " f = " << fun << " P = " << pr << std::endl;
   ++it;
  }
  *sol = x;
  if(it == max_it){
    return -3;
  }
  return it;
}

/*------------------------------------------------------------------------------
  1D Newton's method with backtracking line search
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  2D Newton's method with backtracking line search
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
In this section you will find function for calculating density from pressure
and temperature. This lets you calculate properties in the more standard way as
a function of temperature and pressure instead of density and pressure.

Unfortunatly this is difficult and a bit messy.
*-----------------------------------------------------------------------------*/
s_real delta_p_tau_rf(s_real pr, s_real tau, s_real a, s_real b){
  /*----------------------------------------------------------------------------
  Bracketing methods, false position and bisection, for finding a better initial
  guess for density when solving for density from temperature and pressure. This
  is only used in particularly difficult areas.  At this point it probably
  overused, but there are places where it is probably necessary.  This is used
  in the impimentation specific xxx_guess.cpp files.

  Args:
    pr: pressure (kPa)
    tau: inverse of reduced pressure Tc/T
    a: first density bound (kg/m3)  | really density, why didn't use delta?
    b: second density bound (kg/m3) | it was easier to think in terms of density
  Returns:
    delta: reduced density at pr and tau (more or less approximate)
  ----------------------------------------------------------------------------*/
  s_real c;
  // If right by critical point guess critical density. (okay this isn't part
  // of a bracketing method but it was conveneint).
  if( fabs(T_c/tau - T_c) < 1e-7 && fabs(pr - P_c) < 1e-4) return 1;

  //Okay, now do bracket
  FuncWrapper f(0, tau, pr);
  f.set_f2n(&p);
  bracket(&f, a, b, &c, MAX_IT_BRACKET, TOL_BRACKET, 1e-6);
  return c;
}

s_real delta_liq(s_real pr, s_real tau, s_real *grad, s_real *hes, int *nit){
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
  s_real val = memoize::get_bin(memoize::DL_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  s_real delta, gradp[2], hesp[3];
  int it;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  //Okay, now do bracket
  FuncWrapper f(0, tau, pr);
  f.set_f2(&p_with_derivs);
  it = halley(&f, LIQUID_DELTA_GUESS, &delta, gradp, hesp, MAX_IT_DELTA, TOL_DELTA_LIQ);
  if(nit) {*nit = it;}
  // Error check, want a number even if phase doesn't exist
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){ //avoid eval errors
    delta = 5.0;
    zero_derivs2(grad, hes);
    memoize::add_bin(memoize::DL_FUNC, pr, tau, delta, grad, hes); //store
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
  }

  grad[0] = 1.0/gradp[0];
  grad[1] = -gradp[1]*grad[0]; //triple product
  hes[0] = -hesp[0]*grad[0]/gradp[0]/gradp[0];
  hes[1] = -(hesp[1] + hesp[0]*grad[1])/gradp[0]/gradp[0];
  hes[2] = -(grad[0]*(hesp[2] + grad[1]*hesp[1]) + gradp[1]*hes[1]);

  memoize::add_bin(memoize::DL_FUNC, pr, tau, delta, grad, hes); //store
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   // function
  return delta;
}

s_real delta_vap(s_real pr, s_real tau, s_real *grad, s_real *hes, int *nit){
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
  s_real val = memoize::get_bin(memoize::DV_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val; // return stored result if available
  s_real delta, gradp[2], hesp[3];
  int it=0;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}
  //Okay, now do bracket
  FuncWrapper f(0, tau, pr);
  f.set_f2(&p_with_derivs);
  it = halley(&f, VAPOR_DELTA_GUESS, &delta, gradp, hesp, MAX_IT_DELTA, TOL_DELTA_VAP);
  if(nit != NULL) *nit = it;
  if(nit) {*nit = it;}
  // Error check, want a number even if phase doesn't exist
  if(std::isnan(delta) || delta < 1e-12 || delta > 5.0){ //avoid eval errors
    delta = 0.001;
    zero_derivs2(grad, hes);
    memoize::add_bin(memoize::DV_FUNC, pr, tau, delta, grad, hes); // store result
    if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
    if(free_hes) delete[] hes; // function
  }

  grad[0] = 1.0/gradp[0];
  grad[1] = -gradp[1]*grad[0]; //triple product
  hes[0] = -hesp[0]*grad[0]/gradp[0]/gradp[0];
  hes[1] = -(hesp[1] + hesp[0]*grad[1])/gradp[0]/gradp[0];
  hes[2] = -(grad[0]*(hesp[2] + grad[1]*hesp[1]) + gradp[1]*hes[1]);

  memoize::add_bin(memoize::DV_FUNC, pr, tau, delta, grad, hes); // store result
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; // function
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
  if(pr < P_t){ // below the tripple point pressure
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


s_real mem_tau_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes){
    // TESTING ONLY, function for testing memoization
    tau_with_derivs(ht, pr, grad, hes);  //after this it should be memoized
    //return the memoized value
    return memoize::get_bin(memoize::TAU_FUNC, ht, pr, grad, hes);
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
    s_real (*fun_ptr)(s_real, s_real, s_real*, s_real*);

    tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
    if(pr >= P_t){
      hv = hvpt_with_derivs(pr, tau_sat, NULL, NULL);
      hl = hlpt_with_derivs(pr, tau_sat, NULL, NULL);
    }

    if (ht >= hl && ht <= hv){
      zero_derivs2(grad, hes);
      return tau_sat;
    }

    s_real a, b, c, fa, fb, fc;
    bool prev_a=0, prev_b=0;

    if (pr > P_c){
      a = 0.5;
      b = T_c/T_t;
      tau = a;
      fun_ptr = &hlpt_with_derivs;
      //std::cerr << "Liq P >Pc Tsat = " << T_c/tau_sat << std::endl;
    }
    else if (hl > ht && pr > P_t){ // liquid
      a = T_c/T_t;
      b = tau_sat;
      tau = a;
      fun_ptr = &hlpt_with_derivs;
      //std::cerr << "Liq Tsat = " << T_c/tau_sat << std::endl;
    }
    else{ // vapor
      tau = T_c/(T_c/tau_sat + 100);
      a = tau_sat;
      b = tau;
      fun_ptr = &hvpt_with_derivs;
      //std::cerr << "Vap Tsat = " << T_c/tau_sat << std::endl;
    }
    fa = (*fun_ptr)(pr, a, gradh, hesh) - ht;
    fb = (*fun_ptr)(pr, b, gradh, hesh) - ht;
    if (fa*fb < 0){
      for(it=0;it<15;++it){
        c = b - fb*(b - a)/(fb - fa);
        fc = (*fun_ptr)(pr, c, gradh, hesh) - ht;
        if(fc*fa >= 0){
          a = c;
          fa = fc;
          if (prev_a){
            fb *= 0.5;
          }
          prev_a = 1;
          prev_b = 0;
        }
        else{
          b = c;
          fb = fc;
          if (prev_b){
            fa *= 0.5;
          }
          prev_a = 0;
          prev_b = 1;
        }
        if (fabs(fa) < tol) {tau=a; break;}
        if (fabs(fb) < tol) {tau=b; break;}
        tau = (a + b)/2.0;
        if(fabs(b - a) < 1e-5) {break;}
        //std::cerr << it << " fa = " << fa << " fb = " << fb;
        //std::cerr << " T = " << T_c/c << std::endl;
      }
    }
    it = 0;
    //std::cerr << "Tinit = " << T_c/tau << std::endl;
    fun = (*fun_ptr)(pr, tau, gradh, hesh) - ht;
    while(fabs(fun) > tol && it < max_it){
      tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
      fun = (*fun_ptr)(pr, tau, gradh, hesh) - ht;
      //std::cerr << it << " f = " << fun << " T = " << T_c/tau << std::endl;
      ++it;
    }

    if(tau < 0.0 || tau > TAU_HIGH){
        std::cerr << "WARNING: External Helmholtz EOS low temperature clip, h= ";
        std::cerr << ht << " P= " << pr << " T= " << T_c/tau << " Tsat= ";
        std::cerr << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    else if(tau < TAU_LOW){
        std::cerr << "WARNING: External Helmholtz EOS high temperature clip, h= ";
        std::cerr << ht << " P= " << pr << " T= " << T_c/tau << " Tsat= ";
        std::cerr << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }

    grad[0] = 1.0/gradh[1];
    grad[1] = -grad[0]*gradh[0];
    hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
    hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
    hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
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
  s_real (*fun_ptr)(s_real, s_real, s_real*, s_real*);

  tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
  if(pr >= P_t){
    sv = svpt_with_derivs(pr, tau_sat, NULL, NULL);
    sl = slpt_with_derivs(pr, tau_sat, NULL, NULL);
  }

  if (st >= sl && st <= sv){
    zero_derivs2(grad, hes);
    return tau_sat;
  }

  s_real a, b, c, fa, fb, fc;
  bool prev_a=0, prev_b=0;

  if (pr > P_c){
    a = 0.5;
    b = T_c/T_t;
    tau = a;
    fun_ptr = &slpt_with_derivs;
    //std::cerr << "Liq P >Pc Tsat = " << T_c/tau_sat << std::endl;
  }
  else if (sl > st && pr > P_t){ // liquid
    a = T_c/T_t;
    b = tau_sat;
    tau = a;
    fun_ptr = &slpt_with_derivs;
    //std::cerr << "Liq Tsat = " << T_c/tau_sat << std::endl;
  }
  else{ // vapor
    tau = T_c/(T_c/tau_sat + 100);
    a = tau_sat;
    b = tau;
    fun_ptr = &svpt_with_derivs;
    //std::cerr << "Vap Tsat = " << T_c/tau_sat << std::endl;
  }
  fa = (*fun_ptr)(pr, a, gradh, hesh) - st;
  fb = (*fun_ptr)(pr, b, gradh, hesh) - st;
  if (fa*fb < 0){
    for(it=0;it<15;++it){
      c = b - fb*(b - a)/(fb - fa);
      fc = (*fun_ptr)(pr, c, gradh, hesh) - st;
      if(fc*fa >= 0){
        a = c;
        fa = fc;
        if (prev_a){
          fb *= 0.5;
        }
        prev_a = 1;
        prev_b = 0;
      }
      else{
        b = c;
        fb = fc;
        if (prev_b){
          fa *= 0.5;
        }
        prev_a = 0;
        prev_b = 1;
      }
      if (fabs(fa) < tol) {tau=a; break;}
      if (fabs(fb) < tol) {tau=b; break;}
      tau = (a + b)/2.0;
      if(fabs(b - a) < 1e-5) {break;}
      //std::cerr << it << " fa = " << fa << " fb = " << fb;
      //std::cerr << " T = " << T_c/c << std::endl;
    }
  }
  it = 0;
  //std::cerr << "Tinit = " << T_c/tau << std::endl;
  fun = (*fun_ptr)(pr, tau, gradh, hesh) - st;
  while(fabs(fun) > tol && it < max_it){
    tau = tau - fun*gradh[1]/(gradh[1]*gradh[1] - 0.5*fun*hesh[2]);
    fun = (*fun_ptr)(pr, tau, gradh, hesh) - st;
    //std::cerr << it << " f = " << fun << " T = " << T_c/tau << std::endl;
    ++it;
  }

  if(tau < 0.0 || tau > TAU_HIGH){
      std::cerr << "WARNING: External Helmholtz EOS low temperature clip, s= ";
      std::cerr << st << " P= " << pr << " T= " << T_c/tau << " Tsat= ";
      std::cerr << T_c/tau_sat << std::endl;
      return 0.0/0.0;
  }
  else if(tau < TAU_LOW){
      std::cerr << "WARNING: External Helmholtz EOS high temperature clip, s= ";
      std::cerr << st << " P= " << pr << " T= " << T_c/tau << " Tsat= ";
      std::cerr << T_c/tau_sat << std::endl;
      return 0.0/0.0;
  }

  grad[0] = 1.0/gradh[1];
  grad[1] = -grad[0]*gradh[0];
  hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
  hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
  hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
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

    s_real tau_sat, uv=1.0, ul=1.0, tau, gradh[2], hesh[3], a, b;
    f_ptr2 fun_ptr;

    tau_sat = sat_tau_with_derivs(pr, NULL, NULL);
    if(pr >= P_t){
      uv = uvpt_with_derivs(pr, tau_sat, NULL, NULL);
      ul = ulpt_with_derivs(pr, tau_sat, NULL, NULL);
    }

    if (ut >= ul && ut <= uv){
      zero_derivs2(grad, hes);
      memoize::add_bin(memoize::TAU_INTEN_FUNC, ut, pr, tau_sat, grad, hes);
      if(free_grad) delete[] grad;
      if(free_hes) delete[] hes;
      return tau_sat;
    }

    if (pr > P_c){
      a = 0.5;
      b = T_c/T_t;
      tau = a;
      fun_ptr = &ulpt_with_derivs;
      //std::cerr << "Liq P >Pc Tsat = " << T_c/tau_sat << std::endl;
    }
    else if (ul > ut && pr > P_t){ // liquid
      a = T_c/T_t;
      b = tau_sat;
      tau = a;
      fun_ptr = &ulpt_with_derivs;
      //std::cerr << "Liq Tsat = " << T_c/tau_sat << std::endl;
    }
    else{ // vapor
      tau = T_c/(T_c/tau_sat + 100);
      a = tau_sat;
      b = tau;
      fun_ptr = &uvpt_with_derivs;
      //std::cerr << "Vap Tsat = " << T_c/tau_sat << std::endl;
    }

    FuncWrapper f(1, pr, ut);
    f.set_f2(fun_ptr);
    bracket(&f, a, b, &tau, 10, 1e-5, 1e-5);
    halley(&f, tau, &tau, gradh, hesh, 15, 1e-11);

    if(tau < 0.0 || tau > TAU_HIGH){
        std::cerr << "WARNING: External Helmholtz EOS low temperature clip, u= ";
        std::cerr << ut << " P= " << pr << " T= " << T_c/tau << " Tsat= ";
        std::cerr << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }
    else if(tau < TAU_LOW){
        std::cerr << "WARNING: External Helmholtz EOS high temperature clip, u= ";
        std::cerr << ut << " P= " << pr << " T= " << T_c/tau << " Tsat= ";
        std::cerr << T_c/tau_sat << std::endl;
        return 0.0/0.0;
    }

    grad[0] = 1.0/gradh[1];
    grad[1] = -grad[0]*gradh[0];
    hes[0] = -grad[0]*grad[0]*grad[0]*hesh[2];
    hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[2]*grad[1]);
    hes[2] = -hes[1]*gradh[0] - grad[0]*(hesh[0] + hesh[1]*grad[1]);
    memoize::add_bin(memoize::TAU_INTEN_FUNC, ut, pr, tau, grad, hes);
    // If we alocated grad and hes here, free them
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
    return tau;
}


s_real p_from_stau_with_derivs(s_real st, s_real tau, s_real *grad, s_real *hes){
    // STEP 1 -- Check for memoized answer
    s_real val = memoize::get_bin(memoize::P_ENTR_FUNC, st, tau, grad, hes);
    if(!std::isnan(val)) return val;

    // Define vars
    s_real p_sat, sv=1.0, sl=1.0, a, b, pr, gradh[2], hesh[3], T=T_c/tau;
    f_ptr2 fun_ptr=NULL;
    //Even if you aren't asking for derivatives, calculate them for memo
    bool free_grad = 0, free_hes = 0;
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

    // STEP 2 -- Determine phase and return easy 2-phase answer if applicable
    if(T >= T_t){
      sv = s_with_derivs(sat_delta_vap(tau), tau, NULL, NULL);
      sl = s_with_derivs(sat_delta_liq(tau), tau, NULL, NULL);
      //std::cerr << "sl, sv, st " << sl << ", " << sv << ", " << st << std::endl;
    }
    p_sat = sat_p_with_derivs(tau, NULL, NULL);
    if (st >= sl && st <= sv){
      zero_derivs2(grad, hes);
      memoize::add_bin(memoize::P_ENTR_FUNC, st, tau, p_sat, grad, hes);
      if(free_grad) delete[] grad;
      if(free_hes) delete[] hes;
      return p_sat;
    }

    // STEP 3 -- Determine initial guess or range to look for guess
    FuncWrapper f(0, tau, st);

    if (sl > st && T > T_t && T < T_c){ // liquid
      a = p_sat;
      b = P_max;
      fun_ptr = &slpt_with_derivs;
      //std::cerr << "Liq Psat = " << p_sat << std::endl;
    }
    else{ // vapor
      a = P_min;
      if(T > T_c){
        b = P_max;
      }
      else{
        b = p_sat;
      }
      fun_ptr = &svpt_with_derivs;
      //std::cerr << "Vap Psat = " << p_sat << std::endl;
    }

    //STEP 4 -- Solve 4a backet to get closer 4b solve
    f.set_f2(fun_ptr);
    bracket(&f, a, b, &pr, 15, 1e-5, 1e-5);
    halley(&f, pr, &pr, gradh, hesh, 15, 1e-11);

    //STEP 5 -- check answer for sanity
    if(pr > P_HIGH){
      std::cerr << "WARNING: External Helmholtz EOS high pressure clip, s= ";
      std::cerr << st << " P= " << pr << " T= " << T;
      std::cerr << " Psat= " << p_sat << std::endl;
      return 0.0/0.0;
    }
    else if(pr <= P_LOW){
      std::cerr << "WARNING: External Helmholtz EOS low pressure clip, s= ";
      std::cerr << st << " P= " << pr << " T= " << T;
      std::cerr << " Psat= " << p_sat << std::endl;
      return 0.0/0.0;
    }

    //STEP 6 -- calculate derivatives
    grad[0] = 1.0/gradh[0];
    grad[1] = -grad[0]*gradh[1];
    hes[0] = -grad[0]*grad[0]*grad[0]*hesh[0];
    hes[1] = -grad[0]*grad[0]*(hesh[1] + hesh[0]*grad[1]);
    hes[2] = -hes[1]*gradh[1] - grad[0]*(hesh[2] + hesh[1]*grad[1]);

    //STEP 7 -- store the answer
    memoize::add_bin(memoize::P_ENTR_FUNC, st, tau, pr, grad, hes);

    // If we allocated grad and hes here, free them
    if(free_grad) delete[] grad;
    if(free_hes) delete[] hes;
    return pr;
}

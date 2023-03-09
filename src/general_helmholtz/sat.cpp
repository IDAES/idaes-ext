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
Author: John Eslick
File sat.cpp

 This file contains the functions to solve the saturation curve.  The method
 is described in the following paper:

 Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
     State from Helmholtz Energy Equations of State." Journal of Thermal
     Science and Technology, 3(3), 442-451.
--------------------------------------------------------------------------------*/
//#define BOOST_MATH_INSTRUMENT
#include<boost/math/special_functions/relative_difference.hpp>
#include<boost/math/tools/roots.hpp>
#include<math.h>
#include"phi.h"
#include"props.h"
#include"config.h"
#include"sat.h"
#include<iostream>

prop_memo_table12 memo_table_sat_delta_l;
prop_memo_table12 memo_table_sat_delta_v;
prop_memo_table12 memo_table_sat_p;
prop_memo_table12 memo_table_sat_tau;

inline double J(double delta, f12_struct *phi){
  // Term from Akasaka method for saturation state
  return delta*(1 + delta*phi->f_1);
}

inline double K(double delta, f12_struct *phi){
  // Term from Akasaka method for saturation state
  return delta*phi->f_1 + phi->f + log(delta);
}

inline double J_delta(double delta, f12_struct *phi){
  // Derivative term from Akasaka method for saturation state
  return 1.0 + 2.0*delta*phi->f_1 + delta*delta*phi->f_11;
}

inline double K_delta(double delta, f12_struct *phi){
  // Derivative term from Akasaka method for saturation state
  return 2.0*phi->f_1 + delta*phi->f_11 + 1.0/delta;
}

int sat(uint comp, double tau, double *delta_l, double *delta_v){
  using boost::math::epsilon_difference;
  double Kdiff=1, Jdiff=1, Jv, Jl, Kv, Kl, det, dJv, dJl, dKv, dKl, detla_v_step, detla_l_step;
  unsigned int n=0;
  f12_struct phir_v, phir_l;
  double tol = 1e-12;
  double sat_gamma = 1.00;
  unsigned int max_iter = 10;

  if(tau <= tau_c(comp)){
    // So close to the critical point assume we are at the critical point, if
    // over the critical point this will also use the critical density
    *delta_l = delta_c(comp)*(1+1e-9);
    *delta_v = delta_c(comp)*(1-1e-9);
  }
  else{
    *delta_l = sat_delta_l_approx(comp, tau); // DELTA_LIQ_SAT_GUESS;
    *delta_v = sat_delta_v_approx(comp, tau); // DELTA_VAP_SAT_GUESS;
    if (tau < tau_c(comp)*1.001){
      sat_gamma = 0.25;
      max_iter = 20;
      tol = 1e-10;
    }
    //std::cout << "delta_l guess = " << *delta_l << " delta_v guess = " << *delta_v << std::endl;
    //std::cout << "rho_l guess = " << *delta_l*cdata[comp].rhoc << " rho_v guess = " << *delta_v*cdata[comp].rhoc << std::endl;

    while(n < max_iter){
      phir_v = phi_resi_for_sat(comp, *delta_v, tau);
      phir_l = phi_resi_for_sat(comp, *delta_l, tau);
      Jv = J(*delta_v, &phir_v);
      Jl = J(*delta_l, &phir_l);
      Kv = K(*delta_v, &phir_v);
      Kl = K(*delta_l, &phir_l);
      Jdiff = Jv - Jl;
      Kdiff = Kv - Kl;
      if (fabs(Jdiff/Jl) + fabs(Kdiff/Kl) <= tol){
        break;
      }
      dJv = J_delta(*delta_v, &phir_v);
      dJl = J_delta(*delta_l, &phir_l);
      dKv = K_delta(*delta_v, &phir_v);
      dKl = K_delta(*delta_l, &phir_l);
      det = dJv*dKl - dJl*dKv;
      ++n; // Count iterations

      detla_l_step = sat_gamma*(Kdiff*dJv - Jdiff*dKv)/det;
      detla_v_step = sat_gamma*(Kdiff*dJl - Jdiff*dKl)/det;
      if(detla_l_step + *delta_l <= delta_c(comp) || detla_v_step + *delta_v >= delta_c(comp)){
        // unfortunatly it seems these things can happen near the critical temperature
        // the derivatives blow up and it gets really hard to solve.  This should only
        // catch the cases right near the critical point, but it is possible that
        // it could hide other problems
        // std::cout << "sat curve solve stopped on out of bounds. tau = " << tau << std::endl;
        break; // should be every close
      }
      *delta_l += detla_l_step;
      *delta_v += detla_v_step;
      //std::cout << n << " delta_l = " << *delta_l << " delta_v = " << *delta_v << std::endl;
    }
  }
  return n;
}

int sat_delta_with_derivs(
  uint comp,
  double tau,
  f12_struct *delta_l_vec_ptr,
  f12_struct *delta_v_vec_ptr,
  f12_struct *p_vec_ptr
){
/*
Calculate grad and hes for and memoize

 To calculate the derivatives use the phase-equilibrium condition and the
 fact that the pressures are equal to implicitly calculate the derivatives. Tau
 is a constant.  There are a lot of intermediate varaibles here to hopfully
 break down the derivative calculations into easy(ish) to understand parts. The
 compiler optimization will hopfully take care of the inefficincy of the
 excessive number of variables.

*/
  double delta_l;
  double delta_v;

  // Calculate saturated delta_l and delta_v from tau.
  unsigned int n = sat(comp, tau, &delta_l, &delta_v);

  f23_struct phir_ls = phi_resi(comp, delta_l, tau);
  double phir_l_d = phir_ls.f_1;
  double phir_l_dd = phir_ls.f_11;
  double phir_l_ddd = phir_ls.f_111;
  double phir_l_t = phir_ls.f_2;
  double phir_l_dt = phir_ls.f_12;
  double phir_l_ddt = phir_ls.f_112;
  double phir_l_tt = phir_ls.f_22;
  double phir_l_dtt = phir_ls.f_122;
  f23_struct phir_vs = phi_resi(comp, delta_v, tau);
  double phir_v_d = phir_vs.f_1;
  double phir_v_dd = phir_vs.f_11;
  double phir_v_ddd = phir_vs.f_111;
  double phir_v_t = phir_vs.f_2;
  double phir_v_dt = phir_vs.f_12;
  double phir_v_ddt = phir_vs.f_112;
  double phir_v_tt = phir_vs.f_22;
  double phir_v_dtt = phir_vs.f_122;

  // A is the expression for p/(R*T*rho).  We'll use it to simplify writing
  // some expressions.
  double AL = delta_l * phir_l_d + 1;
  double AL_d = phir_l_d + delta_l*phir_l_dd;
  double AL_t = delta_l*phir_l_dt;
  double AL_dt = phir_l_dt + delta_l*phir_l_ddt;
  double AL_dd = 2*phir_l_dd + delta_l*phir_l_ddd;
  double AL_tt = delta_l*phir_l_dtt;
  double AV = delta_v * phir_v_d + 1;
  double AV_d = phir_v_d + delta_v*phir_v_dd;
  double AV_t = delta_v*phir_v_dt;
  double AV_dt = phir_v_dt + delta_v*phir_v_ddt;
  double AV_dd = 2*phir_v_dd + delta_v*phir_v_ddd;
  double AV_tt = delta_v*phir_v_dtt;

  /*
  From equal pressure:

    delta_l * AL = delta_v * AV

  we can differntiate with respect to tau and get:

    delta_l_t = delta_v_t * BV/BL + (CV - CL)/BL

  where BL, BV, CV and CL are defined below:
  */
  double BL = AL + delta_l * AL_d;
  double BV = AV + delta_v * AV_d;
  double CL = delta_l * AL_t;
  double CV = delta_v * AV_t;

  /*
  From the equlibrum condition we get

    AV - AL + ln(delta_v) - ln(delta_l) = phir_l - phir_v

  we can differentiate that with respect to tau and with a bunch of
  intermediates below and combining with the pressure equation derivatives
  find that:

    delta_v_t = H/G
    delta_l_t = delta_v_t * BV/BL + (CV - CL)/BL

  */
  double DL = phir_l_t + AL_t;
  double DV = phir_v_t + AV_t;
  double EL = AL_d + 1.0/delta_l + phir_l_d;
  double EV = AV_d + 1.0/delta_v + phir_v_d;
  double FL = 1.0/BL;
  double G = EV - EL*BV*FL;
  double H = DL - DV + EL*FL*(CV - CL);

  double delta_v_t = H/G;
  double delta_l_t = delta_v_t * BV/BL + (CV - CL)/BL;

  /*
  Now just differentiate the first detivatives again with respect to tau and
  we're done.
  */
  double BV_t = AV_t + delta_v * AV_dt;
  double BL_t = AL_t + delta_l * AL_dt;
  double BV_d = 2*AV_d + delta_v * AV_dd;
  double BL_d = 2*AL_d + delta_l * AL_dd;
  double EV_t = AV_dt + phir_v_dt;
  double EL_t = AL_dt + phir_l_dt;
  double EV_d = AV_dd - 1.0/delta_v/delta_v + phir_v_dd;
  double EL_d = AL_dd - 1.0/delta_l/delta_l + phir_l_dd;
  double CV_t = delta_v * AV_tt;
  double CL_t = delta_l * AL_tt;
  double CL_d = delta_l * AL_dt + AL_t;
  double CV_d = delta_v * AV_dt + AV_t;
  double DV_t = phir_v_tt + AV_tt;
  double DL_t = phir_l_tt + AL_tt;
  double DV_d = phir_v_dt + AV_dt;
  double DL_d = phir_l_dt + AL_dt;
  double FL_t = -BL_t/BL/BL;
  double FL_d = -BL_d/BL/BL;

  // Since up to here we've been using partials with respect to tau and delta
  // but really this is only a function of tau, we'll switch to total derivatives
  // with respect to tau.
  double G_total_t =
    EV_t +
    EV_d * delta_v_t -
    (EL_t + EL_d * delta_l_t) * BV * FL -
    EL * (BV_t + BV_d * delta_v_t) * FL -
    EL * BV * (FL_t + FL_d*delta_l_t);
  double H_total_t =
    DL_t +
    DL_d * delta_l_t -
    DV_t -
    DV_d * delta_v_t +
    (EL_t + EL_d * delta_l_t) * FL * (CV - CL) +
    (FL_t + FL_d * delta_l_t) * EL * (CV - CL) +
    EL * FL * (CV_t - CL_t + CV_d * delta_v_t - CL_d * delta_l_t);

  double delta_v_tt = H_total_t/G - H*G_total_t/G/G;
  double delta_l_tt =
    delta_v_tt * BV * FL +
    delta_v_t * (BV_t + BV_d * delta_v_t) * FL +
    delta_v*BV*(FL_t + FL_d*delta_l_t) +
    (FL_t + FL_d*delta_l_t)*(CV - CL) +
    FL * (CV_t - CL_t + CV_d*delta_v_t - CL_d*delta_l_t);

  /*
  While we're at it, may as well calculate Psat.  We'll use delta_v to do it.

  */
  f22_struct p_vec;
  pressure2(comp, delta_v, tau, &p_vec);

  double psat = p_vec.f;
  double psat_t =
    p_vec.f_2 +
    p_vec.f_1 * delta_v_t;
  double psat_tt =
    p_vec.f_22 +
    2*p_vec.f_12 * delta_v_t +
    p_vec.f_1 * delta_v_tt +
    p_vec.f_11 * delta_v_t * delta_v_t;

  delta_l_vec_ptr->f = delta_l;
  delta_l_vec_ptr->f_1 = delta_l_t;
  delta_l_vec_ptr->f_11 = delta_l_tt;
  delta_v_vec_ptr->f = delta_v;
  delta_v_vec_ptr->f_1 = delta_v_t;
  delta_v_vec_ptr->f_11 = delta_v_tt;
  p_vec_ptr->f = psat;
  p_vec_ptr->f_1 = psat_t;
  p_vec_ptr->f_11 = psat_tt;
  return n;
}

int cache_sat_delta_with_derivs(
    uint comp,
    double tau)
{

  if(memo_table_sat_delta_l.size() > MAX_MEMO_PROP){
     memo_table_sat_delta_l.clear();
     memo_table_sat_delta_v.clear();
     memo_table_sat_p.clear();
  }

  f12_struct *delta_l_vec_ptr = &memo_table_sat_delta_l[std::make_tuple(comp, tau)];
  f12_struct *delta_v_vec_ptr = &memo_table_sat_delta_v[std::make_tuple(comp, tau)];
  f12_struct *p_vec_ptr = &memo_table_sat_p[std::make_tuple(comp, tau)];  

  if (tau <= tau_c(comp)){
    tau = tau_c(comp);
    delta_l_vec_ptr = &memo_table_sat_delta_l[std::make_tuple(comp, tau)];
    delta_l_vec_ptr->f = delta_c(comp);
    delta_l_vec_ptr->f_1 = 0;
    delta_l_vec_ptr->f_11 = 0;
    delta_v_vec_ptr = &memo_table_sat_delta_v[std::make_tuple(comp, tau)];
    delta_v_vec_ptr->f = delta_c(comp);
    delta_v_vec_ptr->f_1 = 0;
    delta_v_vec_ptr->f_11 = 0;
    p_vec_ptr = &memo_table_sat_p[std::make_tuple(comp, tau)];
    p_vec_ptr->f = cdata[comp].Pc;
    p_vec_ptr->f_1 = 0;
    p_vec_ptr->f_11 = 0;
    return 0;
  }
  return sat_delta_with_derivs(comp, tau, delta_l_vec_ptr, delta_v_vec_ptr, p_vec_ptr);
}

f12_struct sat_delta_l(uint comp, double tau){
  if (isnan(tau)){
    f12_struct res;
    res.f = nan("");
    res.f_1 = nan("");
    res.f_11 = nan("");
    return res;
  }
  if (tau <= tau_c(comp)){
    tau = tau_c(comp);
  }
  try{
    return memo_table_sat_delta_l.at(std::make_tuple(comp, tau));
  }
  catch(std::out_of_range const&){
  }
  cache_sat_delta_with_derivs(comp, tau);
  return memo_table_sat_delta_l.at(std::make_tuple(comp, tau));
}

f12_struct sat_delta_v(uint comp, double tau){
  if (isnan(tau)){
    f12_struct res;
    res.f = nan("");
    res.f_1 = nan("");
    res.f_11 = nan("");
    return res;
  }
  if (tau <= tau_c(comp)){
    tau = tau_c(comp);
  }
  try{
    return memo_table_sat_delta_v.at(std::make_tuple(comp, tau));
  }
  catch(std::out_of_range const&){
  }
  cache_sat_delta_with_derivs(comp, tau);
  return memo_table_sat_delta_v.at(std::make_tuple(comp, tau));
}

f12_struct sat_p(uint comp, double tau){
  if (isnan(tau)){
    f12_struct res;
    res.f = nan("");
    res.f_1 = nan("");
    res.f_11 = nan("");
    std::cout << "got tau = nan" << std::endl;
    return res;
  }
  if (tau <= tau_c(comp)){
    tau = tau_c(comp);
  }
  try{
    return memo_table_sat_p.at(std::make_tuple(comp, tau));
  }
  catch(std::out_of_range const&){
  }
  cache_sat_delta_with_derivs(comp, tau);
  return memo_table_sat_p.at(std::make_tuple(comp, tau));
}


/*------------------------------------------------------------------------------
 Solve for tau_sat as a function of pressure.
------------------------------------------------------------------------------*/

class deltav_functor
{
private:
    uint _comp;
    double _pressure;
public:
    deltav_functor(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    double operator () (double tau) {
        double delta_l, delta_v;
        sat(this->_comp, tau, &delta_l, &delta_v);
        return (pressure(this->_comp, delta_v, tau) - this->_pressure)/cdata[this->_comp].Pc;
    }
};

class deltav_deriv_functor
{
private:
    uint _comp;
    double _pressure;
public:
    deltav_deriv_functor(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    std::tuple<double, double, double> operator () (double tau) {
        f12_struct dl, dv, p;
        sat_delta_with_derivs(this->_comp, tau, &dl, &dv, &p);
        return std::make_tuple(
          (p.f - this->_pressure) / cdata[this->_comp].Pc,
          p.f_1 / cdata[this->_comp].Pc,
          p.f_11 / cdata[this->_comp].Pc
        );
    }
};

// Solve for tau sat, calculate derivatives and cache results.
f12_struct sat_tau(uint comp, double pr){
  try{
    return memo_table_sat_tau.at(std::make_tuple(comp, pr));
  }
  catch(std::out_of_range const&){
  }
  double tau;
  deltav_functor fp(comp);
  deltav_deriv_functor fp2(comp);
  fp.set_pressure(pr);
  fp2.set_pressure(pr);

  if(memo_table_sat_tau.size() > MAX_MEMO_PROP){
     memo_table_sat_tau.clear();
  }
  f12_struct *res_ptr = &memo_table_sat_tau[std::make_tuple(comp, pr)];
  parameters_struct *dat = &cdata[comp];
  if(pr >= dat->Pc){
    tau = tau_c(comp);
    res_ptr->f = tau;
    res_ptr->f_1 = 0;
    res_ptr->f_11 = 0;
    return *res_ptr;
  }
  using namespace boost::math::tools;
  std::uintmax_t h_it_max=50;
  int digits = std::numeric_limits<double>::digits - 4;

  tau = halley_iterate(fp2, 1.01, 1.0, dat->T_star/dat->T_min, digits, h_it_max);
  f12_struct pres = sat_p(comp, tau);
  res_ptr->f = tau;
  res_ptr->f_1 = 1.0/pres.f_1; 
  res_ptr->f_11 = -res_ptr->f_1*res_ptr->f_1*res_ptr->f_1*pres.f_11;
  return *res_ptr;
}

// Solve for tau sat, calculate derivatives and cache results.
f12_struct sat_t(uint comp, double pr){
  f12_struct res, tau = sat_tau(comp, pr);
  double Tc = cdata[comp].T_star;
  res.f = Tc/tau.f;
  res.f_1 = -Tc*tau.f_1/tau.f/tau.f;
  res.f_11 = -Tc*tau.f_11/tau.f/tau.f + 2*Tc*tau.f_1*tau.f_1/tau.f/tau.f/tau.f;
  return res;
}

// Solve for tau sat, calculate derivatives and cache results.
f12_struct sat_p_t(uint comp, double t){
  double Tc = cdata[comp].T_star;
  f12_struct res, p_res = sat_p(comp, Tc/t);
  res.f = p_res.f;
  res.f_1 = -p_res.f_1 * Tc / t / t;
  res.f_11 = p_res.f_11 * Tc / t / t * Tc / t / t + 2*p_res.f_1 * Tc / t / t / t;
  return res;
}
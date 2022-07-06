/*------------------------------------------------------------------------------
Author: John Eslick
File sat.cpp

 This file contains the functions to solve the saturation curve.  The method
 is described in the following paper:

 Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
     State from Helmholtz Energy Equations of State." Journal of Thermal
     Science and Technology, 3(3), 442-451.
------------------------------------------------------------------------------*/
#include<unordered_map>
#include<vector>
#include<math.h>
#include<iostream>
#include<boost/functional/hash.hpp>
#include"phi.h"
#include"props.h"
#include"config.h"

static std::unordered_map<
  std::tuple<comp_enum, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double>>
> memo_table_sat_delta_l;

static std::unordered_map<
  std::tuple<comp_enum, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double>>
> memo_table_sat_delta_v;

static std::unordered_map<
  std::tuple<comp_enum, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double>>
> memo_table_sat_p;

inline double J(double delta, std::vector<double> *phi){
  // Term from Akasaka method for saturation state
  return delta*(1 + delta*phi->at(1));
}

inline double K(double delta, std::vector<double> *phi){
  // Term from Akasaka method for saturation state
  return delta*phi->at(1) + phi->at(0) + log(delta);
}

inline double J_delta(double delta, std::vector<double> *phi){
  // Derivative term from Akasaka method for saturation state
  return 1.0 + 2.0*delta*phi->at(1) + delta*delta*phi->at(2);
}

inline double K_delta(double delta, std::vector<double> *phi){
  // Derivative term from Akasaka method for saturation state
  return 2.0*phi->at(1) + delta*phi->at(2) + 1.0/delta;
}

unsigned int sat(comp_enum comp, double tau, double *delta_l, double *delta_v){
  double Kdiff=1, Jdiff=1, Jv, Jl, Kv, Kl, det, dJv, dJl, dKv, dKl;
  unsigned int n=0;
  std::vector<double> phir_v, phir_l;

  static const double tol_sat = 1e-11;
  static const double sat_gamma = 1.0;
  static const unsigned int max_iter = 100;

  if(tau - 1 < tol_sat){
    // So close to the critical point assume we are at the pritical point, if
    // over the critical point this will also use the critical density
    *delta_l = 1.0;
    *delta_v = 1.0;
  }
  else{
    // okay so you've decided to solve this thing
    *delta_l = 2.89; // DELTA_LIQ_SAT_GUESS;
    *delta_v = 0.0149; // DELTA_VAP_SAT_GUESS;
    while(n < max_iter){
      phi_real_for_sat(comp, *delta_v, tau, &phir_v);
      phi_real_for_sat(comp, *delta_l, tau, &phir_l);
      Jv = J(*delta_v, &phir_v);
      Jl = J(*delta_l, &phir_l);
      Kv = K(*delta_v, &phir_v);
      Kl = K(*delta_l, &phir_l);
      std::cout << Jv << " " << Jl << " " << Kv << " " << Kl << std::endl;
      Jdiff = Jv - Jl;
      Kdiff = Kv - Kl;
      if (fabs(Jdiff) < tol_sat && fabs(Kdiff) < tol_sat){
        break;
      }
      dJv = J_delta(*delta_v, &phir_v);
      dJl = J_delta(*delta_l, &phir_l);
      dKv = K_delta(*delta_v, &phir_v);
      dKl = K_delta(*delta_l, &phir_l);
      det = dJv*dKl - dJl*dKv;
      ++n; // Count iterations
      *delta_l += sat_gamma*(Kdiff*dJv - Jdiff*dKv)/det;
      *delta_v += sat_gamma*(Kdiff*dJl - Jdiff*dKl)/det;
    }
  }
  return n;
}

int cache_sat_delta_with_derivs(
    comp_enum comp,
    double tau)
{
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

  // Calcualte the real contribution to phi at the sat conditions
  std::vector<double> *phir_l_vec = phi_real(comp, delta_l, tau);
  std::vector<double> *phir_v_vec = phi_real(comp, delta_v, tau);

  // Get what I need of phi and it's derivatives.  This just makes the
  // expressions easier to read and write.
  double phir_l = phir_l_vec->at((unsigned int)deriv4_enum::f);
  double phir_l_d = phir_l_vec->at((unsigned int)deriv4_enum::f_d);
  double phir_l_dd = phir_l_vec->at((unsigned int)deriv4_enum::f_dd);
  double phir_l_ddd = phir_l_vec->at((unsigned int)deriv4_enum::f_ddd);
  double phir_l_t = phir_l_vec->at((unsigned int)deriv4_enum::f_t);
  double phir_l_dt = phir_l_vec->at((unsigned int)deriv4_enum::f_dt);
  double phir_l_ddt = phir_l_vec->at((unsigned int)deriv4_enum::f_ddt);
  double phir_l_tt = phir_l_vec->at((unsigned int)deriv4_enum::f_tt);
  double phir_l_dtt = phir_l_vec->at((unsigned int)deriv4_enum::f_dtt);
  double phir_v = phir_v_vec->at((unsigned int)deriv4_enum::f);
  double phir_v_d = phir_v_vec->at((unsigned int)deriv4_enum::f_d);
  double phir_v_dd = phir_v_vec->at((unsigned int)deriv4_enum::f_dd);
  double phir_v_ddd = phir_v_vec->at((unsigned int)deriv4_enum::f_ddd);
  double phir_v_t = phir_v_vec->at((unsigned int)deriv4_enum::f_t);
  double phir_v_dt = phir_v_vec->at((unsigned int)deriv4_enum::f_dt);
  double phir_v_ddt = phir_v_vec->at((unsigned int)deriv4_enum::f_ddt);
  double phir_v_tt = phir_v_vec->at((unsigned int)deriv4_enum::f_tt);
  double phir_v_dtt = phir_v_vec->at((unsigned int)deriv4_enum::f_dtt);

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
  std::vector<double> p_vec;
  pressure2(comp, delta_v, tau, &p_vec);

  double psat = p_vec.at((unsigned int)deriv2_enum::f);
  double psat_t =
    p_vec.at((unsigned int)deriv2_enum::f_t) +
    p_vec.at((unsigned int)deriv2_enum::f_d) * delta_v_t;
  double psat_tt =
    p_vec.at((unsigned int)deriv2_enum::f_tt) +
    2*p_vec.at((unsigned int)deriv2_enum::f_dt) * delta_v_t +
    p_vec.at((unsigned int)deriv2_enum::f_d) * delta_v_tt +
    p_vec.at((unsigned int)deriv2_enum::f_dd) * delta_v_t * delta_v_t;

  // Now memoize the result.  The memo is where the results are stored as well
  std::vector<double> *delta_l_vec_ptr;
  std::vector<double> *delta_v_vec_ptr;
  std::vector<double> *p_vec_ptr;


  if(memo_table_sat_delta_l.size() > MAX_MEMO_PROP){
     memo_table_sat_delta_l.clear();
     memo_table_sat_delta_v.clear();
  }
  delta_l_vec_ptr = &memo_table_sat_delta_l[std::make_tuple(comp, tau)];
  delta_l_vec_ptr->resize(3);
  delta_l_vec_ptr->at(0) = delta_l;
  delta_l_vec_ptr->at(1) = delta_l_t;
  delta_l_vec_ptr->at(2) = delta_l_tt;
  delta_v_vec_ptr = &memo_table_sat_delta_v[std::make_tuple(comp, tau)];
  delta_v_vec_ptr->resize(3);
  delta_v_vec_ptr->at(0) = delta_v;
  delta_v_vec_ptr->at(1) = delta_v_t;
  delta_v_vec_ptr->at(2) = delta_v_tt;
  p_vec_ptr = &memo_table_sat_p[std::make_tuple(comp, tau)];
  p_vec_ptr->resize(3);
  p_vec_ptr->at(0) = psat;
  p_vec_ptr->at(1) = psat_t;
  p_vec_ptr->at(2) = psat_tt;
  return n;
}

std::vector<double> *sat_delta_l(comp_enum comp, double tau){
  try{
    return &memo_table_sat_delta_l.at(std::make_tuple(comp, tau));
  }
  catch(std::out_of_range){
  }
  cache_sat_delta_with_derivs(comp, tau);
  return &memo_table_sat_delta_l.at(std::make_tuple(comp, tau));
}

std::vector<double> *sat_delta_v(comp_enum comp, double tau){
  try{
    return &memo_table_sat_delta_v.at(std::make_tuple(comp, tau));
  }
  catch(std::out_of_range){
  }
  cache_sat_delta_with_derivs(comp, tau);
  return &memo_table_sat_delta_v.at(std::make_tuple(comp, tau));
}

std::vector<double> *sat_p(comp_enum comp, double tau){
  try{
    return &memo_table_sat_p.at(std::make_tuple(comp, tau));
  }
  catch(std::out_of_range){
  }
  cache_sat_delta_with_derivs(comp, tau);
  return &memo_table_sat_p.at(std::make_tuple(comp, tau));
}

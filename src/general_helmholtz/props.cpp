#include"param.h"
#include"props.h"
#include"phi.h"

/*
This contains the basic properties as functions of delta and tau, where
delta = rho/rho_c and tau = Tc/T.

I used a lot of intermediate variables here to make it easier to check
and understand. The compiler optimization should optimize them out.
*/
enum property_cache_enum {
  p = 1,
  u = 2,
  s = 3,
  h = 4,
  cv = 5,
  cp = 6,
  w = 7,
};


double pressure(comp_enum comp, double delta, double tau){
  std::vector<double> *y = phi_real(comp, delta, tau);
  double phir_d = y->at((unsigned int)deriv4_enum::f_d);
  double u = 1 + delta*phir_d;
  return rhoc[comp]*Tc[comp]*R[comp]*delta/tau*u;
}

void pressure1(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_real(comp, delta, tau);
  double res[3];
  // while the intermediate varaibles look inefficent the copiler optimization
  // should fix it, and hopfully .
  double c = rhoc[comp]*Tc[comp]*R[comp];
  double phir_d = y->at((unsigned int)deriv4_enum::f_d);
  double phir_dd = y->at((unsigned int)deriv4_enum::f_dd);
  double phir_dt = y->at((unsigned int)deriv4_enum::f_dt);
  double u = 1 + delta*phir_d;
  double u_d = phir_d + delta*phir_dd;
  double u_t = delta*phir_dt;
  res[(unsigned int)deriv1_enum::f] = c*delta/tau*u; // pressure
  res[(unsigned int)deriv1_enum::f_d] = c*(u/tau + delta/tau*u_d);
  res[(unsigned int)deriv1_enum::f_t] = c*(-delta/tau/tau*u + delta/tau*u_t);
  out->assign(res, res+3);
}

void pressure2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_real(comp, delta, tau);
  double res[6];
  // while the intermediate varaibles look inefficent the copiler optimization
  // should fix it, and hopfully .
  double c = rhoc[comp]*Tc[comp]*R[comp];
  double phir_d = y->at((unsigned int)deriv4_enum::f_d);
  double phir_dd = y->at((unsigned int)deriv4_enum::f_dd);
  double phir_ddd = y->at((unsigned int)deriv4_enum::f_ddd);
  double phir_dt = y->at((unsigned int)deriv4_enum::f_dt);
  double phir_ddt = y->at((unsigned int)deriv4_enum::f_ddt);
  double phir_dtt = y->at((unsigned int)deriv4_enum::f_dtt);
  double u = 1 + delta*phir_d;
  double u_d = phir_d + delta*phir_dd;
  double u_dd = 2*phir_dd + delta*phir_ddd;
  double u_t = delta*phir_dt;
  double u_dt = phir_dt + delta*phir_ddt;
  double u_tt = delta*phir_dtt;
  res[(unsigned int)deriv2_enum::f] = c*delta/tau*u; // pressure
  res[(unsigned int)deriv2_enum::f_d] = c*(u/tau + delta/tau*u_d);
  res[(unsigned int)deriv2_enum::f_dd] = c*(2.0*u_d/tau + delta/tau*u_dd);
  res[(unsigned int)deriv2_enum::f_t] = c*(-delta/tau/tau*u + delta/tau*u_t);
  res[(unsigned int)deriv2_enum::f_dt] = c*(-u/tau/tau + u_t/tau + -delta/tau/tau*u_d + delta/tau*u_dt);
  res[(unsigned int)deriv2_enum::f_tt] = c*(2.0*delta/tau/tau/tau*u + -2.0*delta/tau/tau*u_t + delta/tau*u_tt);
  out->assign(res, res+6);
}

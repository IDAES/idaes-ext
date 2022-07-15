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

#include<unordered_map>
#include<boost/functional/hash.hpp>
#include"param.h"
#include"props.h"
#include"phi.h"

/*------------------------------------------------------------------------------
Author: John Eslick
File props.cpp

 This file contains basic property calculations as a function of delta and tau
 where delta = rho/rho_c and tau = T/T_c. The memo2_{prop} functions also
 memoize the property calculations with first and second derivatives.
------------------------------------------------------------------------------*/

// Memoization tables for property calculations

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_pressure2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_internal_energy2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_entropy2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_enthalpy2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_gibbs2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_helmholtz2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_isochoric_heat_capacity2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_isobaric_heat_capacity2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_speed_of_sound2;

double pressure(comp_enum comp, double delta, double tau){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  double phir_d = y->at((unsigned int)deriv4_enum::f_d);
  double u = 1 + delta*phir_d;
  return rhoc[comp]*Tc[comp]*R[comp]*delta/tau*u;
}

void pressure1(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  double res[3];
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
  std::vector<double> *y = phi_resi(comp, delta, tau);
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

void internal_energy2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = Tc[comp]*R[comp];
  double phii_t = yi->at((unsigned int)deriv4_enum::f_t);
  double phii_dt = yi->at((unsigned int)deriv4_enum::f_dt);
  double phii_tt = yi->at((unsigned int)deriv4_enum::f_tt);
  double phii_ddt = yi->at((unsigned int)deriv4_enum::f_ddt);
  double phii_dtt = yi->at((unsigned int)deriv4_enum::f_dtt);
  double phii_ttt = yi->at((unsigned int)deriv4_enum::f_ttt);
  double phir_t = yr->at((unsigned int)deriv4_enum::f_t);
  double phir_dt = yr->at((unsigned int)deriv4_enum::f_dt);
  double phir_tt = yr->at((unsigned int)deriv4_enum::f_tt);
  double phir_ddt = yr->at((unsigned int)deriv4_enum::f_ddt);
  double phir_dtt = yr->at((unsigned int)deriv4_enum::f_dtt);
  double phir_ttt = yr->at((unsigned int)deriv4_enum::f_ttt);

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  res[(unsigned int)deriv2_enum::f] = c*z;
  res[(unsigned int)deriv2_enum::f_d] = c*z_d;
  res[(unsigned int)deriv2_enum::f_dd] = c*z_dd;
  res[(unsigned int)deriv2_enum::f_t] = c*z_t;
  res[(unsigned int)deriv2_enum::f_dt] = c*z_dt;
  res[(unsigned int)deriv2_enum::f_tt] = c*z_tt;
  out->assign(res, res+6);
}

double entropy(comp_enum comp, double delta, double tau){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);
  double c = R[comp];
  double phii = yi->at((unsigned int)deriv4_enum::f);
  double phir = yr->at((unsigned int)deriv4_enum::f);
  double phii_t = yi->at((unsigned int)deriv4_enum::f_t);
  double phir_t = yr->at((unsigned int)deriv4_enum::f_t);
  double z = phii_t + phir_t;
  return c*(tau*z - phii - phir);
}

void entropy2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = R[comp];
  double phii = yi->at((unsigned int)deriv4_enum::f);
  double phii_d = yi->at((unsigned int)deriv4_enum::f_d);
  double phii_dd = yi->at((unsigned int)deriv4_enum::f_dd);
  double phii_t = yi->at((unsigned int)deriv4_enum::f_t);
  double phii_dt = yi->at((unsigned int)deriv4_enum::f_dt);
  double phii_tt = yi->at((unsigned int)deriv4_enum::f_tt);
  double phii_ddt = yi->at((unsigned int)deriv4_enum::f_ddt);
  double phii_dtt = yi->at((unsigned int)deriv4_enum::f_dtt);
  double phii_ttt = yi->at((unsigned int)deriv4_enum::f_ttt);
  double phir = yr->at((unsigned int)deriv4_enum::f);
  double phir_d = yr->at((unsigned int)deriv4_enum::f_d);
  double phir_dd = yr->at((unsigned int)deriv4_enum::f_dd);
  double phir_t = yr->at((unsigned int)deriv4_enum::f_t);
  double phir_dt = yr->at((unsigned int)deriv4_enum::f_dt);
  double phir_tt = yr->at((unsigned int)deriv4_enum::f_tt);
  double phir_ddt = yr->at((unsigned int)deriv4_enum::f_ddt);
  double phir_dtt = yr->at((unsigned int)deriv4_enum::f_dtt);
  double phir_ttt = yr->at((unsigned int)deriv4_enum::f_ttt);

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  res[(unsigned int)deriv2_enum::f] = c*(tau*z - phii - phir);
  res[(unsigned int)deriv2_enum::f_d] = c*(tau*z_d - phii_d - phir_d);
  res[(unsigned int)deriv2_enum::f_dd] = c*(tau*z_dd - phii_dd - phir_dd);
  res[(unsigned int)deriv2_enum::f_t] = c*(z + tau*z_t - phii_t - phir_t);
  res[(unsigned int)deriv2_enum::f_dt] = c*(z_d + tau*z_dt - phii_dt - phir_dt);
  res[(unsigned int)deriv2_enum::f_tt] = c*(2*z_t + tau*z_tt - phii_tt - phir_tt);
  out->assign(res, res+6);
}

double enthalpy(comp_enum comp, double delta, double tau){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);
  double c = R[comp]*Tc[comp];
  double phii_t = yi->at((unsigned int)deriv4_enum::f_t);
  double phir_t = yr->at((unsigned int)deriv4_enum::f_t);
  double phir_d = yr->at((unsigned int)deriv4_enum::f_d);
  double z = phii_t + phir_t;
  double x = delta/tau;
  return c*(1.0/tau + z + x * phir_d);
}

void enthalpy2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = R[comp]*Tc[comp];
  double phii_t = yi->at((unsigned int)deriv4_enum::f_t);
  double phii_dt = yi->at((unsigned int)deriv4_enum::f_dt);
  double phii_tt = yi->at((unsigned int)deriv4_enum::f_tt);
  double phii_ddt = yi->at((unsigned int)deriv4_enum::f_ddt);
  double phii_dtt = yi->at((unsigned int)deriv4_enum::f_dtt);
  double phii_ttt = yi->at((unsigned int)deriv4_enum::f_ttt);
  double phir_t = yr->at((unsigned int)deriv4_enum::f_t);
  double phir_d = yr->at((unsigned int)deriv4_enum::f_d);
  double phir_dd = yr->at((unsigned int)deriv4_enum::f_dd);
  double phir_ddd = yr->at((unsigned int)deriv4_enum::f_ddd);
  double phir_dt = yr->at((unsigned int)deriv4_enum::f_dt);
  double phir_tt = yr->at((unsigned int)deriv4_enum::f_tt);
  double phir_ddt = yr->at((unsigned int)deriv4_enum::f_ddt);
  double phir_dtt = yr->at((unsigned int)deriv4_enum::f_dtt);
  double phir_ttt = yr->at((unsigned int)deriv4_enum::f_ttt);

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  double x = delta/tau;
  double x_d = 1.0/tau;
  double x_t = -delta/tau/tau;
  double x_dt = -1.0/tau/tau;
  double x_tt = 2*delta/tau/tau/tau;

  res[(unsigned int)deriv2_enum::f] = c*(1.0/tau + z + x * phir_d);
  res[(unsigned int)deriv2_enum::f_d] = c*(z_d + x_d*phir_d + x*phir_dd);
  res[(unsigned int)deriv2_enum::f_dd] = c*(z_dd + 2*x_d*phir_dd + x*phir_ddd);
  res[(unsigned int)deriv2_enum::f_t] =
    c*(-1.0/tau/tau + z_t + x_t*phir_d + x*phir_dt);
  res[(unsigned int)deriv2_enum::f_dt] =
    c*(z_dt + x_dt*phir_d + x_d*phir_dt + x_t*phir_dd + x*phir_ddt);
  res[(unsigned int)deriv2_enum::f_tt] =
    c*(2.0/tau/tau/tau + z_tt + x_tt*phir_d + 2*x_t*phir_dt + x*phir_dtt);
  out->assign(res, res+6);
}




std::vector<double> *memo2_pressure(comp_enum comp, double delta, double tau){
  try{
    return &memo_table_pressure2.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_pressure2.size() > MAX_MEMO_PROP) memo_table_pressure2.clear();
  yvec_ptr = &memo_table_pressure2[std::make_tuple(comp, delta, tau)];
  pressure2(comp, delta, tau, yvec_ptr);
  return yvec_ptr;
}

std::vector<double> *memo2_internal_energy(comp_enum comp, double delta, double tau){
  try{
    return &memo_table_internal_energy2.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_internal_energy2.size() > MAX_MEMO_PROP) memo_table_internal_energy2.clear();
  yvec_ptr = &memo_table_internal_energy2[std::make_tuple(comp, delta, tau)];
  internal_energy2(comp, delta, tau, yvec_ptr);
  return yvec_ptr;
}

std::vector<double> *memo2_entropy(comp_enum comp, double delta, double tau){
  try{
    return &memo_table_entropy2.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_entropy2.size() > MAX_MEMO_PROP) memo_table_entropy2.clear();
  yvec_ptr = &memo_table_entropy2[std::make_tuple(comp, delta, tau)];
  entropy2(comp, delta, tau, yvec_ptr);
  return yvec_ptr;
}

std::vector<double> *memo2_enthalpy(comp_enum comp, double delta, double tau){
  try{
    return &memo_table_enthalpy2.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_enthalpy2.size() > MAX_MEMO_PROP) memo_table_enthalpy2.clear();
  yvec_ptr = &memo_table_enthalpy2[std::make_tuple(comp, delta, tau)];
  enthalpy2(comp, delta, tau, yvec_ptr);
  return yvec_ptr;
}

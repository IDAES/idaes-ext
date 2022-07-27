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

prop_memo_table2 memo_table_pressure2;
prop_memo_table2 memo_table_internal_energy2;
prop_memo_table2 memo_table_entropy2;
prop_memo_table2 memo_table_enthalpy2;
prop_memo_table2 memo_table_gibbs2;
prop_memo_table2 memo_table_helmholtz2;
prop_memo_table2 memo_table_isochoric_heat_capacity2;
prop_memo_table2 memo_table_isobaric_heat_capacity2;
prop_memo_table2 memo_table_speed_of_sound2;
prop_memo_table2 memo_table_phi_ideal2;
prop_memo_table2 memo_table_phi_resi2;
prop_memo_table2 memo_table_phi_ideal_d2;
prop_memo_table2 memo_table_phi_resi_d2;
prop_memo_table2 memo_table_phi_ideal_t2;
prop_memo_table2 memo_table_phi_resi_t2;
prop_memo_table2 memo_table_phi_ideal_dd2;
prop_memo_table2 memo_table_phi_resi_dd2;
prop_memo_table2 memo_table_phi_ideal_dt2;
prop_memo_table2 memo_table_phi_resi_dt2;
prop_memo_table2 memo_table_phi_ideal_tt2;
prop_memo_table2 memo_table_phi_resi_tt2;

static std::vector<double> nan_vec2 = {
  nan(""),
  nan(""),
  nan(""),
  nan(""),
  nan(""),
  nan("")
};

double pressure(comp_enum comp, double delta, double tau){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  double phir_d = y->at(f4_1);
  double u = 1 + delta*phir_d;
  return param::rhoc[comp]*param::Tc[comp]*param::R[comp]*delta/tau*u;
}

void pressure1(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  double res[3];
  double c = param::rhoc[comp]*param::Tc[comp]*param::R[comp];
  double phir_d = y->at(f4_1);
  double phir_dd = y->at(f4_11);
  double phir_dt = y->at(f4_12);
  double u = 1 + delta*phir_d;
  double u_d = phir_d + delta*phir_dd;
  double u_t = delta*phir_dt;
  res[f1] = c*delta/tau*u; // pressure
  res[f1_1] = c*(u/tau + delta/tau*u_d);
  res[f1_2] = c*(-delta/tau/tau*u + delta/tau*u_t);
  out->assign(res, res+3);
}

void pressure2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  double res[6];
  // while the intermediate varaibles look inefficent the copiler optimization
  // should fix it, and hopfully .
  double c = param::rhoc[comp]*param::Tc[comp]*param::R[comp];
  double phir_d = y->at(f4_1);
  double phir_dd = y->at(f4_11);
  double phir_ddd = y->at(f4_111);
  double phir_dt = y->at(f4_12);
  double phir_ddt = y->at(f4_112);
  double phir_dtt = y->at(f4_122);
  double u = 1 + delta*phir_d;
  double u_d = phir_d + delta*phir_dd;
  double u_dd = 2*phir_dd + delta*phir_ddd;
  double u_t = delta*phir_dt;
  double u_dt = phir_dt + delta*phir_ddt;
  double u_tt = delta*phir_dtt;
  res[f2] = c*delta/tau*u; // pressure
  res[f2_1] = c*(u/tau + delta/tau*u_d);
  res[f2_11] = c*(2.0*u_d/tau + delta/tau*u_dd);
  res[f2_2] = c*(-delta/tau/tau*u + delta/tau*u_t);
  res[f2_12] = c*(-u/tau/tau + u_t/tau + -delta/tau/tau*u_d + delta/tau*u_dt);
  res[f2_22] = c*(2.0*delta/tau/tau/tau*u + -2.0*delta/tau/tau*u_t + delta/tau*u_tt);
  out->assign(res, res+6);
}

double internal_energy(comp_enum comp, double delta, double tau){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);
  double c = param::Tc[comp]*param::R[comp];
  double phii_t = yi->at(f4_2);
  double phir_t = yr->at(f4_2);
  double z = phii_t + phir_t;
  return c*z;
}

void internal_energy2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::Tc[comp]*param::R[comp];
  double phii_t = yi->at(f4_2);
  double phii_dt = yi->at(f4_12);
  double phii_tt = yi->at(f4_22);
  double phii_ddt = yi->at(f4_112);
  double phii_dtt = yi->at(f4_122);
  double phii_ttt = yi->at(f4_222);
  double phir_t = yr->at(f4_2);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);
  double phir_ddt = yr->at(f4_112);
  double phir_dtt = yr->at(f4_122);
  double phir_ttt = yr->at(f4_222);

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  res[f2] = c*z;
  res[f2_1] = c*z_d;
  res[f2_11] = c*z_dd;
  res[f2_2] = c*z_t;
  res[f2_12] = c*z_dt;
  res[f2_22] = c*z_tt;
  out->assign(res, res+6);
}

double entropy(comp_enum comp, double delta, double tau){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);
  double c = param::R[comp];
  double phii = yi->at(f4);
  double phir = yr->at(f4);
  double phii_t = yi->at(f4_2);
  double phir_t = yr->at(f4_2);
  double z = phii_t + phir_t;
  return c*(tau*z - phii - phir);
}

void entropy2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp];
  double phii = yi->at(f4);
  double phii_d = yi->at(f4_1);
  double phii_dd = yi->at(f4_11);
  double phii_t = yi->at(f4_2);
  double phii_dt = yi->at(f4_12);
  double phii_tt = yi->at(f4_22);
  double phii_ddt = yi->at(f4_112);
  double phii_dtt = yi->at(f4_122);
  double phii_ttt = yi->at(f4_222);
  double phir = yr->at(f4);
  double phir_d = yr->at(f4_1);
  double phir_dd = yr->at(f4_11);
  double phir_t = yr->at(f4_2);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);
  double phir_ddt = yr->at(f4_112);
  double phir_dtt = yr->at(f4_122);
  double phir_ttt = yr->at(f4_222);

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  res[f2] = c*(tau*z - phii - phir);
  res[f2_1] = c*(tau*z_d - phii_d - phir_d);
  res[f2_11] = c*(tau*z_dd - phii_dd - phir_dd);
  res[f2_2] = c*(z + tau*z_t - phii_t - phir_t);
  res[f2_12] = c*(z_d + tau*z_dt - phii_dt - phir_dt);
  res[f2_22] = c*(2*z_t + tau*z_tt - phii_tt - phir_tt);
  out->assign(res, res+6);
}

double enthalpy(comp_enum comp, double delta, double tau){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);
  double c = param::R[comp]*param::Tc[comp];
  double phii_t = yi->at(f4_2);
  double phir_t = yr->at(f4_2);
  double phir_d = yr->at(f4_1);
  double z = phii_t + phir_t;
  double x = delta/tau;
  return c*(1.0/tau + z + x * phir_d);
}

void enthalpy2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp]*param::Tc[comp];
  double phii_t = yi->at(f4_2);
  double phii_dt = yi->at(f4_12);
  double phii_tt = yi->at(f4_22);
  double phii_ddt = yi->at(f4_112);
  double phii_dtt = yi->at(f4_122);
  double phii_ttt = yi->at(f4_222);
  double phir_t = yr->at(f4_2);
  double phir_d = yr->at(f4_1);
  double phir_dd = yr->at(f4_11);
  double phir_ddd = yr->at(f4_111);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);
  double phir_ddt = yr->at(f4_112);
  double phir_dtt = yr->at(f4_122);
  double phir_ttt = yr->at(f4_222);

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

  res[f2] = c*(1.0/tau + z + x * phir_d);
  res[f2_1] = c*(z_d + x_d*phir_d + x*phir_dd);
  res[f2_11] = c*(z_dd + 2*x_d*phir_dd + x*phir_ddd);
  res[f2_2] =
    c*(-1.0/tau/tau + z_t + x_t*phir_d + x*phir_dt);
  res[f2_12] =
    c*(z_dt + x_dt*phir_d + x_d*phir_dt + x_t*phir_dd + x*phir_ddt);
  res[f2_22] =
    c*(2.0/tau/tau/tau + z_tt + x_tt*phir_d + 2*x_t*phir_dt + x*phir_dtt);
  out->assign(res, res+6);
}

void gibbs2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp]*param::Tc[comp];
  double phii = yi->at(f4);
  double phii_d = yi->at(f4_1);
  double phii_dd = yi->at(f4_11);
  double phii_t = yi->at(f4_2);
  double phii_dt = yi->at(f4_12);
  double phii_tt = yi->at(f4_22);
  double phii_dtt = yi->at(f4_122);
  double phir = yr->at(f4);
  double phir_t = yr->at(f4_2);
  double phir_d = yr->at(f4_1);
  double phir_dd = yr->at(f4_11);
  double phir_ddd = yr->at(f4_111);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);
  double phir_ddt = yr->at(f4_112);
  double phir_dtt = yr->at(f4_122);

  res[f2] = c/tau*(1 + delta*phir_d + phii + phir);
  res[f2_1] = c/tau*(phir_d + delta*phir_dd + phii_d + phir_d);
  res[f2_11] = c/tau*(2*phir_dd + delta*phir_ddd + phii_dd + phir_dd);
  res[f2_2] = -1/tau*res[f2] + c/tau*(delta*phir_dt + phii_t + phir_t);
  res[f2_12] = -1/tau*res[f2_1] + c/tau*(phir_dt + delta*phir_ddt + phii_dt + phir_dt);
  res[f2_22] = 2*res[f2]/tau/tau - 2*c/tau/tau*(delta*phir_dt + phii_t + phir_t) +
    c/tau*(delta*phir_dtt + phii_tt + phir_tt);
  out->assign(res, res+6);
}

void helmholtz2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp]*param::Tc[comp];
  double phii = yi->at(f4);
  double phii_d = yi->at(f4_1);
  double phii_dd = yi->at(f4_11);
  double phii_t = yi->at(f4_2);
  double phii_dt = yi->at(f4_12);
  double phii_tt = yi->at(f4_22);
  double phir = yr->at(f4);
  double phir_t = yr->at(f4_2);
  double phir_d = yr->at(f4_1);
  double phir_dd = yr->at(f4_11);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);

  res[f2] = c/tau*(phii + phir);
  res[f2_1] = c/tau*(phii_d + phir_d);
  res[f2_11] = c/tau*(phii_dd + phir_dd);
  res[f2_2] = -1/tau*res[f2] + c/tau*(phii_t + phir_t);
  res[f2_12] = -1/tau*res[f2_1] + c/tau*(phii_dt + phir_dt);
  res[f2_22] = 1/tau/tau*res[f2] - 1/tau*res[f2_2] - c/tau/tau*(phii_t + phir_t) + c/tau*(phii_tt + phir_tt);
  out->assign(res, res+6);
}

void phi_ideal2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_ideal(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4);
  out->at(f2_1) = y->at(f4_1);
  out->at(f2_2) = y->at(f4_2);
  out->at(f2_11) = y->at(f4_11);
  out->at(f2_12) = y->at(f4_12);
  out->at(f2_22) = y->at(f4_22);
}

void phi_ideal_d2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_ideal(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_1);
  out->at(f2_1) = y->at(f4_11);
  out->at(f2_2) = y->at(f4_12);
  out->at(f2_11) = y->at(f4_111);
  out->at(f2_12) = y->at(f4_112);
  out->at(f2_22) = y->at(f4_122);
}

void phi_ideal_t2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_ideal(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_2);
  out->at(f2_1) = y->at(f4_12);
  out->at(f2_2) = y->at(f4_22);
  out->at(f2_11) = y->at(f4_112);
  out->at(f2_12) = y->at(f4_122);
  out->at(f2_22) = y->at(f4_222);
}

void phi_ideal_dd2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_ideal(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_11);
  out->at(f2_1) = y->at(f4_111);
  out->at(f2_2) = y->at(f4_112);
  out->at(f2_11) = y->at(f4_1111);
  out->at(f2_12) = y->at(f4_1112);
  out->at(f2_22) = y->at(f4_1122);
}

void phi_ideal_dt2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_ideal(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_12);
  out->at(f2_1) = y->at(f4_112);
  out->at(f2_2) = y->at(f4_122);
  out->at(f2_11) = y->at(f4_1112);
  out->at(f2_12) = y->at(f4_1122);
  out->at(f2_22) = y->at(f4_1222);
}

void phi_ideal_tt2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_ideal(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_22);
  out->at(f2_1) = y->at(f4_122);
  out->at(f2_2) = y->at(f4_222);
  out->at(f2_11) = y->at(f4_1122);
  out->at(f2_12) = y->at(f4_1222);
  out->at(f2_22) = y->at(f4_2222);
}

void phi_resi2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4);
  out->at(f2_1) = y->at(f4_1);
  out->at(f2_2) = y->at(f4_2);
  out->at(f2_11) = y->at(f4_11);
  out->at(f2_12) = y->at(f4_12);
  out->at(f2_22) = y->at(f4_22);
}

void phi_resi_d2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_1);
  out->at(f2_1) = y->at(f4_11);
  out->at(f2_2) = y->at(f4_12);
  out->at(f2_11) = y->at(f4_111);
  out->at(f2_12) = y->at(f4_112);
  out->at(f2_22) = y->at(f4_122);
}

void phi_resi_t2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_2);
  out->at(f2_1) = y->at(f4_12);
  out->at(f2_2) = y->at(f4_22);
  out->at(f2_11) = y->at(f4_112);
  out->at(f2_12) = y->at(f4_122);
  out->at(f2_22) = y->at(f4_222);
}

void phi_resi_dd2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_11);
  out->at(f2_1) = y->at(f4_111);
  out->at(f2_2) = y->at(f4_112);
  out->at(f2_11) = y->at(f4_1111);
  out->at(f2_12) = y->at(f4_1112);
  out->at(f2_22) = y->at(f4_1122);
}

void phi_resi_dt2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_12);
  out->at(f2_1) = y->at(f4_112);
  out->at(f2_2) = y->at(f4_122);
  out->at(f2_11) = y->at(f4_1112);
  out->at(f2_12) = y->at(f4_1122);
  out->at(f2_22) = y->at(f4_1222);
}

void phi_resi_tt2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *y = phi_resi(comp, delta, tau);
  out->resize(6);
  out->at(f2) = y->at(f4_22);
  out->at(f2_1) = y->at(f4_122);
  out->at(f2_2) = y->at(f4_222);
  out->at(f2_11) = y->at(f4_1122);
  out->at(f2_12) = y->at(f4_1222);
  out->at(f2_22) = y->at(f4_2222);
}

MEMO2_FUNCTION(memo2_pressure, pressure2, memo_table_pressure2)
MEMO2_FUNCTION(memo2_internal_energy, internal_energy2, memo_table_internal_energy2)
MEMO2_FUNCTION(memo2_entropy, entropy2, memo_table_entropy2)
MEMO2_FUNCTION(memo2_enthalpy, enthalpy2, memo_table_enthalpy2)
MEMO2_FUNCTION(memo2_gibbs, gibbs2, memo_table_gibbs2)
MEMO2_FUNCTION(memo2_helmholtz, helmholtz2, memo_table_helmholtz2)

MEMO2_FUNCTION(memo2_phi_ideal, phi_ideal2, memo_table_phi_ideal2)
MEMO2_FUNCTION(memo2_phi_ideal_d, phi_ideal_d2, memo_table_phi_ideal_d2)
MEMO2_FUNCTION(memo2_phi_ideal_dd, phi_ideal_dd2, memo_table_phi_ideal_dd2)
MEMO2_FUNCTION(memo2_phi_ideal_t, phi_ideal_t2, memo_table_phi_ideal_t2)
MEMO2_FUNCTION(memo2_phi_ideal_dt, phi_ideal_dt2, memo_table_phi_ideal_dt2)
MEMO2_FUNCTION(memo2_phi_ideal_tt, phi_ideal_tt2, memo_table_phi_ideal_tt2)

MEMO2_FUNCTION(memo2_phi_resi, phi_resi2, memo_table_phi_resi2)
MEMO2_FUNCTION(memo2_phi_resi_d, phi_resi_d2, memo_table_phi_resi_d2)
MEMO2_FUNCTION(memo2_phi_resi_dd, phi_resi_dd2, memo_table_phi_resi_dd2)
MEMO2_FUNCTION(memo2_phi_resi_t, phi_resi_t2, memo_table_phi_resi_t2)
MEMO2_FUNCTION(memo2_phi_resi_dt, phi_resi_dt2, memo_table_phi_resi_dt2)
MEMO2_FUNCTION(memo2_phi_resi_tt, phi_resi_tt2, memo_table_phi_resi_tt2)

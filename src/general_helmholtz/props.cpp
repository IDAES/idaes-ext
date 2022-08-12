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
  res[f2_2] =  c*(z +   tau*z_t -  phii_t -  phir_t);
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

void isochoric_heat_capacity2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp];
  double phii_tt = yi->at(f4_22);
  double phii_dtt = yi->at(f4_122);
  double phii_ddtt = yi->at(f4_1122);
  double phii_dttt = yi->at(f4_1222);
  double phii_ttt = yi->at(f4_222);
  double phii_tttt = yi->at(f4_2222);
  double phir_tt = yr->at(f4_22);
  double phir_dtt = yr->at(f4_122);
  double phir_ddtt = yr->at(f4_1122);
  double phir_dttt = yr->at(f4_1222);
  double phir_ttt = yr->at(f4_222);
  double phir_tttt = yr->at(f4_2222);

  res[f2] = -c*tau*tau*(phii_tt + phir_tt);
  res[f2_1] = -c*tau*tau*(phii_dtt + phir_dtt);
  res[f2_11] = -c*tau*tau*(phii_ddtt + phir_ddtt);
  res[f2_2] = -c*(2*tau*(phii_tt + phir_tt) + tau*tau*(phii_ttt + phir_ttt));
  res[f2_12] = -c*(2*tau*(phii_dtt + phir_dtt) + tau*tau*(phii_dttt + phir_dttt));
  res[f2_22] = -c*(2*(phii_tt + phir_tt) + 4*tau*(phii_ttt + phir_ttt) + tau*tau*(phii_tttt + phir_tttt));
  out->assign(res, res+6);
}

void isobaric_heat_capacity2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp];
  double phii_tt = yi->at(f4_22);
  double phii_dtt = yi->at(f4_122);
  double phii_ddtt = yi->at(f4_1122);
  double phii_dttt = yi->at(f4_1222);
  double phii_ttt = yi->at(f4_222);
  double phii_tttt = yi->at(f4_2222);
  double phir_d = yr->at(f4_1);
  double phir_dd = yr->at(f4_11);
  double phir_ddd = yr->at(f4_111);
  double phir_dddd = yr->at(f4_1111);
  double phir_dddt = yr->at(f4_1112);
  double phir_ddt = yr->at(f4_112);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);
  double phir_dtt = yr->at(f4_122);
  double phir_ddtt = yr->at(f4_1122);
  double phir_dttt = yr->at(f4_1222);
  double phir_ttt = yr->at(f4_222);
  double phir_tttt = yr->at(f4_2222);

  double y = -tau*tau*(phii_tt + phir_tt);
  double y_d = -tau*tau*(phii_dtt + phir_dtt);
  double y_dd = -tau*tau*(phii_ddtt + phir_ddtt);
  double y_t = -1*(2*tau*(phii_tt + phir_tt) + tau*tau*(phii_ttt + phir_ttt));
  double y_dt = -1*(2*tau*(phii_dtt + phir_dtt) + tau*tau*(phii_dttt + phir_dttt));
  double y_tt = -1*(2*(phii_tt + phir_tt) + 4*tau*(phii_ttt + phir_ttt) + tau*tau*(phii_tttt + phir_tttt));

  double x = 1 + delta*phir_d - delta*tau*phir_dt;
  double x_d = phir_d + delta*phir_dd - tau*phir_dt - delta*tau*phir_ddt;
  double x_t = -delta*tau*phir_dtt;
  double x_dd = 2*phir_dd + delta*phir_ddd - 2*tau*phir_ddt - delta*tau*phir_dddt;
  double x_dt = -tau*phir_dtt - delta*tau*phir_ddtt;
  double x_tt = -delta*phir_dtt - delta*tau*phir_dttt;

  double z = 1 + 2*delta*phir_d + delta*delta*phir_dd;
  double z_d = 2*phir_d + 4*delta*phir_dd + delta*delta*phir_ddd;
  double z_t = 2*delta*phir_dt + delta*delta*phir_ddt;
  double z_dd = 6*phir_dd + 6*delta*phir_ddd + delta*delta*phir_dddd;
  double z_dt = 2*phir_dt + 4*delta*phir_ddt + delta*delta*phir_dddt;
  double z_tt = 2*delta*phir_dtt + delta*delta*phir_ddtt;

  res[f2] = c*(y + x*x/z);
  res[f2_1] = c*(y_d + 2*x/z*x_d - x*x/z/z*z_d);
  res[f2_11] = c*(y_dd + (2/z*x_d - 2*x/z/z*z_d)*x_d + 2*x/z*x_dd + (-2*x/z/z*x_d + 2*x*x/z/z/z*z_d)*z_d - x*x/z/z*z_dd);
  res[f2_2] = c*(y_t + 2*x/z*x_t - x*x/z/z*z_t);
  res[f2_12] = c*(y_dt + (2/z*x_t - 2*x/z/z*z_t)*x_d + 2*x/z*x_dt + (-2*x/z/z*x_t + 2*x*x/z/z/z*z_t)*z_d - x*x/z/z*z_dt);
  res[f2_22] = c*(y_tt + (2/z*x_t - 2*x/z/z*z_t)*x_t + 2*x/z*x_tt + (-2*x/z/z*x_t + 2*x*x/z/z/z*z_t)*z_t - x*x/z/z*z_tt);
  out->assign(res, res+6);
}

void speed_of_sound2(comp_enum comp, double delta, double tau, std::vector<double> *out){
  std::vector<double> *yr = phi_resi(comp, delta, tau);
  std::vector<double> *yi = phi_ideal(comp, delta, tau);

  double res[6];
  double c = param::R[comp]*param::Tc[comp]*1000; // the 100 is because when the units shae out you get w^2 [=] km*m/s^2 so convert to m^2/s^2
  double phii_tt = yi->at(f4_22);
  double phii_dtt = yi->at(f4_122);
  double phii_ddtt = yi->at(f4_1122);
  double phii_dttt = yi->at(f4_1222);
  double phii_ttt = yi->at(f4_222);
  double phii_tttt = yi->at(f4_2222);
  double phir_d = yr->at(f4_1);
  double phir_dd = yr->at(f4_11);
  double phir_ddd = yr->at(f4_111);
  double phir_dddd = yr->at(f4_1111);
  double phir_dddt = yr->at(f4_1112);
  double phir_ddt = yr->at(f4_112);
  double phir_dt = yr->at(f4_12);
  double phir_tt = yr->at(f4_22);
  double phir_dtt = yr->at(f4_122);
  double phir_ddtt = yr->at(f4_1122);
  double phir_dttt = yr->at(f4_1222);
  double phir_ttt = yr->at(f4_222);
  double phir_tttt = yr->at(f4_2222);

  double y = -tau*tau*(phii_tt + phir_tt);
  double y_d = -tau*tau*(phii_dtt + phir_dtt);
  double y_dd = -tau*tau*(phii_ddtt + phir_ddtt);
  double y_t = -1*(2*tau*(phii_tt + phir_tt) + tau*tau*(phii_ttt + phir_ttt));
  double y_dt = -1*(2*tau*(phii_dtt + phir_dtt) + tau*tau*(phii_dttt + phir_dttt));
  double y_tt = -1*(2*(phii_tt + phir_tt) + 4*tau*(phii_ttt + phir_ttt) + tau*tau*(phii_tttt + phir_tttt));

  double x = 1 + delta*phir_d - delta*tau*phir_dt;
  double x_d = phir_d + delta*phir_dd - tau*phir_dt - delta*tau*phir_ddt;
  double x_t = -delta*tau*phir_dtt;
  double x_dd = 2*phir_dd + delta*phir_ddd - 2*tau*phir_ddt - delta*tau*phir_dddt;
  double x_dt = -tau*phir_dtt - delta*tau*phir_ddtt;
  double x_tt = -delta*phir_dtt - delta*tau*phir_dttt;

  double z = 1 + 2*delta*phir_d + delta*delta*phir_dd;
  double z_d = 2*phir_d + 4*delta*phir_dd + delta*delta*phir_ddd;
  double z_t = 2*delta*phir_dt + delta*delta*phir_ddt;
  double z_dd = 6*phir_dd + 6*delta*phir_ddd + delta*delta*phir_dddd;
  double z_dt = 2*phir_dt + 4*delta*phir_ddt + delta*delta*phir_dddt;
  double z_tt = 2*delta*phir_dtt + delta*delta*phir_ddtt;

  double w2 = c/tau*(z + x*x/y);
  double w2_d = c/tau*(z_d + 2*x/y*x_d - x*x/y/y*y_d);
  double w2_t = -w2/tau + c/tau*(z_t + 2*x/y*x_t - x*x/y/y*y_t);
  double w2_dd = c/tau*(z_dd + 2/y*x_d*x_d - 4*x/y/y*x_d*y_d + 2*x/y*x_dd + 2*x*x/y/y/y*y_d*y_d - x*x/y/y*y_dd);
  double w2_dt = -w2_d/tau   + c/tau*(z_dt + 2/y*x_d*x_t - 2*x/y/y*x_t*y_d - 2*x/y/y*x_d*y_t + 2*x/y*x_dt + 2*x*x/y/y/y*y_d*y_t - x*x/y/y*y_dt);
  double w2_tt = -2*w2_t/tau + c/tau*(z_tt + 2/y*x_t*x_t - 4*x/y/y*x_t*y_t + 2*x/y*x_tt + 2*x*x/y/y/y*y_t*y_t - x*x/y/y*y_tt);

  res[f2] = sqrt(w2);
  res[f2_1] = 0.5*pow(w2, -0.5) * w2_d;
  res[f2_11] = 0.5*pow(w2, -0.5) * w2_dd - 0.25*pow(w2, -1.5) * w2_d*w2_d;
  res[f2_2] = 0.5*pow(w2, -0.5) * w2_t;
  res[f2_12] = 0.5*pow(w2, -0.5) * w2_dt - 0.25*pow(w2, -1.5) * w2_d*w2_t;
  res[f2_22] = 0.5*pow(w2, -0.5) * w2_tt - 0.25*pow(w2, -1.5) * w2_t*w2_t;
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
MEMO2_FUNCTION(memo2_isochoric_heat_capacity, isochoric_heat_capacity2, memo_table_isochoric_heat_capacity2)
MEMO2_FUNCTION(memo2_isobaric_heat_capacity, isobaric_heat_capacity2, memo_table_isobaric_heat_capacity2)
MEMO2_FUNCTION(memo2_speed_of_sound, speed_of_sound2, memo_table_speed_of_sound2)


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

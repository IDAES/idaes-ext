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
#include<iostream>

/*------------------------------------------------------------------------------
Author: John Eslick
File props.cpp

 This file contains basic property calculations as a function of delta and tau
 where delta = rho/rho_c and tau = T/T_c. The memo2_{prop} functions also
 memoize the property calculations with first and second derivatives.
------------------------------------------------------------------------------*/

// Memoization tables for property calculations

prop_memo_table22 memo_table_pressure2;
prop_memo_table22 memo_table_internal_energy2;
prop_memo_table22 memo_table_entropy2;
prop_memo_table22 memo_table_enthalpy2;
prop_memo_table22 memo_table_gibbs2;
prop_memo_table22 memo_table_helmholtz2;
prop_memo_table22 memo_table_isochoric_heat_capacity2;
prop_memo_table22 memo_table_isobaric_heat_capacity2;
prop_memo_table22 memo_table_speed_of_sound2;

prop_memo_table22 memo_table_phi_ideal2;
prop_memo_table22 memo_table_phi_resi2;
prop_memo_table22 memo_table_phi_ideal_d2;
prop_memo_table22 memo_table_phi_resi_d2;
prop_memo_table22 memo_table_phi_ideal_t2;
prop_memo_table22 memo_table_phi_resi_t2;
prop_memo_table22 memo_table_phi_ideal_dd2;
prop_memo_table22 memo_table_phi_resi_dd2;
prop_memo_table22 memo_table_phi_ideal_dt2;
prop_memo_table22 memo_table_phi_resi_dt2;
prop_memo_table22 memo_table_phi_ideal_tt2;
prop_memo_table22 memo_table_phi_resi_tt2;


double pressure(uint comp, double delta, double tau){
  f23_struct phir = phi_resi(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  return dat->rho_star*dat->T_star*dat->R*delta/tau*(1 + delta*phir.f_1);
}

void pressure1(uint comp, double delta, double tau, f21_struct *out){
  f23_struct y = phi_resi(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->rho_star*dat->T_star*dat->R;
  double phir_d = y.f_1;
  double phir_dd = y.f_11;
  double phir_dt = y.f_12;
  double u = 1 + delta*phir_d;
  double u_d = phir_d + delta*phir_dd;
  double u_t = delta*phir_dt;
  out->f = c*delta/tau*u; // pressure
  out->f_1 = c*(u/tau + delta/tau*u_d);
  out->f_2 = c*(-delta/tau/tau*u + delta/tau*u_t);
}

void pressure2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_resi(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  // while the intermediate varaibles look inefficent the copiler optimization
  // should fix it, and hopfully .
  double c = dat->rho_star*dat->T_star*dat->R;
  double phir_d = y.f_1;
  double phir_dd = y.f_11;
  double phir_ddd = y.f_111;
  double phir_dt = y.f_12;
  double phir_ddt = y.f_112;
  double phir_dtt = y.f_122;
  double u = 1 + delta*phir_d;
  double u_d = phir_d + delta*phir_dd;
  double u_dd = 2*phir_dd + delta*phir_ddd;
  double u_t = delta*phir_dt;
  double u_dt = phir_dt + delta*phir_ddt;
  double u_tt = delta*phir_dtt;
  out->f = c*delta/tau*u; // pressure
  out->f_1 = c*(u/tau + delta/tau*u_d);
  out->f_11 = c*(2.0*u_d/tau + delta/tau*u_dd);
  out->f_2 = c*(-delta/tau/tau*u + delta/tau*u_t);
  out->f_12 = c*(-u/tau/tau + u_t/tau + -delta/tau/tau*u_d + delta/tau*u_dt);
  out->f_22 = c*(2.0*delta/tau/tau/tau*u + -2.0*delta/tau/tau*u_t + delta/tau*u_tt);
}

double internal_energy(uint comp, double delta, double tau){
  f23_struct phir = phi_resi(comp, delta, tau);
  f23_struct phii = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->T_star*dat->R;
  double z = phii.f_2 + phir.f_2;
  return c*z;
}

void internal_energy2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->T_star*dat->R;
  double phii_t = yi.f_2;
  double phii_dt = yi.f_12;
  double phii_tt = yi.f_22;
  double phii_ddt = yi.f_112;
  double phii_dtt = yi.f_122;
  double phii_ttt = yi.f_222;
  double phir_t = yr.f_2;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;
  double phir_ddt = yr.f_112;
  double phir_dtt = yr.f_122;
  double phir_ttt = yr.f_222;

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  out->f = c*z;
  out->f_1 = c*z_d;
  out->f_11 = c*z_dd;
  out->f_2 = c*z_t;
  out->f_12 = c*z_dt;
  out->f_22 = c*z_tt;
}

double entropy(uint comp, double delta, double tau){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R;
  double phii = yi.f;
  double phir = yr.f;
  double phii_t = yi.f_2;
  double phir_t = yr.f_2;
  double z = phii_t + phir_t;
  return c*(tau*z - phii - phir);
}

void entropy2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R;
  double phii = yi.f;
  double phii_d = yi.f_1;
  double phii_dd = yi.f_11;
  double phii_t = yi.f_2;
  double phii_dt = yi.f_12;
  double phii_tt = yi.f_22;
  double phii_ddt = yi.f_112;
  double phii_dtt = yi.f_122;
  double phii_ttt = yi.f_222;
  double phir = yr.f;
  double phir_d = yr.f_1;
  double phir_dd = yr.f_11;
  double phir_t = yr.f_2;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;
  double phir_ddt = yr.f_112;
  double phir_dtt = yr.f_122;
  double phir_ttt = yr.f_222;

  double z = phii_t + phir_t;
  double z_d = phii_dt + phir_dt;
  double z_dd = phii_ddt + phir_ddt;
  double z_t = phii_tt + phir_tt;
  double z_dt = phii_dtt + phir_dtt;
  double z_tt = phii_ttt + phir_ttt;
  out->f = c*(tau*z - phii - phir);
  out->f_1 = c*(tau*z_d - phii_d - phir_d);
  out->f_11 = c*(tau*z_dd - phii_dd - phir_dd);
  out->f_2 =  c*(z +   tau*z_t -  phii_t -  phir_t);
  out->f_12 = c*(z_d + tau*z_dt - phii_dt - phir_dt);
  out->f_22 = c*(2*z_t + tau*z_tt - phii_tt - phir_tt);
}

double enthalpy(uint comp, double delta, double tau){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R*dat->T_star;
  double phii_t = yi.f_2;
  double phir_t = yr.f_2;
  double phir_d = yr.f_1;
  double z = phii_t + phir_t;
  double x = delta/tau;
  return c*(1.0/tau + z + x * phir_d);
}

void enthalpy2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R*dat->T_star;
  double phii_t = yi.f_2;
  double phii_dt = yi.f_12;
  double phii_tt = yi.f_22;
  double phii_ddt = yi.f_112;
  double phii_dtt = yi.f_122;
  double phii_ttt = yi.f_222;
  double phir_t = yr.f_2;
  double phir_d = yr.f_1;
  double phir_dd = yr.f_11;
  double phir_ddd = yr.f_111;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;
  double phir_ddt = yr.f_112;
  double phir_dtt = yr.f_122;
  double phir_ttt = yr.f_222;

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

  out->f = c*(1.0/tau + z + x * phir_d);
  out->f_1 = c*(z_d + x_d*phir_d + x*phir_dd);
  out->f_11 = c*(z_dd + 2*x_d*phir_dd + x*phir_ddd);
  out->f_2 = c*(-1.0/tau/tau + z_t + x_t*phir_d + x*phir_dt);
  out->f_12 = c*(z_dt + x_dt*phir_d + x_d*phir_dt + x_t*phir_dd + x*phir_ddt);
  out->f_22 = c*(2.0/tau/tau/tau + z_tt + x_tt*phir_d + 2*x_t*phir_dt + x*phir_dtt);
}

void gibbs2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R*dat->T_star;
  double phii = yi.f;
  double phii_d = yi.f_1;
  double phii_dd = yi.f_11;
  double phii_t = yi.f_2;
  double phii_dt = yi.f_12;
  double phii_tt = yi.f_22;
  double phir = yr.f;
  double phir_t = yr.f_2;
  double phir_d = yr.f_1;
  double phir_dd = yr.f_11;
  double phir_ddd = yr.f_111;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;
  double phir_ddt = yr.f_112;
  double phir_dtt = yr.f_122;

  out->f = c/tau*(1 + delta*phir_d + phii + phir);
  out->f_1 = c/tau*(phir_d + delta*phir_dd + phii_d + phir_d);
  out->f_11 = c/tau*(2*phir_dd + delta*phir_ddd + phii_dd + phir_dd);
  out->f_2 = -1/tau*out->f + c/tau*(delta*phir_dt + phii_t + phir_t);
  out->f_12 = -1/tau*out->f_1 + c/tau*(phir_dt + delta*phir_ddt + phii_dt + phir_dt);
  out->f_22 = 2*out->f/tau/tau - 2*c/tau/tau*(delta*phir_dt + phii_t + phir_t) +
    c/tau*(delta*phir_dtt + phii_tt + phir_tt);
}

void helmholtz2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct yr = phi_resi(comp, delta, tau);
  f23_struct yi = phi_ideal(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R*dat->T_star;
  double phii = yi.f;
  double phii_d = yi.f_1;
  double phii_dd = yi.f_11;
  double phii_t = yi.f_2;
  double phii_dt = yi.f_12;
  double phii_tt = yi.f_22;
  double phir = yr.f;
  double phir_t = yr.f_2;
  double phir_d = yr.f_1;
  double phir_dd = yr.f_11;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;

  out->f = c/tau*(phii + phir);
  out->f_1 = c/tau*(phii_d + phir_d);
  out->f_11 = c/tau*(phii_dd + phir_dd);
  out->f_2 = -1/tau*out->f + c/tau*(phii_t + phir_t);
  out->f_12 = -1/tau*out->f_1 + c/tau*(phii_dt + phir_dt);
  out->f_22 = 1/tau/tau*out->f - 1/tau*out->f_2 - c/tau/tau*(phii_t + phir_t) + c/tau*(phii_tt + phir_tt);
}

void isochoric_heat_capacity2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct yr = phi_resi4(comp, delta, tau);
  f24_struct yi = phi_ideal4(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R;
  double phii_tt = yi.f_22;
  double phii_dtt = yi.f_122;
  double phii_ddtt = yi.f_1122;
  double phii_dttt = yi.f_1222;
  double phii_ttt = yi.f_222;
  double phii_tttt = yi.f_2222;
  double phir_tt = yr.f_22;
  double phir_dtt = yr.f_122;
  double phir_ddtt = yr.f_1122;
  double phir_dttt = yr.f_1222;
  double phir_ttt = yr.f_222;
  double phir_tttt = yr.f_2222;

  out->f = -c*tau*tau*(phii_tt + phir_tt);
  out->f_1 = -c*tau*tau*(phii_dtt + phir_dtt);
  out->f_11 = -c*tau*tau*(phii_ddtt + phir_ddtt);
  out->f_2 = -c*(2*tau*(phii_tt + phir_tt) + tau*tau*(phii_ttt + phir_ttt));
  out->f_12 = -c*(2*tau*(phii_dtt + phir_dtt) + tau*tau*(phii_dttt + phir_dttt));
  out->f_22 = -c*(2*(phii_tt + phir_tt) + 4*tau*(phii_ttt + phir_ttt) + tau*tau*(phii_tttt + phir_tttt));
}

void isobaric_heat_capacity2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct yr = phi_resi4(comp, delta, tau);
  f24_struct yi = phi_ideal4(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];
  double c = dat->R;
  double phii_tt = yi.f_22;
  double phii_dtt = yi.f_122;
  double phii_ddtt = yi.f_1122;
  double phii_dttt = yi.f_1222;
  double phii_ttt = yi.f_222;
  double phii_tttt = yi.f_2222;
  double phir_d = yr.f_1;
  double phir_dd = yr.f_11;
  double phir_ddd = yr.f_111;
  double phir_dddd = yr.f_1111;
  double phir_dddt = yr.f_1112;
  double phir_ddt = yr.f_112;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;
  double phir_dtt = yr.f_122;
  double phir_ddtt = yr.f_1122;
  double phir_dttt = yr.f_1222;
  double phir_ttt = yr.f_222;
  double phir_tttt = yr.f_2222;

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

  out->f = c*(y + x*x/z);
  out->f_1 = c*(y_d + 2*x/z*x_d - x*x/z/z*z_d);
  out->f_11 = c*(y_dd + (2/z*x_d - 2*x/z/z*z_d)*x_d + 2*x/z*x_dd + (-2*x/z/z*x_d + 2*x*x/z/z/z*z_d)*z_d - x*x/z/z*z_dd);
  out->f_2 = c*(y_t + 2*x/z*x_t - x*x/z/z*z_t);
  out->f_12 = c*(y_dt + (2/z*x_t - 2*x/z/z*z_t)*x_d + 2*x/z*x_dt + (-2*x/z/z*x_t + 2*x*x/z/z/z*z_t)*z_d - x*x/z/z*z_dt);
  out->f_22 = c*(y_tt + (2/z*x_t - 2*x/z/z*z_t)*x_t + 2*x/z*x_tt + (-2*x/z/z*x_t + 2*x*x/z/z/z*z_t)*z_t - x*x/z/z*z_tt);
}

void speed_of_sound2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct yr = phi_resi4(comp, delta, tau);
  f24_struct yi = phi_ideal4(comp, delta, tau);
  parameters_struct *dat = &cdata[comp];

  double c = dat->R*dat->T_star*1000; // the 1000 is because when the units shake out you get w^2 [=] km*m/s^2 so convert to m^2/s^2
  double phii_tt = yi.f_22;
  double phii_dtt = yi.f_122;
  double phii_ddtt = yi.f_1122;
  double phii_dttt = yi.f_1222;
  double phii_ttt = yi.f_222;
  double phii_tttt = yi.f_2222;
  double phir_d = yr.f_1;
  double phir_dd = yr.f_11;
  double phir_ddd = yr.f_111;
  double phir_dddd = yr.f_1111;
  double phir_dddt = yr.f_1112;
  double phir_ddt = yr.f_112;
  double phir_dt = yr.f_12;
  double phir_tt = yr.f_22;
  double phir_dtt = yr.f_122;
  double phir_ddtt = yr.f_1122;
  double phir_dttt = yr.f_1222;
  double phir_ttt = yr.f_222;
  double phir_tttt = yr.f_2222;

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

  out->f = sqrt(w2);
  out->f_1 = 0.5*pow(w2, -0.5) * w2_d;
  out->f_11 = 0.5*pow(w2, -0.5) * w2_dd - 0.25*pow(w2, -1.5) * w2_d*w2_d;
  out->f_2 = 0.5*pow(w2, -0.5) * w2_t;
  out->f_12 = 0.5*pow(w2, -0.5) * w2_dt - 0.25*pow(w2, -1.5) * w2_d*w2_t;
  out->f_22 = 0.5*pow(w2, -0.5) * w2_tt - 0.25*pow(w2, -1.5) * w2_t*w2_t;
}

void phi_ideal2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_ideal(comp, delta, tau);
  out->f = y.f;
  out->f_1 = y.f_1;
  out->f_2 = y.f_2;
  out->f_11 = y.f_11;
  out->f_12 = y.f_12;
  out->f_22 = y.f_22;
}

void phi_ideal_d2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_ideal(comp, delta, tau);
  out->f = y.f_1;
  out->f_1 = y.f_11;
  out->f_2 = y.f_12;
  out->f_11 = y.f_111;
  out->f_12 = y.f_112;
  out->f_22 = y.f_122;
}

void phi_ideal_t2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_ideal(comp, delta, tau);
  out->f = y.f_2;
  out->f_1 = y.f_12;
  out->f_2 = y.f_22;
  out->f_11 = y.f_112;
  out->f_12 = y.f_122;
  out->f_22 = y.f_222;
}

void phi_ideal_dd2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct y = phi_ideal4(comp, delta, tau);
  out->f = y.f_11;
  out->f_1 = y.f_111;
  out->f_2 = y.f_112;
  out->f_11 = y.f_1111;
  out->f_12 = y.f_1112;
  out->f_22 = y.f_1122;
}

void phi_ideal_dt2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct y = phi_ideal4(comp, delta, tau);
  out->f = y.f_12;
  out->f_1 = y.f_112;
  out->f_2 = y.f_122;
  out->f_11 = y.f_1112;
  out->f_12 = y.f_1122;
  out->f_22 = y.f_1222;
}

void phi_ideal_tt2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct y = phi_ideal4(comp, delta, tau);
  out->f = y.f_22;
  out->f_1 = y.f_122;
  out->f_2 = y.f_222;
  out->f_11 = y.f_1122;
  out->f_12 = y.f_1222;
  out->f_22 = y.f_2222;
}

void phi_resi2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_resi(comp, delta, tau);
  out->f = y.f;
  out->f_1 = y.f_1;
  out->f_2 = y.f_2;
  out->f_11 = y.f_11;
  out->f_12 = y.f_12;
  out->f_22 = y.f_22;
}

void phi_resi_d2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_resi(comp, delta, tau);
  out->f = y.f_1;
  out->f_1 = y.f_11;
  out->f_2 = y.f_12;
  out->f_11 = y.f_111;
  out->f_12 = y.f_112;
  out->f_22 = y.f_122;
}

void phi_resi_t2(uint comp, double delta, double tau, f22_struct *out){
  f23_struct y = phi_resi(comp, delta, tau);
  out->f = y.f_2;
  out->f_1 = y.f_12;
  out->f_2 = y.f_22;
  out->f_11 = y.f_112;
  out->f_12 = y.f_122;
  out->f_22 = y.f_222;
}

void phi_resi_dd2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct y = phi_resi4(comp, delta, tau);
  out->f = y.f_11;
  out->f_1 = y.f_111;
  out->f_2 = y.f_112;
  out->f_11 = y.f_1111;
  out->f_12 = y.f_1112;
  out->f_22 = y.f_1122;
}

void phi_resi_dt2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct y = phi_resi4(comp, delta, tau);
  out->f = y.f_12;
  out->f_1 = y.f_112;
  out->f_2 = y.f_122;
  out->f_11 = y.f_1112;
  out->f_12 = y.f_1122;
  out->f_22 = y.f_1222;
}

void phi_resi_tt2(uint comp, double delta, double tau, f22_struct *out){
  f24_struct y = phi_resi4(comp, delta, tau);
  out->f = y.f_22;
  out->f_1 = y.f_122;
  out->f_2 = y.f_222;
  out->f_11 = y.f_1122;
  out->f_12 = y.f_1222;
  out->f_22 = y.f_2222;
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
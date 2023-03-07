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

#include"test_util.h"
#include"read_data.h"
#include"phi.h"
#include"props.h"
#include"props_hp.h"
#include"props_sp.h"
#include"props_up.h"
#include"props_tp.h"
#include"sat.h"
#include"delta.h"
#include"state.h"
#include<iostream>
#include<math.h>
#include<chrono>

//
//
//
// Check if two floats are relatively almost the same
//
//
inline bool rel_same(double x1, double x2, double tol){
  if(fabs(x1) < tol) return fabs(x1 - x2) < tol;
  return fabs((x1 - x2)/x1) < tol;
}

//
//
//
// Check derivatives of unary functions with finite difference
//
//
int fd1(
  test_fptr1 func,
  uint comp,
  double x,
  double h,
  double tv,
  double tol,
  bool dbg){
  
  if (std::isnan(x)) return 0;

  f12_struct y0, y1, y1b, y;

  y0 = func(comp, x);
  y1 = func(comp, x + h);
  y1b = func(comp, x - h);

  y.f = y0.f;
  y.f_1 = (y1.f - y1b.f)/2.0/h;
  y.f_11 = (y1.f_1 - y1b.f_1)/2.0/h;

  if(dbg){
    std::cout << "test_value = " << tv << std::endl;
    std::cout << "f1 = " << y0.f << std::endl;
    std::cout << "f1_1 = " << y0.f_1 << " f.d. approx = " << y.f_1 << std::endl;
    std::cout << "f1_11 = " << y0.f_11 << " f.d. approx = " << y.f_11 << std::endl;
  }
  if(!rel_same(y.f, tv, tol*10.0)) return 1;
  if(!rel_same(y0.f_1, y.f_1, tol*10.0)) return 2;
  if(!rel_same(y0.f_11, y.f_11, tol*10.0)) return 3;
  return 0;
}

//
//
//
// Check derivatives of binary functions with finite difference
//
//
int fd2(
  test_fptr2 func,
  uint comp,
  double x1,
  double x2,
  double h1,
  double h2,
  double tv,
  double tol,
  bool dbg){

  if (std::isnan(x1) || std::isnan(x2)){
    return 0;
  }

  f22_struct y0 = func(comp, x1, x2);
  f22_struct y1 = func(comp, x1 + h1, x2);
  f22_struct y2 = func(comp, x1, x2 + h2);
  f22_struct y1b = func(comp, x1 - h1, x2);
  f22_struct y2b = func(comp, x1, x2 - h2);

  f22_struct y; //center finite difference approx
  f22_struct yf; //forward finite difference approx
  f22_struct yb; //backward finite difference approx

  y.f = y0.f;
  y.f_1 = (y1.f - y1b.f)/h1/2.0;
  y.f_2 = (y2.f - y2b.f)/h2/2.0;
  y.f_11 = (y1.f_1 - y1b.f_1)/h1/2.0;
  y.f_12 = (y2.f_1 - y2b.f_1)/h2/2.0;
  y.f_22 = (y2.f_2 - y2b.f_2)/h2/2.0;

  yf.f = y0.f;
  yf.f_1 = (y1.f - y0.f)/h1;
  yf.f_2 = (y2.f - y0.f)/h2;
  yf.f_11 = (y1.f_1 - y0.f_1)/h1;
  yf.f_12 = (y2.f_1 - y0.f_1)/h2;
  yf.f_22 = (y2.f_2 - y0.f_2)/h2;

  yb.f = y0.f;
  yb.f_1 = (y0.f - y1b.f)/h1;
  yb.f_2 = (y0.f - y2b.f)/h2;
  yb.f_11 = (y0.f_1 - y1b.f_1)/h1;
  yb.f_12 = (y0.f_1 - y2b.f_1)/h2;
  yb.f_22 = (y0.f_2 - y2b.f_2)/h2;

  if(dbg){
    std::cout << std::endl;
    std::cout << "test value = " << tv << std::endl;
    std::cout << "f = " << y0.f << std::endl;
    std::cout << "f_1 = " << y0.f_1 << " fd: " << y.f_1 << " " << yb.f_1 << " " << yf.f_1 << std::endl;
    std::cout << "f_2 = " << y0.f_2 << " fd: " << y.f_2 << " " << yb.f_2 << " " << yf.f_2 << std::endl;
    std::cout << "f_11 = " << y0.f_11 << " fd: " << y.f_11 << " " << yb.f_11 << " " << yf.f_11 << std::endl;
    std::cout << "f_12 = " << y0.f_12 << " fd: " << y.f_12 << " " << yb.f_12 << " " << yf.f_12 << std::endl;
    std::cout << "f_22 = " << y0.f_22 << " fd: " << y.f_22 << " " << yb.f_22 << " " << yf.f_22 << std::endl;
  }

  if (!std::isnan(tv)){
    if (!rel_same(y0.f, tv, tol)) return 1;
  }
  if (rel_same(y.f_1, yf.f_1, 1e-3) && rel_same(y.f_1, yb.f_1, 1e-3)){ // trust fd approx?
    if (!rel_same(y0.f_1, y.f_1, tol)) return 2;
  }
  if (rel_same(y.f_2, yf.f_2, 1e-3) && rel_same(y.f_2, yb.f_2, 1e-3)){ // trust fd approx?
    if (!rel_same(y0.f_2, y.f_2, tol)) return 3;
  }
  if (rel_same(y.f_11, yf.f_11, 1e-3) && rel_same(y.f_11, yb.f_11, 1e-3)){ // trust fd approx?
    if (!rel_same(y0.f_11, y.f_11, tol)) return 4;
  }
  if (rel_same(y.f_12, yf.f_12, 1e-3) && rel_same(y.f_12, yb.f_12, 1e-3)){ // trust fd approx?
    if (!rel_same(y0.f_12, y.f_12, tol)) return 5;
  }
  if (rel_same(y.f_22, yf.f_22, 1e-3) && rel_same(y.f_22, yb.f_22, 1e-3)){ // trust fd approx?
    if (!rel_same(y0.f_22, y.f_22, tol)) return 6;
  }
  return 0;
}

//
//
// Check basic property function values and derivatives
//
//
uint test_basic_properties(uint comp, std::string comp_str, test_data::data_set_enum data_set, double u_off, double h_off, double s_off){
  int err = 0;
  unsigned long i;
  double tau, delta;
  parameters_struct *pdat = &cdata[comp];
  std::vector< std::vector<double> > dat = read_data(comp_str, data_set, u_off, h_off, s_off);
    
  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration =  stop - start;

  TEST_FUNCTION_OF_DELTA_TAU("P", memo2_pressure, dat[i][test_data::P_col]*1000, 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("s", memo2_entropy, dat[i][test_data::s_col], 1e-1, 0)
  TEST_FUNCTION_OF_DELTA_TAU("h", memo2_enthalpy, dat[i][test_data::h_col], 1e-1, 0)
  TEST_FUNCTION_OF_DELTA_TAU("u", memo2_internal_energy, dat[i][test_data::u_col], 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("cv", memo2_isochoric_heat_capacity, dat[i][test_data::cv_col], 1e-1, 0)
  TEST_FUNCTION_OF_DELTA_TAU("cp", memo2_isobaric_heat_capacity, dat[i][test_data::cp_col], 1e-1, 1)
  TEST_FUNCTION_OF_DELTA_TAU("w", memo2_speed_of_sound, dat[i][test_data::w_col], 1e-1, 0)
  TEST_FUNCTION_OF_DELTA_TAU("g", memo2_gibbs, dat[i][test_data::h_col] - h_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), 1e-1, 1)
  TEST_FUNCTION_OF_DELTA_TAU("f", memo2_helmholtz, dat[i][test_data::u_col] - u_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), 1e-1, 1)

  // If the properties are right phis are right, just check derivatives.
  TEST_FUNCTION_OF_DELTA_TAU("phii", memo2_phi_ideal, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phii_d", memo2_phi_ideal_d, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phii_t", memo2_phi_ideal_t, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phii_dd", memo2_phi_ideal_dd, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phii_dt", memo2_phi_ideal_dt, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phii_tt", memo2_phi_ideal_tt, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phir", memo2_phi_resi, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phir_d", memo2_phi_resi_d, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phir_t", memo2_phi_resi_t, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phir_dd", memo2_phi_resi_dd, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phir_dt", memo2_phi_resi_dt, nan("no check"), 1e-2, 0)
  TEST_FUNCTION_OF_DELTA_TAU("phir_tt", memo2_phi_resi_tt, nan("no check"), 1e-2, 0)

  // Transport properties
  if (pdat->expr_map[expr_idx::viscosity_idx]!=1000){
      TEST_FUNCTION_OF_DELTA_TAU("viscosity", memo2_viscosity, dat[i][test_data::visc_col], 1e-2, 0)
  }
  // For thermal conductivity, for a lot of reasons, (e.g. different correlation, no critical 
  // correction, ...) these may or may not be super accurate, just check they are roughly correct.
  // ignore errors but print more than 33% error for manual inspection. 
  if (pdat->expr_map[expr_idx::viscosity_idx]!=1000){
      TEST_FUNCTION_OF_DELTA_TAU("thermal conductivity", memo2_thermal_conductivity, dat[i][test_data::tc_col], 0.33, 1)
  }
  return 0;
}


//
//
// Check saturation curve values and derivatives
//
//
uint test_sat_curve(uint comp, std::string comp_str, double u_off, double h_off, double s_off){
  std::vector< std::vector<double> > sat_liq_data, sat_vap_data;
  sort_sat(comp_str, test_data::saturated_set, &sat_liq_data, &sat_vap_data);
  uint err = 0;
  unsigned long i;
  double tau, pressure, delta;
  parameters_struct *pdat = &cdata[comp];
  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration =  stop - start;

  TEST_FUNCTION_SAT("p_sat", "tau", sat_p, sat_liq_data[i][test_data::P_col]*1000, pdat->T_star/sat_liq_data[i][test_data::T_col], 1e-3, 1e-7)
  TEST_FUNCTION_SAT("tau_sat", "p", sat_tau, pdat->T_star/sat_liq_data[i][test_data::T_col], sat_liq_data[i][test_data::P_col]*1000, 1e-3, 1e-5)
  TEST_FUNCTION_SAT("p_sat", "T", sat_p_t, sat_liq_data[i][test_data::P_col]*1000, sat_liq_data[i][test_data::T_col], 1e-3, 1e-7)
  TEST_FUNCTION_SAT("T_sat", "p", sat_t, sat_liq_data[i][test_data::T_col], sat_liq_data[i][test_data::P_col]*1000, 1e-3, 1e-5)

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    rho_liq_sat(" << comp_str << ", tau) ";
  for(i=0; i<sat_liq_data.size(); ++i){
    tau = pdat->T_star/sat_liq_data[i][test_data::T_col];
    delta = sat_liq_data[i][test_data::rho_col]/pdat->rho_star;
    err = fd1(sat_delta_l, comp, tau, 1e-7, delta, 1, 0);
    if(err){
      std::cout << err;
      return err;
    }

  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  std::cout << "Passed " << 3*sat_liq_data.size() << " points in " << duration.count() << "s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    rho_vap_sat(" << comp_str << ", tau) ";
  for(i=0; i<sat_liq_data.size(); ++i){
    tau = pdat->T_star/sat_liq_data[i][test_data::T_col];
    delta = sat_vap_data[i][test_data::rho_col]/pdat->rho_star;
    err = fd1(sat_delta_v, comp, tau, 1e-8, delta, 1e-3, 0);
    if(err){
      std::cout << err;
      return err;
    }

  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  std::cout << "Passed " << 3*sat_liq_data.size() << " points in " << duration.count() << "s" << std::endl;

  return 0;
}

//
//
// Check delta(P, tau) functions and derivatives
//
//
uint test_delta_function(uint comp, std::string comp_str, test_data::data_set_enum data_set, double u_off, double h_off, double s_off){
  int err = 0;
  unsigned long i;
  double tau, delta, P, Psat;
  parameters_struct *pdat = &cdata[comp];
  std::vector< std::vector<double> > dat = read_data(comp_str, data_set);
  test_fptr2 delta_func;

  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration;

  if(data_set == test_data::vapor_set){
    TEST_FUNCTION_OF_P_TAU("delta_vapor", memo2_delta_vapor, dat[i][test_data::rho_col]/pdat->rho_star, 1e-2)
  }
  else{
    TEST_FUNCTION_OF_P_TAU("delta_liquid", memo2_delta_liquid, dat[i][test_data::rho_col]/pdat->rho_star, 1e-2)
  }
  return 0;
}


//
//
// Check state var change 
//
//
uint test_state(uint comp, std::string comp_str, test_data::data_set_enum data_set, double u_off, double h_off, double s_off){
  int err = 0;
  unsigned long i;
  double tau, delta, P, Psat, enth, entr, inte, tv;
  parameters_struct *pdat = &cdata[comp];
  std::vector< std::vector<double> > dat = read_data(comp_str, data_set);
  test_fptr2 h_func, s_func, u_func;

  auto start = std::chrono::high_resolution_clock::now();
  if(data_set == test_data::vapor_set){
    std::cout << "    h_vapor(" << comp_str << ", P, tau) ";
    h_func = memo2_enthalpy_vapor;
    s_func = memo2_entropy_vapor;
    u_func = memo2_internal_energy_vapor;
  }
  else{ // liquid or supercritical use the liquid call
    std::cout << "    h_liquid(" << comp_str << ", P, tau) ";
    h_func = memo2_enthalpy_liquid;
    s_func = memo2_entropy_liquid;
    u_func = memo2_internal_energy_liquid;
  }
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    enth = dat[i][test_data::h_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    if((P >= Psat && h_func == memo2_enthalpy_liquid) || (P <= Psat && h_func == memo2_enthalpy_vapor)){ 
        // make sure the phase is correct, this can be a little off due to lack of sig figs.
        if(fabs(P - pdat->Pc) < 0.1 && fabs(tau - tau_c(comp)) < 0.001) continue;
        err = fd2(h_func, comp, P, tau, 1e-3, 1e-8, enth - h_off, 1e-2, 0);
        if(err){
            std::cout << err;
            std::cout << " rho = " << delta*pdat->rho_star << ", P = " << P << " T = " << pdat->T_star/tau << std::endl;
            fd2(h_func, comp, P, tau, 1e-3, 1e-8, enth - h_off, 1e-2, 1);
            return err;
        }
    }
  }
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  if(data_set == test_data::vapor_set){
    std::cout << "    s_vapor(" << comp_str << ", P, tau) ";
  }
  else{ // liquid or supercritical use the liquid call
    std::cout << "    s_liquid(" << comp_str << ", P, tau) ";
  }
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    entr = dat[i][test_data::s_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    if((P >= Psat && h_func == memo2_enthalpy_liquid) || (P <= Psat && h_func == memo2_enthalpy_vapor)){ 
        // make sure the phase is correct, this can be a little off due to lack of sig figs.
        err = fd2(s_func, comp, P, tau, 1e-2, 1e-6, entr - s_off, 1e-1, 0);
        if(err){
            std::cout << err;
            std::cout << " rho = " << delta*pdat->rho_star << ", P = " << P << " T = " << pdat->T_star/tau << " tau = " << tau << std::endl;
            fd2(s_func, comp, P, tau, 1e-2, 1e-6, entr - s_off, 1e-1, 1);
            return err;
        }
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;
  
  start = std::chrono::high_resolution_clock::now();
  if(data_set == test_data::vapor_set){
    std::cout << "    u_vapor(" << comp_str << ", P, tau) ";
  }
  else{ // liquid or supercritical use the liquid call
    std::cout << "    u_liquid(" << comp_str << ", P, tau) ";
  }
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    inte = dat[i][test_data::u_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    if((P >= Psat && h_func == memo2_enthalpy_liquid) || (P <= Psat && h_func == memo2_enthalpy_vapor)){ 
        // make sure the phase is correct, this can be a little off due to lack of sig figs.
        if(fabs(P - pdat->Pc) < 0.1 && fabs(tau - tau_c(comp)) < 0.001) continue;
        err = fd2(u_func, comp, P, tau, 1e-3, 1e-8, inte - u_off, 1e-2, 0);
        if(err){
            std::cout << err;
            std::cout << " rho = " << delta*pdat->rho_star << ", P = " << P << " T = " << pdat->T_star/tau << std::endl;
            fd2(u_func, comp, P, tau, 1e-3, 1e-8, inte - u_off, 1e-2, 1);
            return err;
        }
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;  

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    tau(" << comp_str << ", h, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    enth = dat[i][test_data::h_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    err = fd2(memo2_tau_hp, comp, enth - h_off, P, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << err;
        std::cout << " rho = " << delta*pdat->rho_star << ", P = " << P << " T = " << pdat->T_star/tau << " tau = " << tau << std::endl;
        fd2(memo2_tau_hp, comp, enth - h_off, P, 1e-3, 1e-8, tau, 1e-2, 1);
        return err;
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;  

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    tau(" << comp_str << ", s, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    entr = dat[i][test_data::s_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    err = fd2(memo2_tau_sp, comp, entr  - s_off, P, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << err;
        std::cout << " rho = " << delta*pdat->rho_star << ", P = " << P << " T = " << pdat->T_star/tau << std::endl;
        fd2(memo2_tau_sp, comp, entr  - s_off, P, 1e-2, 1e-7, tau, 1e-2, 1);
        return err;
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    tau(" << comp_str << ", u, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    inte = dat[i][test_data::u_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    err = fd2(memo2_tau_up, comp, inte  - u_off, P, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << err;
        std::cout << " rho = " << delta*pdat->rho_star << ", P = " << P << " T = " << pdat->T_star/tau << std::endl;
        fd2(memo2_tau_up, comp, inte  - u_off, P, 1e-3, 1e-8, tau, 1e-2, 1);
        return err;
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    vf(" << comp_str << ", h, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    enth = dat[i][test_data::h_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    if(data_set == test_data::vapor_set){
        tv = 1.0;
    }
    else{
        tv = 0.0;
    }
    err = fd2(memo2_vf_hp, comp, enth - h_off, P, 1e-3, 1e-8, tv, 1e-2, 0);
    if(err){
        fd2(memo2_vf_hp, comp, enth - h_off, P, 1e-3, 1e-8, tv, 1e-2, 1);
        return err;
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    vf(" << comp_str << ", s, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    entr = dat[i][test_data::s_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    if(data_set == test_data::vapor_set){
        tv = 1.0;
    }
    else{
        tv = 0.0;
    }
    err = fd2(memo2_vf_sp, comp, entr - s_off, P, 1e-3, 1e-8, tv, 1e-2, 0);
    if(err){
        fd2(memo2_vf_sp, comp, entr - s_off, P, 1e-3, 1e-8, tv, 1e-2, 1);
        return err;
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::cout << "    vf(" << comp_str << ", u, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = pdat->T_star/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/pdat->rho_star;
    inte = dat[i][test_data::u_col];
    P = dat[i][test_data::P_col]*1000.0;
    Psat = sat_p(comp, tau).f;
    if(data_set == test_data::vapor_set){
        tv = 1.0;
    }
    else{
        tv = 0.0;
    }
    err = fd2(memo2_vf_up, comp, inte - u_off, P, 1e-2, 1e-8, tv, 1e-2, 0);
    if(err){
        fd2(memo2_vf_up, comp, inte - u_off, P, 1e-2, 1e-8, tv, 1e-2, 1);
        return err;
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration =  stop - start;
  // The dat size is multiplied by 5 since there are 4 extra points evaluated for finite difference tests.
  std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

  TEST_FUNCTION_OF_STATE_VARS("T", "h", memo2_temperature_hp, dat[i][test_data::T_col], dat[i][test_data::h_col] - h_off, 1e-2)
  TEST_FUNCTION_OF_STATE_VARS("u", "h", memo2_internal_energy_hp, dat[i][test_data::u_col] - u_off, dat[i][test_data::h_col] - h_off, 1e-2)
  TEST_FUNCTION_OF_STATE_VARS("s", "h", memo2_entropy_hp, dat[i][test_data::s_col] - s_off, dat[i][test_data::h_col] - h_off, 5e-2)
  TEST_FUNCTION_OF_STATE_VARS("f", "h", memo2_helmholtz_hp, dat[i][test_data::u_col] - u_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("g", "h", memo2_gibbs_hp, dat[i][test_data::h_col] - h_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("cv", "h", memo2_isochoric_heat_capacity_hp, dat[i][test_data::cv_col], dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("cp", "h", memo2_isobaric_heat_capacity_hp, dat[i][test_data::cp_col], dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("w", "h", memo2_speed_of_sound_hp, dat[i][test_data::w_col], dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("v", "h", memo2_specific_volume_hp, 1/dat[i][test_data::rho_col], dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("mu", "h", memo2_viscosity_hp, dat[i][test_data::visc_col], dat[i][test_data::h_col] - h_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("lambda", "h", memo2_thermal_conductivity_hp, dat[i][test_data::tc_col], dat[i][test_data::h_col] - h_off, 0.33)

  TEST_FUNCTION_OF_STATE_VARS("T", "s", memo2_temperature_sp, dat[i][test_data::T_col], dat[i][test_data::s_col] - s_off, 1e-2)
  TEST_FUNCTION_OF_STATE_VARS("h", "s", memo2_enthalpy_sp, dat[i][test_data::h_col] - h_off, dat[i][test_data::s_col] - s_off, 1e-2)
  TEST_FUNCTION_OF_STATE_VARS("u", "s", memo2_internal_energy_sp, dat[i][test_data::u_col] - u_off, dat[i][test_data::s_col] - s_off, 5e-2)
  TEST_FUNCTION_OF_STATE_VARS("f", "s", memo2_helmholtz_sp, dat[i][test_data::u_col] - u_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("g", "s", memo2_gibbs_sp, dat[i][test_data::h_col] - h_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("cv", "s", memo2_isochoric_heat_capacity_sp, dat[i][test_data::cv_col], dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("cp", "s", memo2_isobaric_heat_capacity_sp, dat[i][test_data::cp_col], dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("w", "s", memo2_speed_of_sound_sp, dat[i][test_data::w_col], dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("v", "s", memo2_specific_volume_sp, 1/dat[i][test_data::rho_col], dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("mu", "s", memo2_viscosity_sp, dat[i][test_data::visc_col], dat[i][test_data::s_col] - s_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("lambda", "s", memo2_thermal_conductivity_sp, dat[i][test_data::tc_col], dat[i][test_data::s_col] - s_off, 0.33)

  TEST_FUNCTION_OF_STATE_VARS("T", "u", memo2_temperature_up, dat[i][test_data::T_col], dat[i][test_data::u_col] - u_off, 1e-2)
  TEST_FUNCTION_OF_STATE_VARS("h", "u", memo2_enthalpy_up, dat[i][test_data::h_col] - h_off, dat[i][test_data::u_col] - u_off, 1e-2)
  TEST_FUNCTION_OF_STATE_VARS("s", "u", memo2_entropy_up, dat[i][test_data::s_col] - s_off, dat[i][test_data::u_col] - u_off, 5e-2)
  TEST_FUNCTION_OF_STATE_VARS("f", "u", memo2_helmholtz_up, dat[i][test_data::u_col] - u_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("g", "u", memo2_gibbs_up, dat[i][test_data::h_col] - h_off - dat[i][test_data::T_col]*(dat[i][test_data::s_col] - s_off), dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("cv", "u", memo2_isochoric_heat_capacity_up, dat[i][test_data::cv_col], dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("cp", "u", memo2_isobaric_heat_capacity_up, dat[i][test_data::cp_col], dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("w", "u", memo2_speed_of_sound_up, dat[i][test_data::w_col], dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("v", "u", memo2_specific_volume_up, 1/dat[i][test_data::rho_col], dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("mu", "u", memo2_viscosity_up, dat[i][test_data::visc_col], dat[i][test_data::u_col] - u_off, 1e-1)
  TEST_FUNCTION_OF_STATE_VARS("lambda", "u", memo2_thermal_conductivity_up, dat[i][test_data::tc_col], dat[i][test_data::u_col] - u_off, 0.33)

  if(data_set == test_data::vapor_set){
    TEST_FUNCTION_OF_TP("h", "T", memo2_enthalpy_vap_tp, dat[i][test_data::h_col] - h_off, dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("s", "T", memo2_entropy_vap_tp, dat[i][test_data::s_col] - s_off, dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("u", "T", memo2_internal_energy_vap_tp, dat[i][test_data::u_col] - u_off, dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("cv", "T", memo2_isochoric_heat_capacity_vap_tp, dat[i][test_data::cv_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("cp", "T", memo2_isobaric_heat_capacity_vap_tp, dat[i][test_data::cp_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("w", "T", memo2_speed_of_sound_vap_tp, dat[i][test_data::w_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("v", "T", memo2_specific_volume_vap_tp, 1.0/dat[i][test_data::rho_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("viscosity", "T", memo2_viscosity_vap_tp, dat[i][test_data::visc_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity", "T", memo2_thermal_conductivity_vap_tp, dat[i][test_data::tc_col], dat[i][test_data::T_col], 0.33)

    TEST_FUNCTION_OF_TP("s_vap", "h", memo2_entropy_vap_hp, dat[i][test_data::s_col] - s_off, dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("u_vap", "h", memo2_internal_energy_vap_hp, dat[i][test_data::u_col] - u_off, dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("cv_vap", "h", memo2_isochoric_heat_capacity_vap_hp, dat[i][test_data::cv_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("cp_vap", "h", memo2_isobaric_heat_capacity_vap_hp, dat[i][test_data::cp_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("w_vap", "h", memo2_speed_of_sound_vap_hp, dat[i][test_data::w_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("v_vap", "h", memo2_specific_volume_vap_hp, 1.0/dat[i][test_data::rho_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("viscosity_vap", "h", memo2_viscosity_vap_hp, dat[i][test_data::visc_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity_vap", "h", memo2_thermal_conductivity_vap_hp, dat[i][test_data::tc_col], dat[i][test_data::h_col] - h_off, 0.33)
  
    TEST_FUNCTION_OF_TP("h_vap", "s", memo2_enthalpy_vap_sp, dat[i][test_data::h_col] - h_off, dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("u_vap", "s", memo2_internal_energy_vap_sp, dat[i][test_data::u_col] - u_off, dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("cv_vap", "s", memo2_isochoric_heat_capacity_vap_sp, dat[i][test_data::cv_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("cp_vap", "s", memo2_isobaric_heat_capacity_vap_sp, dat[i][test_data::cp_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("w_vap", "s", memo2_speed_of_sound_vap_sp, dat[i][test_data::w_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("v_vap", "s", memo2_specific_volume_vap_sp, 1.0/dat[i][test_data::rho_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("viscosity_vap", "s", memo2_viscosity_vap_sp, dat[i][test_data::visc_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity_vap", "s", memo2_thermal_conductivity_vap_sp, dat[i][test_data::tc_col], dat[i][test_data::s_col] - s_off, 0.33)

    TEST_FUNCTION_OF_TP("h_vap", "u", memo2_enthalpy_vap_up, dat[i][test_data::h_col] - h_off, dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("s_vap", "u", memo2_entropy_vap_up, dat[i][test_data::s_col] - s_off, dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("cv_vap", "u", memo2_isochoric_heat_capacity_vap_up, dat[i][test_data::cv_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("cp_vap", "u", memo2_isobaric_heat_capacity_vap_up, dat[i][test_data::cp_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("w_vap", "u", memo2_speed_of_sound_vap_up, dat[i][test_data::w_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("v_vap", "u", memo2_specific_volume_vap_up, 1.0/dat[i][test_data::rho_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("viscosity_vap", "u", memo2_viscosity_vap_up, dat[i][test_data::visc_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity_vap", "u", memo2_thermal_conductivity_vap_up, dat[i][test_data::tc_col], dat[i][test_data::u_col] - u_off, 0.33)

  }
  else{
    TEST_FUNCTION_OF_TP("h", "T", memo2_enthalpy_liq_tp, dat[i][test_data::h_col] - h_off, dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("s", "T", memo2_entropy_liq_tp, dat[i][test_data::s_col] - s_off, dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("u", "T", memo2_internal_energy_liq_tp, dat[i][test_data::u_col] - u_off, dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("cv", "T", memo2_isochoric_heat_capacity_liq_tp, dat[i][test_data::cv_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("cp", "T", memo2_isobaric_heat_capacity_liq_tp, dat[i][test_data::cp_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("w", "T", memo2_speed_of_sound_liq_tp, dat[i][test_data::w_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("v", "T", memo2_specific_volume_liq_tp, 1.0/dat[i][test_data::rho_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("viscosity", "T", memo2_viscosity_liq_tp, dat[i][test_data::visc_col], dat[i][test_data::T_col], 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity", "T", memo2_thermal_conductivity_liq_tp, dat[i][test_data::tc_col], dat[i][test_data::T_col], 0.33)

    TEST_FUNCTION_OF_TP("s_liq", "h", memo2_entropy_liq_hp, dat[i][test_data::s_col] - s_off, dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("u_liq", "h", memo2_internal_energy_liq_hp, dat[i][test_data::u_col] - u_off, dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("cv_liq", "h", memo2_isochoric_heat_capacity_liq_hp, dat[i][test_data::cv_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("cp_liq", "h", memo2_isobaric_heat_capacity_liq_hp, dat[i][test_data::cp_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("w_liq", "h", memo2_speed_of_sound_liq_hp, dat[i][test_data::w_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("v_liq", "h", memo2_specific_volume_liq_hp, 1.0/dat[i][test_data::rho_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("viscosity_liq", "h", memo2_viscosity_liq_hp, dat[i][test_data::visc_col], dat[i][test_data::h_col] - h_off, 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity_liq", "h", memo2_thermal_conductivity_liq_hp, dat[i][test_data::tc_col], dat[i][test_data::h_col] - h_off, 0.33)

    TEST_FUNCTION_OF_TP("h_vap", "s", memo2_enthalpy_liq_sp, dat[i][test_data::h_col] - h_off, dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("u_liq", "s", memo2_internal_energy_liq_sp, dat[i][test_data::u_col] - u_off, dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("cv_liq", "s", memo2_isochoric_heat_capacity_liq_sp, dat[i][test_data::cv_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("cp_liq", "s", memo2_isobaric_heat_capacity_liq_sp, dat[i][test_data::cp_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("w_liq", "s", memo2_speed_of_sound_liq_sp, dat[i][test_data::w_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("v_liq", "s", memo2_specific_volume_liq_sp, 1.0/dat[i][test_data::rho_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("viscosity_liq", "s", memo2_viscosity_liq_sp, dat[i][test_data::visc_col], dat[i][test_data::s_col] - s_off, 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity_liq", "s", memo2_thermal_conductivity_liq_sp, dat[i][test_data::tc_col], dat[i][test_data::s_col] - s_off, 0.33)

    TEST_FUNCTION_OF_TP("h_liq", "u", memo2_enthalpy_liq_up, dat[i][test_data::h_col] - h_off, dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("s_liq", "u", memo2_entropy_liq_up, dat[i][test_data::s_col] - s_off, dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("cv_liq", "u", memo2_isochoric_heat_capacity_liq_up, dat[i][test_data::cv_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("cp_liq", "u", memo2_isobaric_heat_capacity_liq_up, dat[i][test_data::cp_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("w_liq", "u", memo2_speed_of_sound_liq_up, dat[i][test_data::w_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("v_liq", "u", memo2_specific_volume_liq_up, 1.0/dat[i][test_data::rho_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("viscosity_liq", "u", memo2_viscosity_liq_up, dat[i][test_data::visc_col], dat[i][test_data::u_col] - u_off, 1e-1)
    TEST_FUNCTION_OF_TP("thermal_conductivity_liq", "u", memo2_thermal_conductivity_liq_up, dat[i][test_data::tc_col], dat[i][test_data::u_col] - u_off, 0.33)
  }
  

  return 0;
}

//
//
// Check vapor fractions and other functions along sat curve
//
//
uint test_sat_curve_more(uint comp, std::string comp_str, double u_off, double h_off, double s_off){
  std::vector< std::vector<double> > sat_liq_data, sat_vap_data;
  sort_sat(comp_str, test_data::saturated_set, &sat_liq_data, &sat_vap_data);
  uint err = 0;
  unsigned long i;
  double tau, pressure, delta_l, delta_v;
  parameters_struct *pdat = &cdata[comp];
  f22_struct hl_vec, hv_vec, sl_vec, sv_vec, ul_vec, uv_vec;
  double ht, st, ut;

  std::cout << "Testing saturated and two-phase regions ... ";

  for(i=0; i<sat_liq_data.size(); ++i){
    tau = pdat->T_star/sat_liq_data[i][test_data::T_col];
    pressure = sat_liq_data[i][test_data::P_col]*1000;
    delta_l = sat_liq_data[i][test_data::rho_col]/pdat->rho_star;
    delta_v = sat_vap_data[i][test_data::rho_col]/pdat->rho_star;
    enthalpy2(comp, delta_l, tau, &hl_vec);
    enthalpy2(comp, delta_v, tau, &hv_vec);
    entropy2(comp, delta_l, tau, &sl_vec);
    entropy2(comp, delta_v, tau, &sv_vec);
    internal_energy2(comp, delta_l, tau, &ul_vec);
    internal_energy2(comp, delta_v, tau, &uv_vec);

    // Test H, P
    ht = 0.5*hl_vec.f + 0.5*hv_vec.f;
    err = fd2(memo2_tau_hp, comp, ht, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    err = fd2(memo2_vf_hp, comp, ht, pressure, 1e-3, 1e-8, 0.5, 1e-2, 0);
    if(err){
        fd2(memo2_vf_hp, comp, ht, pressure, 1e-3, 1e-8, 0.5, 1e-2, 1);
        std::cout << "fail vfhp " << pressure << ", " << delta_l << ", " << delta_v << std::endl;
        return 1;
    }
    ht = 0.25*hl_vec.f + 0.75*hv_vec.f;
    err = fd2(memo2_tau_hp, comp, ht, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
     if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    err = fd2(memo2_vf_hp, comp, ht, pressure, 1e-3, 1e-8, 0.75, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    ht = hl_vec.f;
    err = fd2(memo2_tau_hp, comp, ht, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
     if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    ht = hv_vec.f;
    err = fd2(memo2_tau_hp, comp, ht, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
     if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }

    // Test S, P
    st = 0.5*sl_vec.f + 0.5*sv_vec.f;
    err = fd2(memo2_tau_sp, comp, st, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    st = 0.25*sl_vec.f + 0.75*sv_vec.f;
    err = fd2(memo2_tau_sp, comp, st, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
     if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    err = fd2(memo2_vf_sp, comp, st, pressure, 1e-3, 1e-8, 0.75, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    st = sl_vec.f;
    err = fd2(memo2_tau_sp, comp, st, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    st = sv_vec.f;
    err = fd2(memo2_tau_sp, comp, st, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }

    // Test U, P
    ut = 0.5*ul_vec.f + 0.5*uv_vec.f;
    err = fd2(memo2_tau_up, comp, ut, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    ut = 0.25*ul_vec.f + 0.75*uv_vec.f;
    err = fd2(memo2_tau_up, comp, ut, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
     if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    err = fd2(memo2_vf_up, comp, ut, pressure, 1e-3, 1e-8, 0.75, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    ut = ul_vec.f;
    err = fd2(memo2_tau_up, comp, ut, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
    ut = uv_vec.f;
    err = fd2(memo2_tau_up, comp, ut, pressure, 1e-3, 1e-8, tau, 1e-2, 0);
    if(err){
        std::cout << "fail" << std::endl;
        return 1;
    }
  }
  std::cout << "passed" << std::endl;
  return 0;
}

uint run_set_mixed(uint comp, std::string comp_str, double u_off, double h_off, double s_off){
    uint err = 0;

    std::cout << std::endl;
    std::cout << "Test basic " << comp_str << " mixed properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::mixed_set, u_off, h_off, s_off);
    if(err){
        return err;
    }
    return 0;
}

uint run_set_all(uint comp, std::string comp_str, double u_off, double h_off, double s_off){
    uint err = 0;

    std::cout << std::endl;
    std::cout << "Test basic " << comp_str << " vapor properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::vapor_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test basic " << comp_str << " supercritical properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::supercritical_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test basic " << comp_str << " liquid properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::liquid_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " sat. curve" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_sat_curve(comp, comp_str, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " liquid delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(comp, comp_str, test_data::liquid_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " vapor delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(comp, comp_str, test_data::vapor_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " supercritical delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(comp, comp_str, test_data::supercritical_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " functions for state var change on liquid data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(comp, comp_str, test_data::liquid_set, u_off, h_off, s_off);
    if(err){
       return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " functions for state var change on vapor data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(comp, comp_str, test_data::vapor_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " functions for state var change on sc data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(comp, comp_str, test_data::supercritical_set, u_off, h_off, s_off);
    if(err){
        return err;
    }

    std::cout << std::endl;
    std::cout << "Test " << comp_str << " two phase" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_sat_curve_more(comp, comp_str, u_off, h_off, s_off);
    if(err){
        return err;
    }

  return 0;
}









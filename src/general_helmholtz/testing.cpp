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

#include"test_data/read_data.h"
#include"phi.h"
#include"props.h"
#include"sat.h"
#include"solver.h"
#include"delta.h"
#include"state.h"
#include"components/h2o.h"
#include"testing.h"
#include <iostream>
#include <math.h>

inline bool rel_same(double x1, double x2, double tol){
  if(fabs(x1) < tol) return fabs(x1 - x2) < tol;
  return fabs((x1 - x2)/x1) < tol;
}

int fd1(
  test_fptr1 func,
  comp_enum comp,
  double x,
  std::vector<double> *yvec_ptr,
  double h,
  double tv,
  double tol,
  bool dbg){

  std::vector<double> *yvec_ptr0, *yvec_ptr1, *yvec_ptr1b;

  yvec_ptr0 = func(comp, x);
  yvec_ptr1 = func(comp, x + h);
  yvec_ptr1b = func(comp, x - h);

  yvec_ptr->resize(3);
  yvec_ptr->at(0) = yvec_ptr0->at(0);
  yvec_ptr->at(1) = (yvec_ptr1->at(0) - yvec_ptr1b->at(0))/2.0/h;
  yvec_ptr->at(2) = (yvec_ptr1->at(1) - yvec_ptr1b->at(1))/2.0/h;

  if(dbg){
    std::cout << "test_value = " << tv << std::endl;
    std::cout << "f1 = " << yvec_ptr0->at(0) << std::endl;
    std::cout << "f1_1 = " << yvec_ptr0->at(1) << " f.d. approx = " << yvec_ptr->at(1) << std::endl;
    std::cout << "f1_11 = " << yvec_ptr0->at(2) << " f.d. approx = " << yvec_ptr->at(2) << std::endl;
  }

  if(!rel_same(yvec_ptr->at(0), tv, tol)) return 1;
  if(!rel_same(yvec_ptr->at(1), yvec_ptr0->at(1), tol)) return 2;
  if(!rel_same(yvec_ptr->at(2), yvec_ptr0->at(2), tol)) return 3;

  return 0;
}

int fd2(
  test_fptr2 func,
  comp_enum comp,
  double x1,
  double x2,
  std::vector<double> *yvec_ptr,
  double h1,
  double h2,
  double tv,
  double tol,
  bool dbg){

  std::vector<double> *yvec_ptr0, *yvec_ptr1, *yvec_ptr2, *yvec_ptr1b, *yvec_ptr2b;
  std::vector<double> yvecf(6), yvecb(6);

  if (std::isnan(x1) || std::isnan(x2)){
    return 0;
  }

  yvec_ptr0 = func(comp, x1, x2);
  yvec_ptr1 = func(comp, x1 + h1, x2);
  yvec_ptr2 = func(comp, x1, x2 + h2);
  yvec_ptr1b = func(comp, x1 - h1, x2);
  yvec_ptr2b = func(comp, x1, x2 - h2);

  yvec_ptr->resize(6);
  yvec_ptr->at(f2) = yvec_ptr0->at(f2);
  yvec_ptr->at(f2_1) = (yvec_ptr1->at(f2) - yvec_ptr1b->at(f2))/h1/2.0;
  yvec_ptr->at(f2_2) = (yvec_ptr2->at(f2) - yvec_ptr2b->at(f2))/h2/2.0;
  yvec_ptr->at(f2_11) = (yvec_ptr1->at(f2_1) - yvec_ptr1b->at(f2_1))/h1/2.0;
  yvec_ptr->at(f2_12) = (yvec_ptr2->at(f2_1) - yvec_ptr2b->at(f2_1))/h2/2.0;
  yvec_ptr->at(f2_22) = (yvec_ptr2->at(f2_2) - yvec_ptr2b->at(f2_2))/h2/2.0;

  yvecf[f2_1] = (yvec_ptr1->at(f2) - yvec_ptr0->at(f2))/h1;
  yvecf[f2_2] = (yvec_ptr2->at(f2) - yvec_ptr0->at(f2))/h2;
  yvecf[f2_11] = (yvec_ptr1->at(f2_1) - yvec_ptr0->at(f2_1))/h1;
  yvecf[f2_12] = (yvec_ptr2->at(f2_1) - yvec_ptr0->at(f2_1))/h2;
  yvecf[f2_22] = (yvec_ptr2->at(f2_2) - yvec_ptr0->at(f2_2))/h2;

  yvecb[f2_1] = (yvec_ptr0->at(f2) - yvec_ptr1b->at(f2))/h1;
  yvecb[f2_2] = (yvec_ptr0->at(f2) - yvec_ptr2b->at(f2))/h2;
  yvecb[f2_11] = (yvec_ptr0->at(f2_1) - yvec_ptr1b->at(f2_1))/h1;
  yvecb[f2_12] = (yvec_ptr0->at(f2_1) - yvec_ptr2b->at(f2_1))/h2;
  yvecb[f2_22] = (yvec_ptr0->at(f2_2) - yvec_ptr2b->at(f2_2))/h2;

  if(dbg){
    std::cout << "test value = " << tv << std::endl;
    std::cout << "f2 = " << yvec_ptr0->at(0) << std::endl;
    std::cout << "f2_1 = " << yvec_ptr0->at(f2_1) << " f.d. approx = " << yvec_ptr->at(f2_1) << " " << yvecb[f2_1]<< " " << yvecf[f2_1]<< std::endl;
    std::cout << "f2_2 = " << yvec_ptr0->at(f2_2) << " f.d. approx = " << yvec_ptr->at(f2_2) << " " << yvecb[f2_2]<< " " << yvecf[f2_2]<< std::endl;
    std::cout << "f2_11 = " << yvec_ptr0->at(f2_11) << " f.d. approx = " << yvec_ptr->at(f2_11) << " " << yvecb[f2_11]<< " " << yvecf[f2_11]<< std::endl;
    std::cout << "f2_12 = " << yvec_ptr0->at(f2_12) << " f.d. approx = " << yvec_ptr->at(f2_12) << " " << yvecb[f2_12]<< " " << yvecf[f2_12]<< std::endl;
    std::cout << "f2_22 = " << yvec_ptr0->at(f2_22) << " f.d. approx = " << yvec_ptr->at(f2_22) << " " << yvecb[f2_22]<< " " << yvecf[f2_22]<< std::endl;
  }

  if (!std::isnan(tv)){
    if (!rel_same(yvec_ptr0->at(0), tv, tol)) return 1;
  }
  if (rel_same(yvecb[f2_1], yvecf[f2_1], 1e-2)){ // rough fd accuracy
    if (!rel_same(yvec_ptr0->at(f2_1), yvec_ptr->at(f2_1), tol)) return 2;
  }
  if (rel_same(yvecb[f2_2], yvecf[f2_2], 1e-2)){ // rough fd accuracy
    if (!rel_same(yvec_ptr0->at(f2_2), yvec_ptr->at(f2_2), tol)) return 3;
  }
  if (rel_same(yvecb[f2_11], yvecf[f2_11], 1e-2)){ // rough fd accuracy
    if (!rel_same(yvec_ptr0->at(f2_11), yvec_ptr->at(f2_11), tol)) return 4;
  }
  if (rel_same(yvecb[f2_12], yvecf[f2_12], 1e-2)){ // rough fd accuracy
    if (!rel_same(yvec_ptr0->at(f2_12), yvec_ptr->at(f2_12), tol)) return 5;
  }
  if (rel_same(yvecb[f2_22], yvecf[f2_22], 1e-2)){ // rough fd accuracy
    if (!rel_same(yvec_ptr0->at(f2_22), yvec_ptr->at(f2_22), tol)) return 6;
  }
  return 0;
}

double solver_test_function1(double x, void *dat){
  // simple test function with roots 3 and 10
  return (x - 3) * (x - 10);
}

void solver_test_function1g(double x, std::vector<double> *out, void *dat){
  // simple test function with roots 3 and 10 with 1st derivative
  out->resize(2);
  out->at(0) = (x - 3) * (x - 10);
  out->at(1) = (x - 3) + (x - 10);
}

void solver_test_function1gh(double x, std::vector<double> *out, void *dat){
  // simple test function with roots 3 and 10 with 1st and 2nd derivative
  out->resize(3);
  out->at(0) = (x - 3) * (x - 10);
  out->at(1) = (x - 3) + (x - 10);
  out->at(2) = 2;
}


int test_bracket1(bool dbg){
  int err=0, n=0;
  double sol;
  n = bracket(solver_test_function1, 4, 20, &sol, 40, 1e-7, 1e-7);
  if(!rel_same(10.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations: " << n << " solution: " << sol << std::endl;
  }
  n = bracket(solver_test_function1, -10, 5, &sol, 40, 1e-7, 1e-7);
  if(!rel_same(3.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations: " << n << " solution: " << sol << std::endl;
  }
  return err;
}

int test_halley1(bool dbg){
  int err=0, n=0;
  double sol;
  std::vector<double> fg;

  n = halley(solver_test_function1gh, 15.0, &sol, &fg, 20, 1e-7);
  if(!rel_same(10.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations: " << n << " solution: " << sol << std::endl;
  }
  n = halley(solver_test_function1gh, -3, &sol, &fg, 20, 1e-7);
  if(!rel_same(3.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations: " << n << " solution: " << sol << std::endl;
  }
  return err;
}

int test_newton_ls1(bool dbg){
  int err=0, n=0;
  double sol;
  std::vector<double> fg;

  n = newton_ls(solver_test_function1g, solver_test_function1, 15.0, &sol, &fg, 20, 1e-7);
  if(!rel_same(10.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations : " << n << " solution: " << sol << std::endl;
  }
  n = newton_ls(solver_test_function1g, solver_test_function1, -3, &sol, &fg, 20, 1e-7);
  if(!rel_same(3.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations : " << n << " solution: " << sol << std::endl;
  }
  n = newton_ls(solver_test_function1g, NULL, 15.0, &sol, &fg, 20, 1e-7, NULL, 0);
  if(!rel_same(10.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations : " << n << " solution: " << sol << std::endl;
  }
  n = newton_ls(solver_test_function1g, NULL, -3, &sol, &fg, 20, 1e-7, NULL, 0);
  if(!rel_same(3.0, sol, 1e-7)){
    ++err;
  }
  if(dbg){
    std::cout << "iterations : " << n << " solution: " << sol << std::endl;
  }

  return err;
}

uint test_basic_properties(comp_enum comp, test_data::data_set_enum data_set){
  std::vector<double> p_vec_fd;
  int err = 0;
  unsigned long i;
  double tau, delta;
  std::vector< std::vector<double> > dat = read_data(comp, data_set);
  std::string comp_str = comp_enum_table[comp];

  std::cout << "P(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_pressure, comp, delta, tau, &p_vec_fd, 1e-10, 1e-7, dat[i][test_data::P_col]*1000, 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "S(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_entropy, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::s_col], 1e-2, 0);
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "H(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_enthalpy, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::h_col], 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "U(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_internal_energy, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::u_col], 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "cv(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_isochoric_heat_capacity, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::cv_col], 1e-2, 0);
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "cp(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_isobaric_heat_capacity, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::cp_col], 1e-2, 0);
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "w(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_speed_of_sound, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::w_col], 1e-2, 0);
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "g(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_gibbs, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::h_col] - dat[i][test_data::T_col]*dat[i][test_data::s_col], 1e-1, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "f(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_helmholtz, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, dat[i][test_data::u_col] - dat[i][test_data::T_col]*dat[i][test_data::s_col], 1e-1, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phii(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_ideal, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phii_d(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_ideal_d, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phii_t(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_ideal_t, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phii_dd(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    // The tolerances may seem a little loose, but the data doesn't have quite enough sig figs.
    err = fd2(memo2_phi_ideal_dd, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phii_dt(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    // The tolerances may seem a little loose, but the data doesn't have quite enough sig figs.
    err = fd2(memo2_phi_ideal_dt, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phii_tt(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    // The tolerances may seem a little loose, but the data doesn't have quite enough sig figs.
    err = fd2(memo2_phi_ideal_tt, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phir(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_resi, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phir_d(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_resi_d, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phir_t(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_resi_t, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phir_dd(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_resi_dd, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phir_dt(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_resi_dt, comp, delta, tau, &p_vec_fd, 1e-9, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "phir_tt(" << comp_str << ", delta, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    delta = dat[i][test_data::rho_col]/param::rhoc[comp];
    err = fd2(memo2_phi_resi_tt, comp, delta, tau, &p_vec_fd, 1e-9, 1e-6, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  return 0;
}

uint test_sat_curve(comp_enum comp, test_data::data_set_enum data_set){
  std::vector<double> p_vec_fd;
  std::string comp_str = comp_enum_table[comp];
  std::vector< std::vector<double> > sat_liq_data, sat_vap_data;
  sort_sat(h2o, test_data::saturated_set, &sat_liq_data, &sat_vap_data);
  int err = 0;
  unsigned long i;
  double tau, pressure, delta;

  std::cout << "P_sat(" << comp_str << ", tau) ";
  for(i=0; i<sat_liq_data.size(); ++i){
    tau = param::Tc[comp]/sat_liq_data[i][test_data::T_col];
    pressure = sat_liq_data[i][test_data::P_col]*1000;
    err = fd1(sat_p, comp, tau, &p_vec_fd, 1e-8, pressure, 1e-3, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau_sat(" << comp_str << ", P) ";
  for(i=0; i<sat_liq_data.size(); ++i){
    tau = param::Tc[comp]/sat_liq_data[i][test_data::T_col];
    pressure = sat_liq_data[i][test_data::P_col]*1000;
    err = fd1(sat_tau, comp, pressure, &p_vec_fd, 1e-5, tau, 1e-3, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "rho_liq_sat(" << comp_str << ", tau) ";
  for(i=0; i<sat_liq_data.size(); ++i){
    tau = param::Tc[comp]/sat_liq_data[i][test_data::T_col];
    delta = sat_liq_data[i][test_data::rho_col]/param::rhoc[comp];
    err = fd1(sat_delta_l, comp, tau, &p_vec_fd, 1e-8, delta, 1, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;


  std::cout << "rho_vap_sat(" << comp_str << ", tau) ";
  for(i=0; i<sat_liq_data.size(); ++i){
    tau = param::Tc[comp]/sat_liq_data[i][test_data::T_col];
    delta = sat_vap_data[i][test_data::rho_col]/param::rhoc[comp];
    err = fd1(sat_delta_v, comp, tau, &p_vec_fd, 1e-8, delta, 1e-3, 0);
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;
  return 0;
}

uint test_vapor_state_var_change(comp_enum comp, test_data::data_set_enum data_set){
  std::vector<double> p_vec_fd;
  int err = 0;
  unsigned long i;
  double tau, pressure;
  std::vector< std::vector<double> > dat = read_data(comp, data_set);
  std::string comp_str = comp_enum_table[comp];

  std::cout << "delta_v(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_delta_vapor, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::rho_col]/param::rhoc[comp], 1e-2, 0);
    }
    else{ // stay loose at the critical point
      err = fd2(memo2_delta_vapor, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::rho_col]/param::rhoc[comp], 1, 0);
    }
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "h_v(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_enthalpy_vapor, comp, pressure, tau, &p_vec_fd, 1e-3, 1e-5, dat[i][test_data::h_col], 1e-2, 0);
    }
    else{
      err = fd2(memo2_enthalpy_vapor, comp, pressure, tau, &p_vec_fd, 1e-3, 1e-5, dat[i][test_data::h_col], 1e-1, 0);
    }
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "s_v(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.01 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.01){
      err = fd2(memo2_entropy_vapor, comp, pressure, tau, &p_vec_fd, pressure*1e-3, 1e-4, dat[i][test_data::s_col], 1e-2, 0);
    }
    else{
      err = fd2(memo2_entropy_vapor, comp, pressure, tau, &p_vec_fd, pressure*1e-3, 1e-4, dat[i][test_data::s_col], 1e-1 , 0);
    }
    if(err){
      std::cout << err;
      //err = fd2(memo2_entropy_vapor, comp, pressure, tau, &p_vec_fd, pressure*1e-3, 1e-4, dat[i][test_data::s_col], 1e-2, 1);
      // return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "u_v(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_internal_energy_vapor, comp, pressure, tau, &p_vec_fd, 1e-3, 1e-5, dat[i][test_data::u_col], 1e-2, 0);
    }
    else{
      err = fd2(memo2_internal_energy_vapor, comp, pressure, tau, &p_vec_fd, 1e-3, 1e-5, dat[i][test_data::u_col], 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau(" << comp_str << ", h, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_tau_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-2, 0);
    }
    else{
      err = fd2(memo2_tau_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau(" << comp_str << ", s, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.05 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.01){
      err = fd2(memo2_tau_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-1, 1e-5, tau, 1e-1, 0);
    }
    else{
      err = fd2(memo2_tau_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-1, 1e-5, tau, 1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau(" << comp_str << ", u, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.05 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.01){
      err = fd2(memo2_tau_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-2, 0);
    }
    else{
      err = fd2(memo2_tau_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "vf(" << comp_str << ", h, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_vf_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, 1.0, 1e-2, 0);
    }
    else{
      err = fd2(memo2_vf_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, 1.0, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "vf(" << comp_str << ", s, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_vf_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-4, 1e-5, 1.0, 1e-2, 0);
    }
    else{
      err = fd2(memo2_vf_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-4, 1e-5, 1.0, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "vf(" << comp_str << ", u, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_vf_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, 1.0, 1e-2, 0);
    }
    else{
      err = fd2(memo2_vf_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, 1.0, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  return 0;
}

uint test_liquid_state_var_change(comp_enum comp, test_data::data_set_enum data_set){
  std::vector<double> p_vec_fd;
  int err = 0;
  unsigned long i;
  double tau, pressure;
  std::vector< std::vector<double> > dat = read_data(comp, data_set);
  std::string comp_str = comp_enum_table[comp];

  std::cout << "delta_l(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.01 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.01){
      err = fd2(memo2_delta_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::rho_col]/param::rhoc[comp], 1e-2, 0);
    }
    else{ // stay loose at the critical point
      err = fd2(memo2_delta_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::rho_col]/param::rhoc[comp], 1, 0);
    }
    if(err){
      std::cout << err;
      err = fd2(memo2_delta_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::rho_col]/param::rhoc[comp], 1e-2, 1);
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "h_l(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_enthalpy_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::h_col], 1e-2, 0);
    }
    else{
      err = fd2(memo2_enthalpy_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::h_col], 1e-1, 0);
    }
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "s_l(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_entropy_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::s_col], 1e-2, 0);
    }
    else{
      err = fd2(memo2_entropy_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::s_col], 1e-1, 0);
    }
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "u_l(" << comp_str << ", P, tau) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_internal_energy_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::u_col], 1e-2, 0);
    }
    else{
      err = fd2(memo2_internal_energy_liquid, comp, pressure, tau, &p_vec_fd, 1e-4, 1e-5, dat[i][test_data::u_col], 1e-1, 0);
    }
    if(err){
      std::cout << err;
      return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau(" << comp_str << ", h, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_tau_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-2, 0);
    }
    else{
      err = fd2(memo2_tau_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau(" << comp_str << ", s, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.05 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.01){
      err = fd2(memo2_tau_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-2, 0);
    }
    else{
      err = fd2(memo2_tau_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "tau(" << comp_str << ", u, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.05 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.01){
      err = fd2(memo2_tau_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1e-2, 0);
    }
    else{
      err = fd2(memo2_tau_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, tau, 1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "vf(" << comp_str << ", h, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_vf_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, 0.0, 1e-2, 0);
    }
    else{
      err = fd2(memo2_vf_hp, comp, dat[i][test_data::h_col], pressure, &p_vec_fd, 1e-4, 1e-5, 0.0, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "vf(" << comp_str << ", s, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_vf_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-4, 1e-5, 0.0, 1e-2, 0);
    }
    else{
      err = fd2(memo2_vf_sp, comp, dat[i][test_data::s_col], pressure, &p_vec_fd, 1e-4, 1e-5, 0.0, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  std::cout << "vf(" << comp_str << ", u, P) ";
  for(i=0; i<dat.size(); ++i){
    tau = param::Tc[comp]/dat[i][test_data::T_col];
    pressure = dat[i][test_data::P_col]*1000;
    if (fabs(tau - 1.0) > 0.001 || fabs(pressure - param::Pc[comp])/param::Pc[comp] > 0.0001){
      err = fd2(memo2_vf_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, 0.0, 1e-2, 0);
    }
    else{
      err = fd2(memo2_vf_up, comp, dat[i][test_data::u_col], pressure, &p_vec_fd, 1e-4, 1e-5, 0.0, 1e-1, 0);
    }
    if(err){
      std::cout << err;
      //return err;
    }
    else{
      std::cout << ".";
    }
  }
  std::cout << "Passed" << std::endl;

  return 0;
}


int main(){
  std::vector<double> *p_vec_ptr;
  std::vector<double> p_vec_fd;
  int err = 0;
  int i;
  double tau, pressure, delta;

  std::cout << std::endl;
  std::cout << "Basic solver tests" << std::endl;
  std::cout << "----------------------------------" << std::endl;

  err = !test_bracket1(0);
  std::cout << "test_bracket1 passed: " << err << std::endl;

  err = !test_halley1(0);
  std::cout << "test_halley1 passed: " << err << std::endl;

  err = !test_newton_ls1(0);
  std::cout << "test_newton_ls1 passed: " << err << std::endl;

  std::cout << std::endl;
  std::cout << "Test basic h2o vapor properties" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_basic_properties(h2o, test_data::vapor_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test basic h2o liquid properties" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_basic_properties(h2o, test_data::liquid_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test basic h2o supercritical properties" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_basic_properties(h2o, test_data::supercritical_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test h2o sat curve" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::vector< std::vector<double> > sat_liq_data_h2o, sat_vap_data_h2o;
  err = test_sat_curve(h2o, test_data::saturated_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test h2o liquid state var change" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_liquid_state_var_change(h2o, test_data::liquid_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test h2o supercritical state var change" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_liquid_state_var_change(h2o, test_data::supercritical_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test h2o vapor state var change" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_vapor_state_var_change(h2o, test_data::vapor_set);
  if(err){
    exit(err);
  }

  std::cout << std::endl;
  std::cout << "Test basic r1234ze properties" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_basic_properties(r1234ze, test_data::mixed_set);
  if(err){
    exit(err);
  }

  return 0;
}

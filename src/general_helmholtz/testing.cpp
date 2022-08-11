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
  if(fabs(x1) < 1e-8) return fabs(x1 - x2) < 1e-8;
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

  std::vector<double> *yvec_ptr0, *yvec_ptr1;

  yvec_ptr0 = func(comp, x);
  yvec_ptr1 = func(comp, x + h);

  yvec_ptr->resize(3);
  yvec_ptr->at(0) = yvec_ptr0->at(0);
  yvec_ptr->at(1) = (yvec_ptr1->at(0) - yvec_ptr0->at(0))/h;
  yvec_ptr->at(2) = (yvec_ptr1->at(1) - yvec_ptr0->at(1))/h;

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

  yvec_ptr0 = func(comp, x1, x2);
  yvec_ptr1 = func(comp, x1 + h1, x2);
  yvec_ptr2 = func(comp, x1, x2 + h2);
  yvec_ptr1b = func(comp, x1 - h1, x2);
  yvec_ptr2b = func(comp, x1, x2 - h2);

  yvec_ptr->resize(6);
  yvec_ptr->at(f2) = yvec_ptr0->at(f2);
  yvec_ptr->at(f2_1) = (yvec_ptr1->at(f2) - yvec_ptr1b->at(f2))/h1/2;
  yvec_ptr->at(f2_2) = (yvec_ptr2->at(f2) - yvec_ptr2b->at(f2))/h2/2;
  yvec_ptr->at(f2_11) = (yvec_ptr1->at(f2_1) - yvec_ptr1b->at(f2_1))/h1/2;
  yvec_ptr->at(f2_12) = (yvec_ptr2->at(f2_1) - yvec_ptr2b->at(f2_1))/h2/2;
  yvec_ptr->at(f2_22) = (yvec_ptr2->at(f2_2) - yvec_ptr2b->at(f2_2))/h2/2;

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
    err = fd2(memo2_entropy, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, dat[i][test_data::s_col], 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
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
    err = fd2(memo2_enthalpy, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, dat[i][test_data::h_col], 1e-2, 0);
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
    err = fd2(memo2_internal_energy, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, dat[i][test_data::u_col], 1e-2, 0);
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
    err = fd2(memo2_isochoric_heat_capacity, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, dat[i][test_data::cv_col], 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
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
    err = fd2(memo2_isobaric_heat_capacity, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, dat[i][test_data::cp_col], 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
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
    err = fd2(memo2_speed_of_sound, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, dat[i][test_data::w_col], 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
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
    err = fd2(memo2_gibbs, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, dat[i][test_data::h_col] - dat[i][test_data::T_col]*dat[i][test_data::s_col], 1e-1, 0);
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
    err = fd2(memo2_helmholtz, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, dat[i][test_data::u_col] - dat[i][test_data::T_col]*dat[i][test_data::s_col], 1e-1, 0);
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
    err = fd2(memo2_phi_ideal, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_ideal_d, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_ideal_t, comp, delta, tau, &p_vec_fd, 1e-11, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_ideal_dd, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_ideal_tt, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_resi, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_resi_d, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_resi_t, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_resi_dd, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
    if(err){
      std::cout << err;
      return err;
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
    err = fd2(memo2_phi_resi_dt, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
    err = fd2(memo2_phi_resi_tt, comp, delta, tau, &p_vec_fd, 1e-10, 1e-5, nan("no check"), 1e-2, 0);
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
  std::cout << "Test basic r1234ze properties" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  err = test_basic_properties(r1234ze, test_data::mixed_set);
  if(err){
    exit(err);
  }

/*

    err = !fd2(memo2_delta_liquid, comp_enum::h2o, 99.2418352, 647.096/300.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_delta_liquid passed: " << err << std::endl;

    err = !fd2(memo2_enthalpy_liquid, comp_enum::h2o, 99.2418352, 647.096/300.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_enthalpy_liquid passed: " << err << std::endl;

    err = !fd2(memo2_entropy_liquid, comp_enum::h2o, 99.2418352, 647.096/300.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_entropy_liquid passed: " << err << std::endl;

    err = !fd2(memo2_internal_energy_liquid, comp_enum::h2o, 50, 647.096/310.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_internal_energy_liquid passed: " << err << std::endl;

    err = !fd2(memo2_delta_vapor, comp_enum::h2o, 99.9679423, 647.096/500.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_delta_vapor passed: " << err << std::endl;

    err = !fd2(memo2_enthalpy_vapor, comp_enum::h2o, 99.9679423, 647.096/500.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_enthalpy_vapor passed: " << err << std::endl;

    err = !fd2(memo2_entropy_vapor, comp_enum::h2o, 99.9679423, 647.096/500.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_entropy_vapor passed: " << err << std::endl;

    err = !fd2(memo2_internal_energy_vapor, comp_enum::h2o, 99.9679423, 647.096/500.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_internal_energy_vapor passed: " << err << std::endl;

    err = !fd1(sat_p, comp_enum::h2o, 647.096/450, &p_vec_fd, 1e-6, 0);
    std::cout << "sat_p passed: " << err << std::endl;

    err = !fd1(sat_delta_l, comp_enum::h2o, 647.096/450, &p_vec_fd, 1e-6, 0);
    std::cout << "sat_delta_l passed: " << err << std::endl;

    err = !fd1(sat_delta_v, comp_enum::h2o, 647.096/450, &p_vec_fd, 1e-6, 0);
    std::cout << "sat_delta_v passed: " << err << std::endl;

    err = !fd1(sat_tau, comp_enum::h2o, 932.203564, &p_vec_fd, 1e-6, 0);
    std::cout << "sat_tau passed: " << err << std::endl;

    err = !fd2(memo2_tau_hp, comp_enum::h2o, 154.406, 50.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_tau_hp (liquid) passed: " << err << std::endl;

    err = !fd2(memo2_tau_hp, comp_enum::h2o, 1000, 932.203564, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_tau_hp passed (two-phase): " << err << std::endl;

    err = !fd2(memo2_tau_sp, comp_enum::h2o, 0.5301, 50.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_tau_sp passed (liquid): " << err << std::endl;

    err = !fd2(memo2_tau_sp, comp_enum::h2o, 5.0, 932.203564, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_tau_sp passed (two-phase): " << err << std::endl;

    err = !fd2(memo2_tau_up, comp_enum::h2o, 154.355, 50.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_tau_up passed (liquid): " << err << std::endl;

    err = !fd2(memo2_tau_up, comp_enum::h2o, 1200.0, 932.22, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_tau_up passed (two-phase): " << err << std::endl;

    err = !fd2(memo2_vf_hp, comp_enum::h2o, 154.406, 50.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_vf_hp (liquid) passed: " << err << std::endl;

    err = !fd2(memo2_vf_hp, comp_enum::h2o, 1000, 932.203564, &p_vec_fd, 0.01, 0);
    std::cout << "memo2_vf_hp passed (two-phase): " << err << std::endl;

    err = !fd2(memo2_vf_sp, comp_enum::h2o, 0.5301, 50.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_vf_sp passed (liquid): " << err << std::endl;

    err = !fd2(memo2_vf_sp, comp_enum::h2o, 5.0, 932.203564, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_vf_sp passed (two-phase): " << err << std::endl;

    err = !fd2(memo2_vf_up, comp_enum::h2o, 154.355, 50.0, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_vf_up passed (liquid): " << err << std::endl;

    err = !fd2(memo2_vf_up, comp_enum::h2o, 1200.0, 932.22, &p_vec_fd, 1e-4, 0);
    std::cout << "memo2_vf_up passed (two-phase): " << err << std::endl;


  for(tau=param::Tc[r1234ze]/param::Tt[r1234ze]; tau >= 1.0; tau -= 0.01){
    std::cout << param::Tc[r1234ze]/tau << "\t";
    std::cout << sat_delta_l(comp_enum::r1234ze, tau)->at(0)*param::rhoc[r1234ze] << "\t";
    std::cout << sat_delta_v(comp_enum::r1234ze, tau)->at(0)*param::rhoc[r1234ze] << "\t";
    std::cout << sat_p(comp_enum::r1234ze, tau)->at(0) << std::endl;
  }
 */


  return 0;
}

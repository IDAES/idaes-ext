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

#include"phi.h"
#include"props.h"
#include"sat.h"
#include"testing.h"
#include <iostream>
#include <math.h>

inline bool rel_same(double x1, double x2, double tol){
  if(fabs(x1) < 1e-10) return fabs(x1 - x2) < tol;
  return fabs((x1 - x2)/x1) < tol;
}

int fd1(test_fptr1 func, comp_enum comp, double x, std::vector<double> *yvec_ptr, double h, bool dbg){
  std::vector<double> *yvec_ptr0, *yvec_ptr1;

  yvec_ptr0 = func(comp, x);
  yvec_ptr1 = func(comp, x + h);

  yvec_ptr->resize(3);
  yvec_ptr->at(0) = yvec_ptr0->at(0);
  yvec_ptr->at(1) = (yvec_ptr1->at(0) - yvec_ptr0->at(0))/h;
  yvec_ptr->at(2) = (yvec_ptr1->at(1) - yvec_ptr0->at(1))/h;

  if(dbg){
    std::cout << "f = " << yvec_ptr0->at(0) << std::endl;
    std::cout << "f_d = " << yvec_ptr0->at(1) << " f.d. approx = " << yvec_ptr->at(1) << std::endl;
    std::cout << "f_t = " << yvec_ptr0->at(2) << " f.d. approx = " << yvec_ptr->at(2) << std::endl;
  }

  if(!rel_same(yvec_ptr->at(1), yvec_ptr0->at(1), 1e-6)) return 1;
  if(!rel_same(yvec_ptr->at(2), yvec_ptr0->at(2), 1e-6)) return 1;

  return 0;
}

int fd2(test_fptr2 func, comp_enum comp, double x1, double x2, std::vector<double> *yvec_ptr, double h, bool dbg){
  std::vector<double> *yvec_ptr0, *yvec_ptr1, *yvec_ptr2;

  yvec_ptr0 = func(comp, x1, x2);
  yvec_ptr1 = func(comp, x1 + h, x2);
  yvec_ptr2 = func(comp, x1, x2 + h);

  yvec_ptr->resize(6);
  yvec_ptr->at((uint)deriv2_enum::f) = yvec_ptr0->at((uint)deriv2_enum::f);
  yvec_ptr->at((uint)deriv2_enum::f_d) =
    (yvec_ptr1->at((uint)deriv2_enum::f) - yvec_ptr0->at((uint)deriv2_enum::f))/h;
  yvec_ptr->at((uint)deriv2_enum::f_t) =
    (yvec_ptr2->at((uint)deriv2_enum::f) - yvec_ptr0->at((uint)deriv2_enum::f))/h;
  yvec_ptr->at((uint)deriv2_enum::f_dd) =
    (yvec_ptr1->at((uint)deriv2_enum::f_d) - yvec_ptr0->at((uint)deriv2_enum::f_d))/h;
  yvec_ptr->at((uint)deriv2_enum::f_dt) =
    (yvec_ptr2->at((uint)deriv2_enum::f_d) - yvec_ptr0->at((uint)deriv2_enum::f_d))/h;
  yvec_ptr->at((uint)deriv2_enum::f_tt) =
    (yvec_ptr2->at((uint)deriv2_enum::f_t) - yvec_ptr0->at((uint)deriv2_enum::f_t))/h;

  if(dbg){
    std::cout << "f = " << yvec_ptr0->at(0) << std::endl;
    std::cout << "f_d = " << yvec_ptr0->at((uint)deriv2_enum::f_d) << " f.d. approx = " << yvec_ptr->at((uint)deriv2_enum::f_d) << std::endl;
    std::cout << "f_t = " << yvec_ptr0->at((uint)deriv2_enum::f_t) << " f.d. approx = " << yvec_ptr->at((uint)deriv2_enum::f_t) << std::endl;
    std::cout << "f_dd = " << yvec_ptr0->at((uint)deriv2_enum::f_dd) << " f.d. approx = " << yvec_ptr->at((uint)deriv2_enum::f_dd) << std::endl;
    std::cout << "f_dt = " << yvec_ptr0->at((uint)deriv2_enum::f_dt) << " f.d. approx = " << yvec_ptr->at((uint)deriv2_enum::f_dt) << std::endl;
    std::cout << "f_tt = " << yvec_ptr0->at((uint)deriv2_enum::f_tt) << " f.d. approx = " << yvec_ptr->at((uint)deriv2_enum::f_tt) << std::endl;
  }

  if (!rel_same(yvec_ptr0->at((uint)deriv2_enum::f_d), yvec_ptr->at((uint)deriv2_enum::f_d), 1e-4)) return 1;
  if (!rel_same(yvec_ptr0->at((uint)deriv2_enum::f_t), yvec_ptr->at((uint)deriv2_enum::f_t), 1e-4)) return 1;
  if (!rel_same(yvec_ptr0->at((uint)deriv2_enum::f_dd), yvec_ptr->at((uint)deriv2_enum::f_dd), 1e-4)) return 1;
  if (!rel_same(yvec_ptr0->at((uint)deriv2_enum::f_dt), yvec_ptr->at((uint)deriv2_enum::f_dt), 1e-4)) return 1;
  if (!rel_same(yvec_ptr0->at((uint)deriv2_enum::f_tt), yvec_ptr->at((uint)deriv2_enum::f_tt), 1e-4)) return 1;
  return 0;
}


int main(){
  std::vector<double> *p_vec_ptr;
  std::vector<double> p_vec_fd;
  int err = 0;

  err = !fd2(memo2_pressure, comp_enum::h2o, 838.025/322.0, 647.096/500.0, &p_vec_fd, 1e-9, 1);
  std::cout << "memo2_pressure passed: " << err << std::endl;

  err = !fd2(memo2_internal_energy, comp_enum::h2o, 838.025/322.0, 647.096/500.0, &p_vec_fd, 1e-9, 1);
  std::cout << "memo2_internal_energy passed: " << err << std::endl;

  err = !fd2(memo2_entropy, comp_enum::h2o, 838.025/322.0, 647.096/500.0, &p_vec_fd, 1e-9, 1);
  std::cout << "memo2_entropy passed: " << err << std::endl;

  err = !fd2(memo2_enthalpy, comp_enum::h2o, 838.025/322.0, 647.096/500.0, &p_vec_fd, 1e-9, 1);
  std::cout << "memo2_enthalpy passed: " << err << std::endl;


  std::cout << "sat" << std::endl;
  std::vector<double> *delta_l_vec_ptr, *delta_v_vec_ptr;
  std::vector<double> delta_l_fd_ptr, delta_v_fd_ptr, p_fd_ptr;
  p_vec_ptr = sat_p(comp_enum::h2o, 647.096/450);
  delta_l_vec_ptr = sat_delta_l(comp_enum::h2o, 647.096/450);
  delta_v_vec_ptr = sat_delta_v(comp_enum::h2o, 647.096/450);
  fd1(sat_p, comp_enum::h2o, 647.096/450, &p_fd_ptr, 1e-9, 0);
  fd1(sat_delta_l, comp_enum::h2o, 647.096/450, &delta_l_fd_ptr, 1e-8, 0);
  fd1(sat_delta_v, comp_enum::h2o, 647.096/450, &delta_v_fd_ptr, 1e-8, 0);

  std::cout << "p = " << p_vec_ptr->at(0) << std::endl;
  std::cout << "p_t = " << p_vec_ptr->at(1) << " fd: " << p_fd_ptr.at(1) << std::endl;
  std::cout << "p_tt = " << p_vec_ptr->at(2) << " fd: " << p_fd_ptr.at(2) << std::endl;

  std::cout << "rho_v " << 322*delta_v_vec_ptr->at(0) << std::endl;
  std::cout << "delta_v_t = " << delta_v_vec_ptr->at(1) << " fd: " << delta_v_fd_ptr.at(1) << std::endl;
  std::cout << "delta_v_tt = " << delta_v_vec_ptr->at(2) << " fd: " << delta_v_fd_ptr.at(2) << std::endl;

  std::cout << "rho_l " << 322*delta_l_vec_ptr->at(0) << std::endl;
  std::cout << "delta_l_t = " << delta_l_vec_ptr->at(1) << " fd: " << delta_l_fd_ptr.at(1) << std::endl;
  std::cout << "delta_l_tt = " << delta_l_vec_ptr->at(2) << " fd: " << delta_l_fd_ptr.at(2) << std::endl;
  /*
  double temperature;
  for(temperature=250; temperature < 647.096; temperature += 1){
    std::cout << temperature << "\t";
    std::cout << pressure(comp_enum::h2o, 900.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 800.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 700.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 600.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 500.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 400.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 300.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 200.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 100.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 10.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 5.0/322.0, 647.096/temperature) << "\t";
    std::cout << pressure(comp_enum::h2o, 1/322.0, 647.096/temperature) << std::endl;
  }
  */
  return 0;
}

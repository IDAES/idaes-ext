#include"phi.h"
#include"props.h"
#include"sat.h"
#include <iostream>

int main(){
  std::cout << "Ideal" << std::endl;
  std::vector<double> *phi_ideal_ptr = phi_ideal(comp_enum::h2o, 838.025/322.0, 647.096/500.0);
  std::cout << "f(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f) << std::endl;
  std::cout << "f_d(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_d) << std::endl;
  std::cout << "f_dd(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_dd) << std::endl;
  std::cout << "f_ddd(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_ddd) << std::endl;
  std::cout << "f_dddd(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_dddd) << std::endl;
  std::cout << "f_t(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_t) << std::endl;
  std::cout << "f_dt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_dt) << std::endl;
  std::cout << "f_ddt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_ddt) << std::endl;
  std::cout << "f_dddt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_dddt) << std::endl;
  std::cout << "f_tt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_tt) << std::endl;
  std::cout << "f_dtt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_dtt) << std::endl;
  std::cout << "f_ddtt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_ddtt) << std::endl;
  std::cout << "f_ttt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_ttt) << std::endl;
  std::cout << "f_dttt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_dttt) << std::endl;
  std::cout << "f_tttt(x) = " << phi_ideal_ptr->at((unsigned int)deriv4_enum::f_tttt) << std::endl;

  std::cout << "Real" << std::endl;
  std::vector<double> *phi_real_ptr = phi_real(comp_enum::h2o, 838.025/322.0, 647.096/500.0);
  std::cout << "f(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f) << std::endl;
  std::cout << "f_d(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_d) << std::endl;
  std::cout << "f_dd(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_dd) << std::endl;
  std::cout << "f_ddd(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_ddd) << std::endl;
  std::cout << "f_dddd(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_dddd) << std::endl;
  std::cout << "f_t(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_t) << std::endl;
  std::cout << "f_dt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_dt) << std::endl;
  std::cout << "f_ddt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_ddt) << std::endl;
  std::cout << "f_dddt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_dddt) << std::endl;
  std::cout << "f_tt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_tt) << std::endl;
  std::cout << "f_dtt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_dtt) << std::endl;
  std::cout << "f_ddtt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_ddtt) << std::endl;
  std::cout << "f_ttt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_ttt) << std::endl;
  std::cout << "f_dttt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_dttt) << std::endl;
  std::cout << "f_tttt(x) = " << phi_real_ptr->at((unsigned int)deriv4_enum::f_tttt) << std::endl;

  std::cout << "Pressure" << std::endl;
  std::cout << pressure(comp_enum::h2o, 838.025/322.0, 647.096/500.0) << std::endl;

  std::cout << "Pressure2" << std::endl;
  std::vector<double> p;
  std::vector<double> p_d;
  std::vector<double> p_t;
  double eps = 1e-8;
  pressure2(comp_enum::h2o, 838.025/322.0, 647.096/500.0, &p);
  pressure2(comp_enum::h2o, 838.025/322.0 + eps, 647.096/500.0, &p_d);
  pressure2(comp_enum::h2o, 838.025/322.0, 647.096/500.0 + eps, &p_t);
  double f_d = (p_d.at((unsigned int)deriv2_enum::f) - p.at((unsigned int)deriv2_enum::f))/eps;
  double f_t = (p_t.at((unsigned int)deriv2_enum::f) - p.at((unsigned int)deriv2_enum::f))/eps;
  double f_dd = (p_d.at((unsigned int)deriv2_enum::f_d) - p.at((unsigned int)deriv2_enum::f_d))/eps;
  double f_dt = (p_d.at((unsigned int)deriv2_enum::f_t) - p.at((unsigned int)deriv2_enum::f_t))/eps;
  double f_tt = (p_t.at((unsigned int)deriv2_enum::f_t) - p.at((unsigned int)deriv2_enum::f_t))/eps;
  std::cout << "f(x) = " << p.at((unsigned int)deriv2_enum::f) << std::endl;
  std::cout << "f_d(x) = " << p.at((unsigned int)deriv2_enum::f_d) << " " <<  f_d << std::endl;
  std::cout << "f_dd(x) = " << p.at((unsigned int)deriv2_enum::f_dd) << " " <<  f_dd << std::endl;
  std::cout << "f_t(x) = " << p.at((unsigned int)deriv2_enum::f_t) << " " <<  f_t << std::endl;
  std::cout << "f_td(x) = " << p.at((unsigned int)deriv2_enum::f_dt) << " " <<  f_dt << std::endl;
  std::cout << "f_tt(x) = " << p.at((unsigned int)deriv2_enum::f_tt) << " " <<  f_tt << std::endl;

  std::cout << "phir for sat" << std::endl;
  std::vector<double> phir;
  phi_real_for_sat(comp_enum::h2o, 900.025/322.0, 647.096/400.0, &phir);
  phi_real_for_sat(comp_enum::h2o, 838.025/322.0, 647.096/500.0, &phir);
  std::cout << "f(x) = " << phir.at(0) << std::endl;
  std::cout << "f_d(x) = " << phir.at(1) << std::endl;
  std::cout << "f_dd(x) = " << phir.at(2) << std::endl;


  std::cout << "sat" << std::endl;
  double delta_l, delta_v;

  std::cout << sat(comp_enum::h2o, 647.096/450, &delta_l, &delta_v) << std::endl;
  std::cout << "rho_v = " << 322*delta_v << std::endl;
  std::cout << "rho_l = " << 322*delta_l << std::endl;

  return 0;
}

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

/*--------------------------------------------------------------------------------
 Provide functions to solve for liquid and vapor density from T and P

 Author: John Eslick
 File: delta.cpp
--------------------------------------------------------------------------------*/
#include<math.h>
#include<iostream>
#include<boost/math/tools/roots.hpp>
#include"props.h"
#include"sat.h"
#include"delta.h"

prop_memo_table22 memo_table_delta_liquid2;
prop_memo_table22 memo_table_delta_vapor2;

class pfunctor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _tau;
public:
    pfunctor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_tau(double t){
      this->_tau = t;
    }
    std::tuple<double, double, double>  operator () (double delta) {
      f22_struct out;
      pressure2(this->_comp, delta, this->_tau, &out);
      return std::make_tuple(
        (out.f - this->_pressure)/cdata[this->_comp].Pc,
        out.f_1/cdata[this->_comp].Pc,
        out.f_11/cdata[this->_comp].Pc
      );
    }
};

double delta_vapor(uint comp, double pr, double tau){
  double delta_sat, delta, delta_guess=0.0;
  double p_sat;
  parameters_struct *dat = &cdata[comp];

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=40;
  int digits = std::numeric_limits<double>::digits - 4;
  pfunctor_deriv fgh = pfunctor_deriv(comp); fgh.set_pressure(pr); fgh.set_tau(tau); 

  if(pr >= dat->Pc){ // over ciritical pressure
    if(tau <= tau_c(comp)){ // supercritical or liquid (there is no vapor over critial pressure)
      return delta_liquid(comp, pr, tau);
    }
  }
  if(tau <= tau_c(comp)){
    p_sat = dat->Pc;
    delta_sat = delta_c(comp);
  }
  else{
    p_sat = sat_p(comp, tau).f;
    delta_sat = sat_delta_v(comp, tau).f;
  }
  try{
    delta_guess = pr/dat->R/(dat->T_star/tau)/dat->rhoc;
    if(delta_guess > delta_sat){
      delta_guess = delta_sat;
    }
    delta = halley_iterate(fgh, delta_guess, 1e-8, delta_sat, digits, h_it_max);
    return delta;
  }
  catch(...){
    try{ // edgy cases
      delta = halley_iterate(fgh, delta_sat, delta_sat*0.75, delta_sat*1.25, digits, h_it_max);
      return delta;
    }
    catch(...){
      delta = delta_guess;
      std::cout << "delta_vapor(p, t) solve failed" << std::endl;
      std::cout << "res: " << std::get<0>(fgh(delta)) << " delta: " << delta;
      std::cout << " delta_guess: " << delta_guess << " delta_sat " << delta_sat;
      std::cout << " tau: " << tau << " p:" << pr << std::endl;
      return delta_guess;
    }
  }
  return delta_sat;
}

double delta_liquid(uint comp, double pr, double tau){
  double delta_sat, delta;
  double p_sat;
  parameters_struct *dat = &cdata[comp];

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=40;
  int digits = std::numeric_limits<double>::digits - 2;
  pfunctor_deriv fgh = pfunctor_deriv(comp); fgh.set_pressure(pr); fgh.set_tau(tau); 

  if(pr >= dat->Pc && tau <= tau_c(comp)){ // super critical
    try{
      delta = halley_iterate(fgh, dat->rho_max/dat->rho_star*0.99, 1e-7, dat->rho_max/dat->rho_star, digits, h_it_max);
      return delta;
    }
    catch(...){
      delta = delta_c(comp);
      std::cout << "delta_liquid(p, t) solve failed in supercritical region" << std::endl;
      std::cout << "res: " << std::get<0>(fgh(delta)) << " delta: " << delta << " tau: " << tau << " p:" << pr << std::endl;
      return delta;
    }
  }
  if(tau <= tau_c(comp)){
    delta_sat = delta_c(comp);
    p_sat = dat->Pc;
  }
  else{
    delta_sat = sat_delta_l(comp, tau).f;
    p_sat = sat_p(comp, tau).f;
  }
  try{ 
    delta = halley_iterate(fgh, dat->rho_max/dat->rho_star, delta_sat, dat->rho_max/dat->rho_star, digits, h_it_max);
    return delta;
  }
  catch(...){
    try{  // edgy cases
      delta = halley_iterate(fgh, delta_sat, 0.5*delta_sat, delta_sat, digits, h_it_max);
      return delta;
    }
    catch(...){
      std::cout << "delta_liquid(p, t) solve failed" << std::endl;
      std::cout << "res: " << std::get<0>(fgh(delta)) << " delta: " << delta << " tau: " << tau << " p:" << pr << std::endl;
    }
  }
  return delta_sat;
}

void delta_liquid2(uint comp, double pr, double tau, f22_struct *out){
  double delta_l = delta_liquid(comp, pr, tau);
  f22_struct pr_vec = memo2_pressure(comp, delta_l, tau); // get derivatives
  out->f = delta_l;
  out->f_1 = 1.0/pr_vec.f_1;
  out->f_2 = -pr_vec.f_2/pr_vec.f_1;
  out->f_11 = -pr_vec.f_11*out->f_1*out->f_1*out->f_1;
  out->f_12 = -(pr_vec.f_12 + pr_vec.f_11*out->f_2)*out->f_1*out->f_1;
  out->f_22 = -(out->f_1*(pr_vec.f_22 + out->f_2*pr_vec.f_12) + pr_vec.f_2*out->f_12);
}

void delta_vapor2(uint comp, double pr, double tau, f22_struct *out){
  double delta_v = delta_vapor(comp, pr, tau);
  f22_struct pr_vec = memo2_pressure(comp, delta_v, tau); // get derivatives
  out->f = delta_v;
  out->f_1 = 1.0/pr_vec.f_1;
  out->f_2 = -pr_vec.f_2/pr_vec.f_1;
  out->f_11 = -pr_vec.f_11*out->f_1*out->f_1*out->f_1;
  out->f_12 = -(pr_vec.f_12 + pr_vec.f_11*out->f_2)*out->f_1*out->f_1;
  out->f_22 = -(out->f_1*(pr_vec.f_22 + out->f_2*pr_vec.f_12) + pr_vec.f_2*out->f_12);
}

f22_struct memo2_delta_liquid(uint comp, double pr, double tau){
  try{
    return memo_table_delta_liquid2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range const&){
  }
  if(memo_table_delta_liquid2.size() > MAX_MEMO_PROP) memo_table_delta_liquid2.clear();
  f22_struct *yvec_ptr = &memo_table_delta_liquid2[std::make_tuple(comp, pr, tau)];
  delta_liquid2(comp, pr, tau, yvec_ptr);
  return *yvec_ptr;
}

f22_struct memo2_delta_vapor(uint comp, double pr, double tau){
  try{
    return memo_table_delta_vapor2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range const&){
  }
  if(memo_table_delta_vapor2.size() > MAX_MEMO_PROP) memo_table_delta_vapor2.clear();
  f22_struct *yvec_ptr = &memo_table_delta_vapor2[std::make_tuple(comp, pr, tau)];
  delta_vapor2(comp, pr, tau, yvec_ptr);
  return *yvec_ptr;
}

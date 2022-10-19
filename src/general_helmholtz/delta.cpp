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

class pfunctor
{
private:
    uint _comp;
    double _pressure;
    double _tau;
public:
    pfunctor(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_tau(double t){
      this->_tau = t;
    }
    double operator () (double delta) {
        return (pressure(this->_comp, delta, this->_tau) - this->_pressure)/cdata[this->_comp].Pc;
    }
};

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
  double delta_sat;
  double p_sat;
  parameters_struct *dat = &cdata[comp];

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=80;
  int digits = std::numeric_limits<double>::digits - 4;
  pfunctor f = pfunctor(comp); f.set_pressure(pr); f.set_tau(tau);
  pfunctor_deriv fgh = pfunctor_deriv(comp); fgh.set_pressure(pr); fgh.set_tau(tau); 

  if(pr > dat->Pc){ // over ciritical pressure
    std::vector<double> out;
    if(tau >= tau_c(comp)){ // supercritical
      return delta_liquid(comp, pr, tau);
    }
    else{ // high pressure vapor
      double delta_guess = pr/dat->R/(dat->T_star/tau)/dat->rhoc;
      return halley_iterate(fgh, delta_guess, 1e-6, delta_c(comp) * 1.001, digits, h_it_max);
    }
  }
  delta_sat = sat_delta_v(comp, tau).f;
  p_sat = sat_p(comp, tau).f;
  if(pr <= p_sat){ // vapor
    double delta_guess = pr/dat->R/(dat->T_star/tau)/dat->rhoc;
    return halley_iterate(fgh, delta_guess, 1e-8, delta_sat, digits, h_it_max);
  }
  try{ // liquid or right on the edge
    return halley_iterate(fgh, delta_sat, delta_sat*0.5, delta_sat*2.0, digits, h_it_max);
  }
  catch(...){
  }
  return delta_sat;
}

double delta_liquid(uint comp, double pr, double tau){
  double delta_sat;
  double p_sat;
  parameters_struct *dat = &cdata[comp];

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=80;
  int digits = std::numeric_limits<double>::digits - 4;
  pfunctor f = pfunctor(comp); f.set_pressure(pr); f.set_tau(tau);
  pfunctor_deriv fgh = pfunctor_deriv(comp); fgh.set_pressure(pr); fgh.set_tau(tau); 

  if(pr >= dat->Pc && tau <= tau_c(comp)){ // super critical
    return halley_iterate(fgh, delta_c(comp), 1e-7, dat->rho_max/dat->rho_star, digits, h_it_max);
  }
  delta_sat = sat_delta_l(comp, tau).f;
  p_sat = sat_p(comp, tau).f;
  if(pr >= p_sat){ // liquid
    return halley_iterate(fgh, delta_sat, delta_sat*0.99, dat->rho_max/dat->rho_star, digits, h_it_max);
  }
  try{ // vapor or right on the edge
    //std::cout << "delta sat " << delta_sat << std::endl;
    double delta = halley_iterate(fgh, delta_sat, delta_sat*0.8, delta_sat*1.2, digits, h_it_max);
    //std::cout << "delta = " << delta << " residual " << std::get<0>(fgh(delta)) << std::endl;
    return delta;
  }
  catch(...){
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

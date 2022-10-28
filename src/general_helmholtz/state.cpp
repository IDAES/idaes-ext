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

/*------------------------------------------------------------------------------
Author: John Eslick
File state.cpp

Functions to enable change of state variables. The end goal of this section is
given a set of state varible, calculate T, P, and vapor fraction.  From there
you can calculate delta for each phase and then all the rest of the properties
can be calculated.

The general method for changing state varaibles is (first step is here):
  1) tau = tau(v1, v2), vf=vf(v1, v2), p(v1, v2); so far p is always a state
     variable but it doesn't need to be that case we can add on to support more
  2) delta_v = delta_v(p, tau), delta_l = delta_l(p, tau)
  3) property_v = f(delta_v, tau), property_l = f(delta_l, tau)
  4) calculate mixed phase properies
------------------------------------------------------------------------------*/

#include <boost/math/tools/roots.hpp>
#include "config.h"
#include "state.h"
#include "props.h"
#include "delta.h"
#include "sat.h"
#include <iostream>

struct state_solve_dat{ // Solver wrapper data structure for const args
  uint comp;
  double p;
  double h;
};

/*------------------------------------------------------------------------------
  Memoization tables
------------------------------------------------------------------------------*/

prop_memo_table22 memo_table_enthalpy_vapor2;
prop_memo_table22 memo_table_entropy_vapor2;
prop_memo_table22 memo_table_internal_energy_vapor2;
prop_memo_table22 memo_table_enthalpy_liquid2;
prop_memo_table22 memo_table_entropy_liquid2;
prop_memo_table22 memo_table_internal_energy_liquid2;
prop_memo_table22 memo_table_tau_hp2;
prop_memo_table22 memo_table_tau_sp2;
prop_memo_table22 memo_table_tau_up2;
prop_memo_table22 memo_table_vf_hp2;
prop_memo_table22 memo_table_vf_sp2;
prop_memo_table22 memo_table_vf_up2;

/*------------------------------------------------------------------------------
  Vapor enthalpy as a function of pressure and tau, used to solve T(h, p)
------------------------------------------------------------------------------*/

void enthalpy_vapor2(uint comp, double pr, double tau, f22_struct *out){
  f22_struct delta_v_vec, h_vec;
  delta_vapor2(comp, pr, tau, &delta_v_vec);
  enthalpy2(comp, delta_v_vec.f, tau, &h_vec);
  out->f = h_vec.f;
  out->f_1 = h_vec.f_1*delta_v_vec.f_1;
  out->f_2 = h_vec.f_2 + delta_v_vec.f_2*h_vec.f_1;
  out->f_11 =
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_1 +
    h_vec.f_1*delta_v_vec.f_11;
  out->f_12 =
    h_vec.f_12*delta_v_vec.f_1 +
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_12;
  out->f_22 =
    h_vec.f_22 +
    2*h_vec.f_12*delta_v_vec.f_2 +
    h_vec.f_11*delta_v_vec.f_2*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_22;
}

/*------------------------------------------------------------------------------
  Vapor entropy as a function of pressure and tau, used to solve T(s, p)
------------------------------------------------------------------------------*/

void entropy_vapor2(uint comp, double pr, double tau, f22_struct *out){
  f22_struct delta_v_vec, h_vec;
  delta_v_vec = memo2_delta_vapor(comp, pr, tau);
  h_vec = memo2_entropy(comp, delta_v_vec.f, tau);
  out->f = h_vec.f;
  out->f_1 = h_vec.f_1*delta_v_vec.f_1;
  out->f_2 = h_vec.f_2 + delta_v_vec.f_2*h_vec.f_1;
  out->f_11 =
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_1 +
    h_vec.f_1*delta_v_vec.f_11;
  out->f_12 =
    h_vec.f_12*delta_v_vec.f_1 +
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_12;
  out->f_22 =
    h_vec.f_22 +
    2*h_vec.f_12*delta_v_vec.f_2 +
    h_vec.f_11*delta_v_vec.f_2*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_22;
}

/*------------------------------------------------------------------------------
  Vapor internal energy as a function of pressure and tau, used to solve T(u, p)
------------------------------------------------------------------------------*/

void internal_energy_vapor2(uint comp, double pr, double tau, f22_struct *out){
  f22_struct delta_v_vec, h_vec;
  delta_v_vec = memo2_delta_vapor(comp, pr, tau);
  h_vec = memo2_internal_energy(comp, delta_v_vec.f, tau);
  out->f = h_vec.f;
  out->f_1 = h_vec.f_1*delta_v_vec.f_1;
  out->f_2 = h_vec.f_2 + delta_v_vec.f_2*h_vec.f_1;
  out->f_11 =
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_1 +
    h_vec.f_1*delta_v_vec.f_11;
  out->f_12 =
    h_vec.f_12*delta_v_vec.f_1 +
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_12;
  out->f_22 =
    h_vec.f_22 +
    2*h_vec.f_12*delta_v_vec.f_2 +
    h_vec.f_11*delta_v_vec.f_2*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_22;
}

/*------------------------------------------------------------------------------
  Vapor enthalpy as a function of pressure and tau, used to solve T(h, p)
------------------------------------------------------------------------------*/

double enthalpy_liquid(uint comp, double pr, double tau){
  double delta_l = delta_liquid(comp, pr, tau);
  return enthalpy(comp, delta_l, tau);
}

void enthalpy_liquid2(uint comp, double pr, double tau, f22_struct *out){
  f22_struct delta_v_vec, h_vec;
  delta_liquid2(comp, pr, tau, &delta_v_vec);
  enthalpy2(comp, delta_v_vec.f, tau, &h_vec);
  out->f = h_vec.f;
  out->f_1 = h_vec.f_1*delta_v_vec.f_1;
  out->f_2 = h_vec.f_2 + delta_v_vec.f_2*h_vec.f_1;
  out->f_11 =
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_1 +
    h_vec.f_1*delta_v_vec.f_11;
  out->f_12 =
    h_vec.f_12*delta_v_vec.f_1 +
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_12;
  out->f_22 =
    h_vec.f_22 +
    2*h_vec.f_12*delta_v_vec.f_2 +
    h_vec.f_11*delta_v_vec.f_2*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_22;
}

/*------------------------------------------------------------------------------
  Liquid entropy as a function of pressure and tau, used to solve T(s, p)
------------------------------------------------------------------------------*/

void entropy_liquid2(uint comp, double pr, double tau, f22_struct *out){
  f22_struct delta_v_vec, h_vec;
  delta_v_vec = memo2_delta_liquid(comp, pr, tau);
  h_vec = memo2_entropy(comp, delta_v_vec.f, tau);
  out->f = h_vec.f;
  out->f_1 = h_vec.f_1*delta_v_vec.f_1;
  out->f_2 = h_vec.f_2 + delta_v_vec.f_2*h_vec.f_1;
  out->f_11 =
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_1 +
    h_vec.f_1*delta_v_vec.f_11;
  out->f_12 =
    h_vec.f_12*delta_v_vec.f_1 +
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_12;
  out->f_22 =
    h_vec.f_22 +
    2*h_vec.f_12*delta_v_vec.f_2 +
    h_vec.f_11*delta_v_vec.f_2*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_22;
}

/*------------------------------------------------------------------------------
  Liquid internal energy as a function of pressure and tau, used to solve T(u, p)
------------------------------------------------------------------------------*/
void internal_energy_liquid2(uint comp, double pr, double tau, f22_struct *out){
  f22_struct delta_v_vec, h_vec;
  delta_v_vec = memo2_delta_liquid(comp, pr, tau);
  h_vec = memo2_internal_energy(comp, delta_v_vec.f, tau);
  out->f = h_vec.f;
  out->f_1 = h_vec.f_1*delta_v_vec.f_1;
  out->f_2 = h_vec.f_2 + delta_v_vec.f_2*h_vec.f_1;
  out->f_11 =
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_1 +
    h_vec.f_1*delta_v_vec.f_11;
  out->f_12 =
    h_vec.f_12*delta_v_vec.f_1 +
    h_vec.f_11*delta_v_vec.f_1*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_12;
  out->f_22 =
    h_vec.f_22 +
    2*h_vec.f_12*delta_v_vec.f_2 +
    h_vec.f_11*delta_v_vec.f_2*delta_v_vec.f_2 +
    h_vec.f_1*delta_v_vec.f_22;
}

/*------------------------------------------------------------------------------
  Enthalpy solver functors
------------------------------------------------------------------------------*/

class fthl_functor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _enthalpy;
public:
    fthl_functor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_enthalpy(double t){
      this->_enthalpy = t;
    }
    std::tuple<double, double, double>  operator () (double tau) {
      f22_struct out;
      enthalpy_liquid2(this->_comp, this->_pressure, tau, &out);
      return std::make_tuple(
        out.f - this->_enthalpy,
        out.f_2,
        out.f_22
      );
    }
};

class fthv_functor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _enthalpy;
public:
    fthv_functor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_enthalpy(double t){
      this->_enthalpy = t;
    }
    std::tuple<double, double, double>  operator () (double tau) {
      f22_struct out;
      enthalpy_vapor2(this->_comp, this->_pressure, tau, &out);
      return std::make_tuple(
        out.f - this->_enthalpy,
        out.f_2,
        out.f_22
      );
    }
};

/*------------------------------------------------------------------------------
  Entropy solver functors
------------------------------------------------------------------------------*/

class ftsl_functor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _entropy;
public:
    ftsl_functor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_entropy(double t){
      this->_entropy = t;
    }
    std::tuple<double, double, double>  operator () (double tau) {
      f22_struct out;
      entropy_liquid2(this->_comp, this->_pressure, tau, &out);
      return std::make_tuple(
        out.f - this->_entropy,
        out.f_2,
        out.f_22
      );
    }
};

class ftsv_functor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _entropy;
public:
    ftsv_functor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_entropy(double t){
      this->_entropy = t;
    }
    std::tuple<double, double, double>  operator () (double tau) {
      f22_struct out;
      entropy_vapor2(this->_comp, this->_pressure, tau, &out);
      return std::make_tuple(
        out.f - this->_entropy,
        out.f_2,
        out.f_22
      );
    }
};


/*------------------------------------------------------------------------------
  Internal Energy solver functors
------------------------------------------------------------------------------*/

class ftul_functor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _intenergy;
public:
    ftul_functor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_intenergy(double t){
      this->_intenergy = t;
    }
    std::tuple<double, double, double>  operator () (double tau) {
      f22_struct out;
      internal_energy_liquid2(this->_comp, this->_pressure, tau, &out);
      return std::make_tuple(
        out.f - this->_intenergy,
        out.f_2,
        out.f_22
      );
    }
};

class ftuv_functor_deriv
{
private:
    uint _comp;
    double _pressure;
    double _intenergy;
public:
    ftuv_functor_deriv(uint comp){
      this->_comp = comp;
    }
    void set_pressure(double p){
      this->_pressure = p;
    }
    void set_intenergy(double t){
      this->_intenergy = t;
    }
    std::tuple<double, double, double>  operator () (double tau) {
      f22_struct out;
      internal_energy_vapor2(this->_comp, this->_pressure, tau, &out);
      return std::make_tuple(
        out.f - this->_intenergy,
        out.f_2,
        out.f_22
      );
    }
};

/*------------------------------------------------------------------------------
  tau(h, p) solver with derivatives

  This function first checks that if at the given pressure, the enthalpy is
  in the 2-phase region.  If it is, tau is just tau_sat.  If not it tries to
  classify the region in either liquid, vapor, vapor below the tripple point,
  or supercritical. 
------------------------------------------------------------------------------*/

void tau_hp2(uint comp, double ht, double pr, f22_struct *out){
  f12_struct taus_vec;
  double taus, hvs, hls, tau_hi, tau_lo, tau_fail, hc, tau=0;

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=40;
  int digits = std::numeric_limits<double>::digits - 5;
  parameters_struct *dat = &cdata[comp];
  fthv_functor_deriv fghv = fthv_functor_deriv(comp); fghv.set_pressure(pr); fghv.set_enthalpy(ht); 
  fthl_functor_deriv fghl = fthl_functor_deriv(comp); fghl.set_pressure(pr); fghl.set_enthalpy(ht);
  bool is_vapor = 0;

  taus_vec = sat_tau(comp, pr);
  taus = taus_vec.f;
  if(pr < dat->Pc){ // Could be two phase
    hvs = enthalpy(comp, sat_delta_v(comp, taus).f, taus);
    hls = enthalpy(comp, sat_delta_l(comp, taus).f, taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->f = taus;
      out->f_1 = 0;
      out->f_2 = taus_vec.f_1;
      out->f_11 = 0;
      out->f_12 = 0;
      out->f_22 = taus_vec.f_11;
      return;
    }
  }
  f22_struct hvec, outc;
  if(pr >= dat->Pc){ // liquid or supercritical
    enthalpy_liquid2(comp, pr, tau_c(comp), &outc);
    hc = outc.f;
    if(ht < hc){ // enthalpy < critical enthalpy, liquid (Tmin <= T <= Tc)
      tau_lo = tau_c(comp) * 0.999; 
      tau_hi = dat->T_star/dat->T_min; 
      tau_fail = tau_lo;
      is_vapor = 0;
    }
    else{ // enthalpy >= critival enthalpy, supercritical (Tc <= T <= Tmax)
      tau_hi = tau_c(comp) * 1.001;
      tau_lo = dat->T_star/dat->T_max;
      tau_fail = tau_hi;
      is_vapor = 0;
    }
  } 
  else if(pr < dat->Pt){ // it's vapor (unless it's ice)
    tau_lo = dat->T_star/dat->T_max;
    tau_hi = dat->T_star/dat->Tt;
    tau_fail = tau_hi;
    is_vapor = 1;
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    tau_lo = taus;
    tau_hi = dat->T_star/dat->T_min;  //melting_tau_func(pr);
    tau_fail = tau_hi;
    is_vapor = 0;
  }
  else{ // vapor for sure
    tau_lo = dat->T_star/dat->T_max;
    tau_hi = taus;
    tau_fail = tau_hi;
    is_vapor = 1;
  }

  try{
    if(is_vapor) tau = halley_iterate(fghv, (tau_lo + tau_hi)/2.0, tau_lo, tau_hi, digits, h_it_max);
    else tau = halley_iterate(fghl, (tau_lo + tau_hi)/2.0, tau_lo, tau_hi, digits, h_it_max);
  }
  catch(...){
    tau = tau_fail;
    std::cout << "tau(h = " << ht <<  " kJ/kg, p = " << pr << " kPa) solve failed" << std::endl;
  }
  if(isnan(tau)){
    tau = tau_fail;
    std::cout << "tau(h = " << ht <<  " kJ/kg, p = " << pr << " kPa) solve failed" << std::endl;
  }
  if(is_vapor) hvec = memo2_enthalpy_vapor(comp, pr, tau);
  else hvec = memo2_enthalpy_liquid(comp, pr, tau);

  out->f = tau;
  out->f_1 = 1.0/hvec.f_2; // d/dh
  out->f_2 = -out->f_1 * hvec.f_1; // d/dp, triple product
  out->f_11 = -out->f_1 * out->f_1 * out->f_1 * hvec.f_22;
  out->f_12 = -out->f_1 * out->f_1 *
    (hvec.f_12 + hvec.f_22*out->f_2);
  out->f_22 = -out->f_12*hvec.f_1 -
    out->f_1*(hvec.f_11 + hvec.f_12*out->f_2);
}


/*------------------------------------------------------------------------------
  tau(s, p) solver with derivatives

  This function first checks that if at the given pressure, the entropy is
  in the 2-phase region.  If it is, tau is just tau_sat.  If not it tries to
  classify the region in either liquid, vapor, vapor below the tripple point,
  or supercritical. 
------------------------------------------------------------------------------*/

void tau_sp2(uint comp, double ht, double pr, f22_struct *out){
  f12_struct taus_vec;
  double taus, hvs, hls, tau_hi, tau_lo, tau_fail, hc, tau=0;

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=40;
  int digits = std::numeric_limits<double>::digits - 5;
  parameters_struct *dat = &cdata[comp];
  ftsv_functor_deriv fghv = ftsv_functor_deriv(comp); fghv.set_pressure(pr); fghv.set_entropy(ht); 
  ftsl_functor_deriv fghl = ftsl_functor_deriv(comp); fghl.set_pressure(pr); fghl.set_entropy(ht);
  bool is_vapor = 0;

  taus_vec = sat_tau(comp, pr);
  taus = taus_vec.f;
  if(pr < dat->Pc){ // Could be two phase
    hvs = entropy(comp, sat_delta_v(comp, taus).f, taus);
    hls = entropy(comp, sat_delta_l(comp, taus).f, taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->f = taus;
      out->f_1 = 0;
      out->f_2 = taus_vec.f_1;
      out->f_11 = 0;
      out->f_12 = 0;
      out->f_22 = taus_vec.f_11;
      return;
    }
  }
  f22_struct hvec, outc;
  if(pr >= dat->Pc){ // liquid or supercritical
   entropy_liquid2(comp, pr, tau_c(comp), &outc);
    hc = outc.f;
    if(ht < hc){ // entropy < critical entropy, liquid (Tmin <= T <= Tc)
      tau_lo = tau_c(comp) * 0.999; 
      tau_hi = dat->T_star/dat->T_min; 
      tau_fail = tau_lo;
      is_vapor = 0;
    }
    else{ // entropy >= critival eentropy, supercritical (Tc <= T <= Tmax)
      tau_hi = tau_c(comp) * 1.001;
      tau_lo = dat->T_star/dat->T_max;
      tau_fail = tau_hi;
      is_vapor = 0;
    }
  } 
  else if(pr < dat->Pt){ // it's vapor (unless it's ice)
    tau_lo = dat->T_star/dat->T_max;
    tau_hi = dat->T_star/dat->Tt;
    tau_fail = tau_hi;
    is_vapor = 1;
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    tau_lo = taus;
    tau_hi = dat->T_star/dat->T_min;  //melting_tau_func(pr);
    tau_fail = tau_hi;
    is_vapor = 0;
  }
  else{ // vapor for sure
    tau_lo = dat->T_star/dat->T_max;
    tau_hi = taus;
    tau_fail = tau_hi;
    is_vapor = 1;
  }

  try{
    if(is_vapor) tau = halley_iterate(fghv, (tau_lo + tau_hi)/2.0, tau_lo, tau_hi, digits, h_it_max);
    else tau = halley_iterate(fghl, (tau_lo + tau_hi)/2.0, tau_lo, tau_hi, digits, h_it_max);
  }
  catch(...){
    tau = tau_fail;
    std::cout << "tau(s = " << ht <<  " kJ/kg/K, p = " << pr << " kPa) solve failed" << std::endl;
  }
  if(isnan(tau)){
    tau = tau_fail;
    std::cout << "tau(s = " << ht <<  " kJ/kg/K, p = " << pr << " kPa) solve failed" << std::endl;
  }
  if(is_vapor) hvec = memo2_entropy_vapor(comp, pr, tau);
  else hvec = memo2_entropy_liquid(comp, pr, tau);

  out->f = tau;
  out->f_1 = 1.0/hvec.f_2; // d/dh
  out->f_2 = -out->f_1 * hvec.f_1; // d/dp, triple product
  out->f_11 = -out->f_1 * out->f_1 * out->f_1 * hvec.f_22;
  out->f_12 = -out->f_1 * out->f_1 *
    (hvec.f_12 + hvec.f_22*out->f_2);
  out->f_22 = -out->f_12*hvec.f_1 -
    out->f_1*(hvec.f_11 + hvec.f_12*out->f_2);
}

/*------------------------------------------------------------------------------
  tau(u, p) solver with derivatives

  This function first checks that if at the given pressure, the int. energy is
  in the 2-phase region.  If it is, tau is just tau_sat.  If not it tries to
  classify the region in either liquid, vapor, vapor below the tripple point,
  or supercritical. 
------------------------------------------------------------------------------*/

void tau_up2(uint comp, double ht, double pr, f22_struct *out){
  f12_struct taus_vec;
  double taus, hvs, hls, tau_hi, tau_lo, tau_fail, hc, tau=0;

  using namespace boost::math::tools;
  std::uintmax_t h_it_max=40;
  int digits = std::numeric_limits<double>::digits - 5;
  parameters_struct *dat = &cdata[comp];
  ftuv_functor_deriv fghv = ftuv_functor_deriv(comp); fghv.set_pressure(pr); fghv.set_intenergy(ht); 
  ftul_functor_deriv fghl = ftul_functor_deriv(comp); fghl.set_pressure(pr); fghl.set_intenergy(ht);
  bool is_vapor = 0;

  taus_vec = sat_tau(comp, pr);
  taus = taus_vec.f;
  if(pr < dat->Pc){ // Could be two phase
    hvs = internal_energy(comp, sat_delta_v(comp, taus).f, taus);
    hls = internal_energy(comp, sat_delta_l(comp, taus).f, taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->f = taus;
      out->f_1 = 0;
      out->f_2 = taus_vec.f_1;
      out->f_11 = 0;
      out->f_12 = 0;
      out->f_22 = taus_vec.f_11;
      return;
    }
  }
  f22_struct hvec, outc;
  if(pr >= dat->Pc){ // liquid or supercritical
    internal_energy_liquid2(comp, pr, tau_c(comp), &outc);
    hc = outc.f;
    if(ht < hc){ // enthalpy < critical enthalpy, liquid (Tmin <= T <= Tc)
      tau_lo = tau_c(comp) * 0.999; 
      tau_hi = dat->T_star/dat->T_min; 
      tau_fail = tau_lo;
      is_vapor = 0;
    }
    else{ // enthalpy >= critival enthalpy, supercritical (Tc <= T <= Tmax)
      tau_hi = tau_c(comp) * 1.001;
      tau_lo = dat->T_star/dat->T_max;
      tau_fail = tau_hi;
      is_vapor = 0;
    }
  } 
  else if(pr < dat->Pt){ // it's vapor (unless it's ice)
    tau_lo = dat->T_star/dat->T_max;
    tau_hi = dat->T_star/dat->Tt;
    tau_fail = tau_hi;
    is_vapor = 1;
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    tau_lo = taus;
    tau_hi = dat->T_star/dat->T_min;  //melting_tau_func(pr);
    tau_fail = tau_hi;
    is_vapor = 0;
  }
  else{ // vapor for sure
    tau_lo = dat->T_star/dat->T_max;
    tau_hi = taus;
    tau_fail = tau_hi;
    is_vapor = 1;
  }

  try{
    if(is_vapor) tau = halley_iterate(fghv, (tau_lo + tau_hi)/2.0, tau_lo, tau_hi, digits, h_it_max);
    else tau = halley_iterate(fghl, (tau_lo + tau_hi)/2.0, tau_lo, tau_hi, digits, h_it_max);
  }
  catch(...){
    tau = tau_fail;
    std::cout << "tau(u = " << ht <<  " kJ/kg, p = " << pr << " kPa) solve failed" << std::endl;
  }
  if(isnan(tau)){
    tau = tau_fail;
    std::cout << "tau(u = " << ht <<  " kJ/kg, p = " << pr << " kPa) solve failed" << std::endl;
  }
  if(is_vapor) hvec = memo2_internal_energy_vapor(comp, pr, tau);
  else hvec = memo2_internal_energy_liquid(comp, pr, tau);

  out->f = tau;
  out->f_1 = 1.0/hvec.f_2; // d/dh
  out->f_2 = -out->f_1 * hvec.f_1; // d/dp, triple product
  out->f_11 = -out->f_1 * out->f_1 * out->f_1 * hvec.f_22;
  out->f_12 = -out->f_1 * out->f_1 *
    (hvec.f_12 + hvec.f_22*out->f_2);
  out->f_22 = -out->f_12*hvec.f_1 -
    out->f_1*(hvec.f_11 + hvec.f_12*out->f_2);
}


/*------------------------------------------------------------------------------
  Vapor fraction functions
------------------------------------------------------------------------------*/

void vf_hp2(uint comp, double ht, double pr, f22_struct *out){
  parameters_struct *dat = &cdata[comp];

  if(pr >= dat->Pc){ // consider supercritical fluid liquid
    out->f = 0;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }
  f22_struct hsat_l_vec, hsat_v_vec;
  f12_struct taus_vec;

  taus_vec = sat_tau(comp, pr);
  hsat_l_vec = memo2_enthalpy_liquid(comp, pr, taus_vec.f);
  hsat_v_vec = memo2_enthalpy_vapor(comp, pr, taus_vec.f);
  double hv = hsat_v_vec.f;
  double hl = hsat_l_vec.f;

  if(ht < hl){ // Liquid, vapor fraction is 0 and no change with T, P
    out->f = 0;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  if(ht > hv){ // Vapor, vapor fraction is 1 and no change with T, P
    out->f = 1;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  // if you're here then your in the saturated area
  double vf = (ht - hl)/(hv - hl);
  out->f = vf;

  double dhvdp = hsat_v_vec.f_1 + hsat_v_vec.f_2*taus_vec.f_1;
  double dhldp = hsat_l_vec.f_1 + hsat_l_vec.f_2*taus_vec.f_1;

  out->f_1 = 1.0/(hv - hl);
  out->f_2 = -dhldp/(hv - hl) - vf/(hv-hl)*(dhvdp - dhldp);

  double d2hvdp2 = hsat_v_vec.f_11 +
    2*hsat_v_vec.f_12*taus_vec.f_1 +
    hsat_v_vec.f_22*taus_vec.f_1*taus_vec.f_1 +
    hsat_v_vec.f_2*taus_vec.f_11;
  double d2hldp2 = hsat_l_vec.f_11 +
    2*hsat_l_vec.f_12*taus_vec.f_1 +
    hsat_l_vec.f_22*taus_vec.f_1*taus_vec.f_1 +
    hsat_l_vec.f_2*taus_vec.f_11;
  out->f_11 = 0;
  out->f_12 = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
  out->f_22 = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
            2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
            (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
}

void vf_sp2(uint comp, double ht, double pr, f22_struct *out){
  parameters_struct *dat = &cdata[comp];

  if(pr >= dat->Pc){ // consider supercritical fluid liquid
    out->f = 0;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }
  f22_struct hsat_l_vec, hsat_v_vec;
  f12_struct taus_vec;

  taus_vec = sat_tau(comp, pr);
  hsat_l_vec = memo2_entropy_liquid(comp, pr, taus_vec.f);
  hsat_v_vec = memo2_entropy_vapor(comp, pr, taus_vec.f);
  double hv = hsat_v_vec.f;
  double hl = hsat_l_vec.f;

  if(ht < hl){
    out->f = 0;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  if(ht > hv){
    out->f = 1;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  // if you're here then your in the saturated area
  double vf = (ht - hl)/(hv - hl);
  out->f = vf;

  double dhvdp = hsat_v_vec.f_1 + hsat_v_vec.f_2*taus_vec.f_1;
  double dhldp = hsat_l_vec.f_1 + hsat_l_vec.f_2*taus_vec.f_1;

  out->f_1 = 1.0/(hv - hl);
  out->f_2 = -dhldp/(hv - hl) - vf/(hv-hl)*(dhvdp - dhldp);

  double d2hvdp2 = hsat_v_vec.f_11 +
    2*hsat_v_vec.f_12*taus_vec.f_1 +
    hsat_v_vec.f_22*taus_vec.f_1*taus_vec.f_1 +
    hsat_v_vec.f_2*taus_vec.f_11;
  double d2hldp2 = hsat_l_vec.f_11 +
    2*hsat_l_vec.f_12*taus_vec.f_1 +
    hsat_l_vec.f_22*taus_vec.f_1*taus_vec.f_1 +
    hsat_l_vec.f_2*taus_vec.f_11;
  out->f_11 = 0;
  out->f_12 = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
  out->f_22 = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
            2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
            (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
}

void vf_up2(uint comp, double ht, double pr, f22_struct *out){
  parameters_struct *dat = &cdata[comp];

  if(pr >= dat->Pc){ // consider supercritical fluid liquid
    out->f = 0;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  f22_struct hsat_l_vec, hsat_v_vec;
  f12_struct taus_vec;

  taus_vec = sat_tau(comp, pr);
  hsat_l_vec = memo2_internal_energy_liquid(comp, pr, taus_vec.f);
  hsat_v_vec = memo2_internal_energy_vapor(comp, pr, taus_vec.f);
  double hv = hsat_v_vec.f;
  double hl = hsat_l_vec.f;

  if(ht < hl){
    out->f = 0;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  if(ht > hv){
    out->f = 1;
    out->f_1 = 0;
    out->f_2 = 0;
    out->f_11 = 0;
    out->f_12 = 0;
    out->f_22 = 0;
    return;
  }

  // if you're here then your in the saturated area
  double vf = (ht - hl)/(hv - hl);
  out->f = vf;

  double dhvdp = hsat_v_vec.f_1 + hsat_v_vec.f_2*taus_vec.f_1;
  double dhldp = hsat_l_vec.f_1 + hsat_l_vec.f_2*taus_vec.f_1;

  out->f_1 = 1.0/(hv - hl);
  out->f_2 = -dhldp/(hv - hl) - vf/(hv-hl)*(dhvdp - dhldp);

  double d2hvdp2 = hsat_v_vec.f_11 +
    2*hsat_v_vec.f_12*taus_vec.f_1 +
    hsat_v_vec.f_22*taus_vec.f_1*taus_vec.f_1 +
    hsat_v_vec.f_2*taus_vec.f_11;
  double d2hldp2 = hsat_l_vec.f_11 +
    2*hsat_l_vec.f_12*taus_vec.f_1 +
    hsat_l_vec.f_22*taus_vec.f_1*taus_vec.f_1 +
    hsat_l_vec.f_2*taus_vec.f_11;
  out->f_11 = 0;
  out->f_12 = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
  out->f_22 = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
            2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
            (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
}

/*------------------------------------------------------------------------------
  Memo functions
------------------------------------------------------------------------------*/

MEMO2_FUNCTION(memo2_enthalpy_vapor, enthalpy_vapor2, memo_table_enthalpy_vapor2)
MEMO2_FUNCTION(memo2_entropy_vapor, entropy_vapor2, memo_table_entropy_vapor2)
MEMO2_FUNCTION(memo2_internal_energy_vapor, internal_energy_vapor2, memo_table_internal_energy_vapor2)
MEMO2_FUNCTION(memo2_enthalpy_liquid, enthalpy_liquid2, memo_table_enthalpy_liquid2)
MEMO2_FUNCTION(memo2_entropy_liquid, entropy_liquid2, memo_table_entropy_liquid2)
MEMO2_FUNCTION(memo2_internal_energy_liquid, internal_energy_liquid2, memo_table_internal_energy_liquid2)

MEMO2_FUNCTION(memo2_tau_hp, tau_hp2, memo_table_tau_hp2)
MEMO2_FUNCTION(memo2_tau_sp, tau_sp2, memo_table_tau_sp2)
MEMO2_FUNCTION(memo2_tau_up, tau_up2, memo_table_tau_up2)

MEMO2_FUNCTION(memo2_vf_hp, vf_hp2, memo_table_vf_hp2)
MEMO2_FUNCTION(memo2_vf_sp, vf_sp2, memo_table_vf_sp2)
MEMO2_FUNCTION(memo2_vf_up, vf_up2, memo_table_vf_up2)
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

#include<unordered_map>
#include<boost/functional/hash.hpp>
#include "state.h"
#include "props.h"
#include "solver.h"
#include "delta.h"

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_enthalpy_vapor2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_entropy_vapor2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_internal_energy_vapor2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_enthalpy_liquid2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_entropy_liquid2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_internal_energy_liquid2;

static const uint f_d = (uint)deriv2_enum::f_d;
static const uint f_t = (uint)deriv2_enum::f_t;
static const uint f_dd = (uint)deriv2_enum::f_dd;
static const uint f_dt = (uint)deriv2_enum::f_dt;
static const uint f_tt = (uint)deriv2_enum::f_tt;

static const uint f_p = (uint)deriv2_enum::f_d;
static const uint f_pp = (uint)deriv2_enum::f_dd;
static const uint f_pt = (uint)deriv2_enum::f_dt;

void enthalpy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_vapor(comp, pr, tau);
  h_vec = memo2_enthalpy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f_p) = h_vec->at(f_d)*delta_v_vec->at(f_p);
  out->at(f_t) = h_vec->at(f_t) + delta_v_vec->at(f_t)*h_vec->at(f_d);
  out->at(f_pp) =
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_p) +
    h_vec->at(f_d)*delta_v_vec->at(f_pp);
  out->at(f_pt) =
    h_vec->at(f_dt)*delta_v_vec->at(f_p) +
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_dt);
  out->at(f_tt) =
    h_vec->at(f_tt) +
    2*h_vec->at(f_dt)*delta_v_vec->at(f_t) +
    h_vec->at(f_dd)*delta_v_vec->at(f_t)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_tt);
}

std::vector<double> *memo2_enthalpy_vapor(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_enthalpy_vapor2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_enthalpy_vapor2.size() > MAX_MEMO_PROP) memo_table_enthalpy_vapor2.clear();
  yvec_ptr = &memo_table_enthalpy_vapor2[std::make_tuple(comp, pr, tau)];
  enthalpy_vapor2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

void entropy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_vapor(comp, pr, tau);
  h_vec = memo2_entropy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f_p) = h_vec->at(f_d)*delta_v_vec->at(f_p);
  out->at(f_t) = h_vec->at(f_t) + delta_v_vec->at(f_t)*h_vec->at(f_d);
  out->at(f_pp) =
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_p) +
    h_vec->at(f_d)*delta_v_vec->at(f_pp);
  out->at(f_pt) =
    h_vec->at(f_dt)*delta_v_vec->at(f_p) +
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_dt);
  out->at(f_tt) =
    h_vec->at(f_tt) +
    2*h_vec->at(f_dt)*delta_v_vec->at(f_t) +
    h_vec->at(f_dd)*delta_v_vec->at(f_t)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_tt);
}

std::vector<double> *memo2_entropy_vapor(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_entropy_vapor2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_entropy_vapor2.size() > MAX_MEMO_PROP) memo_table_entropy_vapor2.clear();
  yvec_ptr = &memo_table_entropy_vapor2[std::make_tuple(comp, pr, tau)];
  entropy_vapor2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

void internal_energy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_vapor(comp, pr, tau);
  h_vec = memo2_internal_energy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f_p) = h_vec->at(f_d)*delta_v_vec->at(f_p);
  out->at(f_t) = h_vec->at(f_t) + delta_v_vec->at(f_t)*h_vec->at(f_d);
  out->at(f_pp) =
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_p) +
    h_vec->at(f_d)*delta_v_vec->at(f_pp);
  out->at(f_pt) =
    h_vec->at(f_dt)*delta_v_vec->at(f_p) +
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_dt);
  out->at(f_tt) =
    h_vec->at(f_tt) +
    2*h_vec->at(f_dt)*delta_v_vec->at(f_t) +
    h_vec->at(f_dd)*delta_v_vec->at(f_t)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_tt);
}

std::vector<double> *memo2_internal_energy_vapor(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_internal_energy_vapor2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_internal_energy_vapor2.size() > MAX_MEMO_PROP) memo_table_internal_energy_vapor2.clear();
  yvec_ptr = &memo_table_internal_energy_vapor2[std::make_tuple(comp, pr, tau)];
  internal_energy_vapor2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

void enthalpy_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_liquid(comp, pr, tau);
  h_vec = memo2_enthalpy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f_p) = h_vec->at(f_d)*delta_v_vec->at(f_p);
  out->at(f_t) = h_vec->at(f_t) + delta_v_vec->at(f_t)*h_vec->at(f_d);
  out->at(f_pp) =
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_p) +
    h_vec->at(f_d)*delta_v_vec->at(f_pp);
  out->at(f_pt) =
    h_vec->at(f_dt)*delta_v_vec->at(f_p) +
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_dt);
  out->at(f_tt) =
    h_vec->at(f_tt) +
    2*h_vec->at(f_dt)*delta_v_vec->at(f_t) +
    h_vec->at(f_dd)*delta_v_vec->at(f_t)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_tt);
}

std::vector<double> *memo2_enthalpy_liquid(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_enthalpy_liquid2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_enthalpy_liquid2.size() > MAX_MEMO_PROP) memo_table_enthalpy_liquid2.clear();
  yvec_ptr = &memo_table_enthalpy_liquid2[std::make_tuple(comp, pr, tau)];
  enthalpy_liquid2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

void entropy_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_liquid(comp, pr, tau);
  h_vec = memo2_entropy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f_p) = h_vec->at(f_d)*delta_v_vec->at(f_p);
  out->at(f_t) = h_vec->at(f_t) + delta_v_vec->at(f_t)*h_vec->at(f_d);
  out->at(f_pp) =
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_p) +
    h_vec->at(f_d)*delta_v_vec->at(f_pp);
  out->at(f_pt) =
    h_vec->at(f_dt)*delta_v_vec->at(f_p) +
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_dt);
  out->at(f_tt) =
    h_vec->at(f_tt) +
    2*h_vec->at(f_dt)*delta_v_vec->at(f_t) +
    h_vec->at(f_dd)*delta_v_vec->at(f_t)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_tt);
}

std::vector<double> *memo2_entropy_liquid(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_entropy_liquid2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_entropy_liquid2.size() > MAX_MEMO_PROP) memo_table_entropy_liquid2.clear();
  yvec_ptr = &memo_table_entropy_liquid2[std::make_tuple(comp, pr, tau)];
  entropy_liquid2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

void internal_energy_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_liquid(comp, pr, tau);
  h_vec = memo2_internal_energy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f_p) = h_vec->at(f_d)*delta_v_vec->at(f_p);
  out->at(f_t) = h_vec->at(f_t) + delta_v_vec->at(f_t)*h_vec->at(f_d);
  out->at(f_pp) =
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_p) +
    h_vec->at(f_d)*delta_v_vec->at(f_pp);
  out->at(f_pt) =
    h_vec->at(f_dt)*delta_v_vec->at(f_p) +
    h_vec->at(f_dd)*delta_v_vec->at(f_p)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_dt);
  out->at(f_tt) =
    h_vec->at(f_tt) +
    2*h_vec->at(f_dt)*delta_v_vec->at(f_t) +
    h_vec->at(f_dd)*delta_v_vec->at(f_t)*delta_v_vec->at(f_t) +
    h_vec->at(f_d)*delta_v_vec->at(f_tt);
}

std::vector<double> *memo2_internal_energy_liquid(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_internal_energy_liquid2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_internal_energy_liquid2.size() > MAX_MEMO_PROP) memo_table_internal_energy_liquid2.clear();
  yvec_ptr = &memo_table_internal_energy_liquid2[std::make_tuple(comp, pr, tau)];
  internal_energy_liquid2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

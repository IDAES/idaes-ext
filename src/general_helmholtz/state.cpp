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
#include<iostream>
#include "function_pointers.h"
#include "state.h"
#include "props.h"
#include "param.h"
#include "solver.h"
#include "delta.h"
#include "sat.h"

struct state_solve_dat{ // Solver wrapper data structure for const args
  comp_enum comp;
  double p;
  double h;
};

/*------------------------------------------------------------------------------
  Memoization tables
------------------------------------------------------------------------------*/

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

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_tau_hp2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_tau_sp2;

static const uint f_d = (uint)deriv2_enum::f_d;
static const uint f_t = (uint)deriv2_enum::f_t;
static const uint f_dd = (uint)deriv2_enum::f_dd;
static const uint f_dt = (uint)deriv2_enum::f_dt;
static const uint f_tt = (uint)deriv2_enum::f_tt;

static const uint f_p = (uint)deriv2_enum::f_d;
static const uint f_pp = (uint)deriv2_enum::f_dd;
static const uint f_pt = (uint)deriv2_enum::f_dt;

/*------------------------------------------------------------------------------
  Vapor enthalpy as a function of pressure and tau, used to solve T(h, p)
------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------
  Vapor entropy as a function of pressure and tau, used to solve T(s, p)
------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------
  Vapor internal energy as a function of pressure and tau, used to solve T(u, p)
------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------
  Vapor enthalpy as a function of pressure and tau, used to solve T(h, p)
------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------
  Liquid entropy as a function of pressure and tau, used to solve T(s, p)
------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------
  Liquid internal energy as a function of pressure and tau, used to solve T(u, p)
------------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------------
  Enthalpy solver function wrappers
------------------------------------------------------------------------------*/

double f_thv(double tau, void* dat){
  // This is for bracket solver, so may want to use an underlying function that
  // doesn't bother calculating deriavtives, but I don't currently have such a
  // version.
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> out;
  enthalpy_vapor2(sd->comp, sd->p, tau, &out);
  return out[0] - sd->h;
}

double f_thl(double tau, void* dat){
  // This is for bracket solver, so may want to use an underlying function that
  // doesn't bother calculating deriavtives, but I don't currently have such a
  // version.
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> out;
  enthalpy_liquid2(sd->comp, sd->p, tau, &out);
  return out[0] - sd->h;
}

void f_thv2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  enthalpy_vapor2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f_t]; // here p is constant
  out->at(2) = hvec[f_tt]; // here p is constant
}

void f_thl2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  enthalpy_liquid2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f_t]; // here p is constant
  out->at(2) = hvec[f_tt]; // here p is constant
}

/*------------------------------------------------------------------------------
  tau(h, p) solver with derivatives

  This function first checks that if at the given pressure, the enthalpy is
  in the 2-phase region.  If it is, tau is just tau_sat.  If not it tries to
  classify the region in either liquid, vapor, vapor below the tripple point,
  or supercritical. Once the region is itentified, a bracketing method is used
  to narrow down the initial guess before solving with a Newton method.
------------------------------------------------------------------------------*/

void tau_hp2(comp_enum comp, double ht, double pr, std::vector<double> *out){
  std::vector<double> *taus_vec_ptr;
  double taus, hvs, hls, tau_hi, tau_lo, tau=0;

  taus_vec_ptr = sat_tau(comp, pr);
  taus = taus_vec_ptr->at(0);
  if(pr > Pt[comp] && pr < Pc[comp]){ // Could be two phase
    hvs = enthalpy(comp, sat_delta_v(comp, taus)->at(0), taus);
    hls = enthalpy(comp, sat_delta_l(comp, taus)->at(0), taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->resize(6);
      out->at(0) = taus;
      out->at(f_p) = taus_vec_ptr->at(1);
      out->at(f_t) = 0;
      out->at(f_pp) = taus_vec_ptr->at(2);
      out->at(f_pt) = 0;
      out->at(f_tt) = 0;
      return;
    }
  }
  std::vector<double> *hvec_ptr;
  state_solve_dat sd;
  std::vector<double> hout;
  sd.p = pr;
  sd.h = ht;
  sd.comp = comp;
  if(pr >= Pc[comp]){ // it's liquid or supercritical
    // Assume between melting and Tmax
    //   This is a pretty big temperature range bracketing for a good guess for
    //   Newton method may be slow here.  May want some aux functions to break
    //   it up a bit more
    tau_lo = Tc[comp]/T_max[comp];
    tau_hi = Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_thl, tau_lo, tau_hi, &tau, 20, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_liquid(comp, pr, tau);
  }
  else if(pr < Pt[comp]){ // it's vapor (unless it's ice)
    // Assume between sublimation temperature and Tmax
    tau_lo = Tc[comp]/T_max[comp];
    tau_hi = Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_thv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_vapor(comp, pr, tau);
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    // Assume between melting temperature and sat temperature
    tau_lo = taus;
    tau_hi = Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_thl, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_liquid(comp, pr, tau);
  }
  else{ //if(ht >= hvs){ // vapor for sure
    // Assume between saturation temperature and Tmax
    tau_lo = Tc[comp]/T_max[comp];
    tau_hi = taus;
    int n1 = bracket(f_thv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_vapor(comp, pr, tau);
  }
  out->resize(6);
  out->at(0) = tau;
  // here, for indexing of **out only**:
  //   d will be derivative with resepect to h and
  //   t will be with respect to p
  out->at(f_d) = 1.0/hvec_ptr->at(f_t);           // d/dh
  out->at(f_t) = -out->at(f_d) * hvec_ptr->at(f_p); // d/dp, tripple product
  out->at(f_dd) = -out->at(f_d) * out->at(f_d) * out->at(f_d) * hvec_ptr->at(f_tt);
  out->at(f_dt) = -out->at(f_d) * out->at(f_d) *
    (hvec_ptr->at(f_pt) + hvec_ptr->at(f_tt)*out->at(f_t));
  out->at(f_tt) = -out->at(f_dt)*hvec_ptr->at(f_p) -
    out->at(f_d)*(hvec_ptr->at(f_pp) + hvec_ptr->at(f_pt)*out->at(f_t));
}

std::vector<double> *memo2_tau_hp(comp_enum comp, double ht, double pr){
  try{
    return &memo_table_tau_hp2.at(std::make_tuple(comp, ht, pr));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_tau_hp2.size() > MAX_MEMO_PROP) memo_table_tau_hp2.clear();
  yvec_ptr = &memo_table_tau_hp2[std::make_tuple(comp, ht, pr)];
  tau_hp2(comp, ht, pr, yvec_ptr);
  return yvec_ptr;
}

/*------------------------------------------------------------------------------
  Entropy solver function wrappers
------------------------------------------------------------------------------*/

double f_tsv(double tau, void* dat){
  // This is for bracket solver, so may want to use an underlying function that
  // doesn't bother calculating deriavtives, but I don't currently have such a
  // version.
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> out;
  entropy_vapor2(sd->comp, sd->p, tau, &out);
  return out[0] - sd->h;
}

double f_tsl(double tau, void* dat){
  // This is for bracket solver, so may want to use an underlying function that
  // doesn't bother calculating deriavtives, but I don't currently have such a
  // version.
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> out;
  entropy_liquid2(sd->comp, sd->p, tau, &out);
  return out[0] - sd->h;
}

void f_tsv2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  entropy_vapor2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f_t]; // here p is constant
  out->at(2) = hvec[f_tt]; // here p is constant
}

void f_tsl2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  entropy_liquid2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f_t]; // here p is constant
  out->at(2) = hvec[f_tt]; // here p is constant
}

/*------------------------------------------------------------------------------
  tau(s, p) solver with derivatives

  This function first checks that if at the given pressure, the enthalpy is
  in the 2-phase region.  If it is, tau is just tau_sat.  If not it tries to
  classify the region in either liquid, vapor, vapor below the tripple point,
  or supercritical. Once the region is itentified, a bracketing method is used
  to narrow down the initial guess before solving with a Newton method.
------------------------------------------------------------------------------*/

void tau_sp2(comp_enum comp, double ht, double pr, std::vector<double> *out){
  std::vector<double> *taus_vec_ptr;
  double taus, hvs, hls, tau_hi, tau_lo, tau=0;

  taus_vec_ptr = sat_tau(comp, pr);
  taus = taus_vec_ptr->at(0);
  if(pr > Pt[comp] && pr < Pc[comp]){ // Could be two phase
    hvs = entropy(comp, sat_delta_v(comp, taus)->at(0), taus);
    hls = entropy(comp, sat_delta_l(comp, taus)->at(0), taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->resize(6);
      out->at(0) = taus;
      out->at(f_p) = taus_vec_ptr->at(1);
      out->at(f_t) = 0;
      out->at(f_pp) = taus_vec_ptr->at(2);
      out->at(f_pt) = 0;
      out->at(f_tt) = 0;
      return;
    }
  }
  std::vector<double> *hvec_ptr;
  state_solve_dat sd;
  std::vector<double> hout;
  sd.p = pr;
  sd.h = ht;
  sd.comp = comp;
  if(pr >= Pc[comp]){ // it's liquid or supercritical
    // Assume between melting and Tmax
    //   This is a pretty big temperature range bracketing for a good guess for
    //   Newton method may be slow here.  May want some aux functions to break
    //   it up a bit more
    tau_lo = Tc[comp]/T_max[comp];
    tau_hi = Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tsl, tau_lo, tau_hi, &tau, 20, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tsl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_liquid(comp, pr, tau);
  }
  else if(pr < Pt[comp]){ // it's vapor (unless it's ice)
    // Assume between sublimation temperature and Tmax
    tau_lo = Tc[comp]/T_max[comp];
    tau_hi = Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tsv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tsv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_vapor(comp, pr, tau);
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    // Assume between melting temperature and sat temperature
    std::cout << "liquid" << std::endl;
    tau_lo = taus;
    tau_hi = Tc[comp]/melting_temperature_func[comp](pr);

    std::cout << "detla_liquid(" << Tc[comp]/tau_hi << ") = " << delta_liquid(comp, pr, tau_hi) << std::endl;
    std::cout << "detla_liquid(" << Tc[comp]/tau_lo << ") = " << delta_liquid(comp, pr, tau_lo) << std::endl;

    int n1 = bracket(f_tsl, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    std::cout << tau << std::endl;
    int n2 = halley(f_tsl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_liquid(comp, pr, tau);
  }
  else{ //if(ht >= hvs){ // vapor for sure
    // Assume between saturation temperature and Tmax
    tau_lo = Tc[comp]/T_max[comp];
    tau_hi = taus;
    int n1 = bracket(f_tsv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tsv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_vapor(comp, pr, tau);
  }
  out->resize(6);
  out->at(0) = tau;
  // here, for indexing of **out only**:
  //   d will be derivative with resepect to h and
  //   t will be with respect to p
  out->at(f_d) = 1.0/hvec_ptr->at(f_t);           // d/dh
  out->at(f_t) = -out->at(f_d) * hvec_ptr->at(f_p); // d/dp, tripple product
  out->at(f_dd) = -out->at(f_d) * out->at(f_d) * out->at(f_d) * hvec_ptr->at(f_tt);
  out->at(f_dt) = -out->at(f_d) * out->at(f_d) *
    (hvec_ptr->at(f_pt) + hvec_ptr->at(f_tt)*out->at(f_t));
  out->at(f_tt) = -out->at(f_dt)*hvec_ptr->at(f_p) -
    out->at(f_d)*(hvec_ptr->at(f_pp) + hvec_ptr->at(f_pt)*out->at(f_t));
}

std::vector<double> *memo2_tau_sp(comp_enum comp, double st, double pr){
  try{
    return &memo_table_tau_sp2.at(std::make_tuple(comp, st, pr));
  }
  catch(std::out_of_range){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_tau_sp2.size() > MAX_MEMO_PROP) memo_table_tau_sp2.clear();
  yvec_ptr = &memo_table_tau_sp2[std::make_tuple(comp, st, pr)];
  tau_sp2(comp, st, pr, yvec_ptr);
  return yvec_ptr;
}

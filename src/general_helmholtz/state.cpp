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

#include "components/function_pointers.h"
#include "state.h"
#include "props.h"
#include "solver.h"
#include "delta.h"
#include "sat.h"
#include <iostream>

struct state_solve_dat{ // Solver wrapper data structure for const args
  comp_enum comp;
  double p;
  double h;
};

/*------------------------------------------------------------------------------
  Memoization tables
------------------------------------------------------------------------------*/

prop_memo_table2 memo_table_enthalpy_vapor2;
prop_memo_table2 memo_table_entropy_vapor2;
prop_memo_table2 memo_table_internal_energy_vapor2;
prop_memo_table2 memo_table_enthalpy_liquid2;
prop_memo_table2 memo_table_entropy_liquid2;
prop_memo_table2 memo_table_internal_energy_liquid2;
prop_memo_table2 memo_table_tau_hp2;
prop_memo_table2 memo_table_tau_sp2;
prop_memo_table2 memo_table_tau_up2;
prop_memo_table2 memo_table_vf_hp2;
prop_memo_table2 memo_table_vf_sp2;
prop_memo_table2 memo_table_vf_up2;

/*------------------------------------------------------------------------------
  Vapor enthalpy as a function of pressure and tau, used to solve T(h, p)
------------------------------------------------------------------------------*/

void enthalpy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  std::vector<double> *delta_v_vec, *h_vec;
  delta_v_vec = memo2_delta_vapor(comp, pr, tau);
  h_vec = memo2_enthalpy(comp, delta_v_vec->at(0), tau);
  out->resize(6);
  out->at(0) = h_vec->at(0);
  out->at(f2_1) = h_vec->at(f2_1)*delta_v_vec->at(f2_1);
  out->at(f2_2) = h_vec->at(f2_2) + delta_v_vec->at(f2_2)*h_vec->at(f2_1);
  out->at(f2_11) =
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_11);
  out->at(f2_12) =
    h_vec->at(f2_12)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_12);
  out->at(f2_22) =
    h_vec->at(f2_22) +
    2*h_vec->at(f2_12)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_2)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_22);
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
  out->at(f2_1) = h_vec->at(f2_1)*delta_v_vec->at(f2_1);
  out->at(f2_2) = h_vec->at(f2_2) + delta_v_vec->at(f2_2)*h_vec->at(f2_1);
  out->at(f2_11) =
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_11);
  out->at(f2_12) =
    h_vec->at(f2_12)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_12);
  out->at(f2_22) =
    h_vec->at(f2_22) +
    2*h_vec->at(f2_12)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_2)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_22);
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
  out->at(f2_1) = h_vec->at(f2_1)*delta_v_vec->at(f2_1);
  out->at(f2_2) = h_vec->at(f2_2) + delta_v_vec->at(f2_2)*h_vec->at(f2_1);
  out->at(f2_11) =
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_11);
  out->at(f2_12) =
    h_vec->at(f2_12)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_12);
  out->at(f2_22) =
    h_vec->at(f2_22) +
    2*h_vec->at(f2_12)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_2)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_22);
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
  out->at(f2_1) = h_vec->at(f2_1)*delta_v_vec->at(f2_1);
  out->at(f2_2) = h_vec->at(f2_2) + delta_v_vec->at(f2_2)*h_vec->at(f2_1);
  out->at(f2_11) =
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_11);
  out->at(f2_12) =
    h_vec->at(f2_12)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_12);
  out->at(f2_22) =
    h_vec->at(f2_22) +
    2*h_vec->at(f2_12)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_2)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_22);
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
  out->at(f2_1) = h_vec->at(f2_1)*delta_v_vec->at(f2_1);
  out->at(f2_2) = h_vec->at(f2_2) + delta_v_vec->at(f2_2)*h_vec->at(f2_1);
  out->at(f2_11) =
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_11);
  out->at(f2_12) =
    h_vec->at(f2_12)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_12);
  out->at(f2_22) =
    h_vec->at(f2_22) +
    2*h_vec->at(f2_12)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_2)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_22);
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
  out->at(f2_1) = h_vec->at(f2_1)*delta_v_vec->at(f2_1);
  out->at(f2_2) = h_vec->at(f2_2) + delta_v_vec->at(f2_2)*h_vec->at(f2_1);
  out->at(f2_11) =
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_11);
  out->at(f2_12) =
    h_vec->at(f2_12)*delta_v_vec->at(f2_1) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_1)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_12);
  out->at(f2_22) =
    h_vec->at(f2_22) +
    2*h_vec->at(f2_12)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_11)*delta_v_vec->at(f2_2)*delta_v_vec->at(f2_2) +
    h_vec->at(f2_1)*delta_v_vec->at(f2_22);
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
  out->at(1) = hvec[f2_2]; // here p is constant
  out->at(2) = hvec[f2_22]; // here p is constant
}

void f_thl2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  enthalpy_liquid2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f2_2]; // here p is constant
  out->at(2) = hvec[f2_22]; // here p is constant
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
  if(pr < param::Pc[comp]){ // Could be two phase
    hvs = enthalpy(comp, sat_delta_v(comp, taus)->at(0), taus);
    hls = enthalpy(comp, sat_delta_l(comp, taus)->at(0), taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->resize(6);
      out->at(0) = taus;
      out->at(f2_1) = 0;
      out->at(f2_2) = taus_vec_ptr->at(1);
      out->at(f2_11) = 0;
      out->at(f2_12) = 0;
      out->at(f2_22) = taus_vec_ptr->at(2);
      return;
    }
  }
  std::vector<double> *hvec_ptr;
  state_solve_dat sd;
  std::vector<double> hout;
  sd.p = pr;
  sd.h = ht;
  sd.comp = comp;
  if(pr >= param::Pc[comp]){ // it's liquid or supercritical
    // Assume between melting and Tmax
    //   This is a pretty big temperature range bracketing for a good guess for
    //   Newton method may be slow here.  May want some aux functions to break
    //   it up a bit more
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = param::Tc[comp]/(melting_temperature_func[comp](pr));
    int n1 = bracket(f_thl, tau_lo, tau_hi, &tau, 1, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_liquid(comp, pr, tau);
  }
  else if(pr < param::Pt[comp]){ // it's vapor (unless it's ice)
    // Assume between sublimation temperature and Tmax
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_thv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_vapor(comp, pr, tau);
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    // Assume between melting temperature and sat temperature
    tau_lo = taus;
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_thl, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_thl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_enthalpy_liquid(comp, pr, tau);
  }
  else{ //if(ht >= hvs){ // vapor for sure
    // Assume between saturation temperature and Tmax
    tau_lo = param::Tc[comp]/param::T_max[comp];
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
  out->at(f2_1) = 1.0/hvec_ptr->at(f2_2);           // d/dh
  out->at(f2_2) = -out->at(f2_1) * hvec_ptr->at(f2_1); // d/dp, tripple product
  out->at(f2_11) = -out->at(f2_1) * out->at(f2_1) * out->at(f2_1) * hvec_ptr->at(f2_22);
  out->at(f2_12) = -out->at(f2_1) * out->at(f2_1) *
    (hvec_ptr->at(f2_12) + hvec_ptr->at(f2_22)*out->at(f2_2));
  out->at(f2_22) = -out->at(f2_12)*hvec_ptr->at(f2_1) -
    out->at(f2_1)*(hvec_ptr->at(f2_11) + hvec_ptr->at(f2_12)*out->at(f2_2));
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
  out->at(1) = hvec[f2_2]; // here p is constant
  out->at(2) = hvec[f2_22]; // here p is constant
}

void f_tsl2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  entropy_liquid2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f2_2]; // here p is constant
  out->at(2) = hvec[f2_22]; // here p is constant
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
  if(pr < param::Pc[comp]){ // Could be two phase
    hvs = entropy(comp, sat_delta_v(comp, taus)->at(0), taus);
    hls = entropy(comp, sat_delta_l(comp, taus)->at(0), taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->resize(6);
      out->at(0) = taus;
      out->at(f2_1) = 0;
      out->at(f2_2) = taus_vec_ptr->at(1);
      out->at(f2_11) = 0;
      out->at(f2_12) = 0;
      out->at(f2_22) = taus_vec_ptr->at(2);
      return;
    }
  }
  std::vector<double> *hvec_ptr;
  state_solve_dat sd;
  std::vector<double> hout;
  sd.p = pr;
  sd.h = ht;
  sd.comp = comp;
  if(pr >= param::Pc[comp]){ // it's liquid or supercritical
    // Assume between melting and Tmax
    //   This is a pretty big temperature range bracketing for a good guess for
    //   Newton method may be slow here.  May want some aux functions to break
    //   it up a bit more
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tsl, tau_lo, tau_hi, &tau, 20, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tsl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_liquid(comp, pr, tau);
  }
  else if(pr < param::Pt[comp]){ // it's vapor (unless it's ice)
    // Assume between sublimation temperature and Tmax
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tsv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tsv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_vapor(comp, pr, tau);
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    // Assume between melting temperature and sat temperature
    tau_lo = taus;
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tsl, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tsl2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_entropy_liquid(comp, pr, tau);
  }
  else{ //if(ht >= hvs){ // vapor for sure
    // Assume between saturation temperature and Tmax
    tau_lo = param::Tc[comp]/param::T_max[comp];
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
  out->at(f2_1) = 1.0/hvec_ptr->at(f2_2);           // d/dh
  out->at(f2_2) = -out->at(f2_1) * hvec_ptr->at(f2_1); // d/dp, tripple product
  out->at(f2_11) = -out->at(f2_1) * out->at(f2_1) * out->at(f2_1) * hvec_ptr->at(f2_22);
  out->at(f2_12) = -out->at(f2_1) * out->at(f2_1) *
    (hvec_ptr->at(f2_12) + hvec_ptr->at(f2_22)*out->at(f2_2));
  out->at(f2_22) = -out->at(f2_12)*hvec_ptr->at(f2_1) -
    out->at(f2_1)*(hvec_ptr->at(f2_11) + hvec_ptr->at(f2_12)*out->at(f2_2));
}

/*------------------------------------------------------------------------------
  Internal energy solver function wrappers
------------------------------------------------------------------------------*/

double f_tuv(double tau, void* dat){
  // This is for bracket solver, so may want to use an underlying function that
  // doesn't bother calculating deriavtives, but I don't currently have such a
  // version.
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> out;
  internal_energy_vapor2(sd->comp, sd->p, tau, &out);
  return out[0] - sd->h;
}

double f_tul(double tau, void* dat){
  // This is for bracket solver, so may want to use an underlying function that
  // doesn't bother calculating deriavtives, but I don't currently have such a
  // version.
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> out;
  internal_energy_liquid2(sd->comp, sd->p, tau, &out);
  return out[0] - sd->h;
}

void f_tuv2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  internal_energy_vapor2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f2_2]; // here p is constant
  out->at(2) = hvec[f2_22]; // here p is constant
}

void f_tul2(double tau, std::vector<double> *out, void* dat){
  // This is for the halley method which needs second order derivatives
  state_solve_dat *sd = (state_solve_dat*)dat;
  std::vector<double> hvec;
  internal_energy_liquid2(sd->comp, sd->p, tau, &hvec);
  out->resize(3);
  out->at(0) = hvec[0] - sd->h;
  out->at(1) = hvec[f2_2]; // here p is constant
  out->at(2) = hvec[f2_22]; // here p is constant
}

/*------------------------------------------------------------------------------
  tau(u, p) solver with derivatives

  This function first checks that if at the given pressure, the enthalpy is
  in the 2-phase region.  If it is, tau is just tau_sat.  If not it tries to
  classify the region in either liquid, vapor, vapor below the tripple point,
  or supercritical. Once the region is itentified, a bracketing method is used
  to narrow down the initial guess before solving with a Newton method.
------------------------------------------------------------------------------*/

void tau_up2(comp_enum comp, double ht, double pr, std::vector<double> *out){
  std::vector<double> *taus_vec_ptr;
  double taus, hvs, hls, tau_hi, tau_lo, tau=0;

  taus_vec_ptr = sat_tau(comp, pr);
  taus = taus_vec_ptr->at(0);
  if(pr < param::Pc[comp]){ // Could be two phase
    hvs = internal_energy(comp, sat_delta_v(comp, taus)->at(0), taus);
    hls = internal_energy(comp, sat_delta_l(comp, taus)->at(0), taus);
    if(ht > hls && ht < hvs){ // two-phase
      out->resize(6);
      out->at(0) = taus;
      out->at(f2_1) = 0;
      out->at(f2_2) = taus_vec_ptr->at(1);
      out->at(f2_11) = 0;
      out->at(f2_12) = 0;
      out->at(f2_22) = taus_vec_ptr->at(2);
      return;
    }
  }
  std::vector<double> *hvec_ptr;
  state_solve_dat sd;
  std::vector<double> hout;
  sd.p = pr;
  sd.h = ht;
  sd.comp = comp;
  if(pr >= param::Pc[comp]){ // it's liquid or supercritical
    // Assume between melting and Tmax
    //   This is a pretty big temperature range bracketing for a good guess for
    //   Newton method may be slow here.  May want some aux functions to break
    //   it up a bit more
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tul, tau_lo, tau_hi, &tau, 20, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tul2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_internal_energy_liquid(comp, pr, tau);
  }
  else if(pr < param::Pt[comp]){ // it's vapor (unless it's ice)
    // Assume between sublimation temperature and Tmax
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tuv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tuv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_internal_energy_vapor(comp, pr, tau);
  }
  else if(ht <= hls){ // liquid (unless it's ice)
    // Assume between melting temperature and sat temperature
    tau_lo = taus;
    tau_hi = param::Tc[comp]/melting_temperature_func[comp](pr);
    int n1 = bracket(f_tul, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tul2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_internal_energy_liquid(comp, pr, tau);
  }
  else{ //if(ht >= hvs){ // vapor for sure
    // Assume between saturation temperature and Tmax
    tau_lo = param::Tc[comp]/param::T_max[comp];
    tau_hi = taus;
    int n1 = bracket(f_tuv, tau_lo, tau_hi, &tau, 10, 1e-4, 1e-4, &sd);
    int n2 = halley(f_tuv2, tau, &tau, &hout, 40, 1e-9, &sd);
    hvec_ptr = memo2_internal_energy_vapor(comp, pr, tau);
  }
  out->resize(6);
  out->at(0) = tau;
  // here, for indexing of **out only**:
  //   d will be derivative with resepect to h and
  //   t will be with respect to p
  out->at(f2_1) = 1.0/hvec_ptr->at(f2_2);           // d/dh
  out->at(f2_2) = -out->at(f2_1) * hvec_ptr->at(f2_1); // d/dp, tripple product
  out->at(f2_11) = -out->at(f2_1) * out->at(f2_1) * out->at(f2_1) * hvec_ptr->at(f2_22);
  out->at(f2_12) = -out->at(f2_1) * out->at(f2_1) *
    (hvec_ptr->at(f2_12) + hvec_ptr->at(f2_22)*out->at(f2_2));
  out->at(f2_22) = -out->at(f2_12)*hvec_ptr->at(f2_1) -
    out->at(f2_1)*(hvec_ptr->at(f2_11) + hvec_ptr->at(f2_12)*out->at(f2_2));
}

void vf_hp2(comp_enum comp, double ht, double pr, std::vector<double> *out){
  out->resize(6);
  if(pr >= param::Pc[comp]){ // consider supercritical fluid liquid
    out->at(0) = 0;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }
  if(pr < param::Pt[comp]){  // either vapor or ice, ignoring ice
    out->at(0) = 1;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  std::vector<double> *taus_vec_ptr, *hsat_l_vec_ptr, *hsat_v_vec_ptr;

  taus_vec_ptr = sat_tau(comp, pr);
  hsat_l_vec_ptr = memo2_enthalpy_liquid(comp, pr, taus_vec_ptr->at(0));
  hsat_v_vec_ptr = memo2_enthalpy_vapor(comp, pr, taus_vec_ptr->at(0));
  double hv = hsat_v_vec_ptr->at(0);
  double hl = hsat_l_vec_ptr->at(0);

  if(ht < hl){
    out->at(0) = 0;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  if(ht > hv){
    out->at(0) = 1;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  // if you're here then your in the saturated area
  double vf = (ht - hl)/(hv - hl);
  out->at(0) = vf;

  double dhvdp = hsat_v_vec_ptr->at(f2_1) + hsat_v_vec_ptr->at(f2_2)*taus_vec_ptr->at(1);
  double dhldp = hsat_l_vec_ptr->at(f2_1) + hsat_l_vec_ptr->at(f2_2)*taus_vec_ptr->at(1);

  out->at(f2_1) = 1.0/(hv - hl);
  out->at(f2_2) = -dhldp/(hv - hl) - vf/(hv-hl)*(dhvdp - dhldp);

  double d2hvdp2 = hsat_v_vec_ptr->at(f2_11) +
    2*hsat_v_vec_ptr->at(f2_12)*taus_vec_ptr->at(1) +
    hsat_v_vec_ptr->at(f2_22)*taus_vec_ptr->at(1)*taus_vec_ptr->at(1) +
    hsat_v_vec_ptr->at(f2_2)*taus_vec_ptr->at(2);
  double d2hldp2 = hsat_l_vec_ptr->at(f2_11) +
    2*hsat_l_vec_ptr->at(f2_12)*taus_vec_ptr->at(1) +
    hsat_l_vec_ptr->at(f2_22)*taus_vec_ptr->at(1)*taus_vec_ptr->at(1) +
    hsat_l_vec_ptr->at(f2_2)*taus_vec_ptr->at(2);
  out->at(f2_11) = 0;
  out->at(f2_12) = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
  out->at(f2_22) = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
            2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
            (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
}

void vf_sp2(comp_enum comp, double ht, double pr, std::vector<double> *out){
  out->resize(6);
  if(pr >= param::Pc[comp]){ // consider supercritical fluid liquid
    out->at(0) = 0;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }
  if(pr < param::Pt[comp]){  // either vapor or ice, ignoring ice
    out->at(0) = 1;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  std::vector<double> *taus_vec_ptr, *hsat_l_vec_ptr, *hsat_v_vec_ptr;

  taus_vec_ptr = sat_tau(comp, pr);
  hsat_l_vec_ptr = memo2_entropy_liquid(comp, pr, taus_vec_ptr->at(0));
  hsat_v_vec_ptr = memo2_entropy_vapor(comp, pr, taus_vec_ptr->at(0));
  double hv = hsat_v_vec_ptr->at(0);
  double hl = hsat_l_vec_ptr->at(0);

  if(ht < hl){
    out->at(0) = 0;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  if(ht > hv){
    out->at(0) = 1;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  // if you're here then your in the saturated area
  double vf = (ht - hl)/(hv - hl);
  out->at(0) = vf;

  double dhvdp = hsat_v_vec_ptr->at(f2_1) + hsat_v_vec_ptr->at(f2_2)*taus_vec_ptr->at(1);
  double dhldp = hsat_l_vec_ptr->at(f2_1) + hsat_l_vec_ptr->at(f2_2)*taus_vec_ptr->at(1);

  out->at(f2_1) = 1.0/(hv - hl);
  out->at(f2_2) = -dhldp/(hv - hl) - vf/(hv-hl)*(dhvdp - dhldp);

  double d2hvdp2 = hsat_v_vec_ptr->at(f2_11) +
    2*hsat_v_vec_ptr->at(f2_12)*taus_vec_ptr->at(1) +
    hsat_v_vec_ptr->at(f2_22)*taus_vec_ptr->at(1)*taus_vec_ptr->at(1) +
    hsat_v_vec_ptr->at(f2_2)*taus_vec_ptr->at(2);
  double d2hldp2 = hsat_l_vec_ptr->at(f2_11) +
    2*hsat_l_vec_ptr->at(f2_12)*taus_vec_ptr->at(1) +
    hsat_l_vec_ptr->at(f2_22)*taus_vec_ptr->at(1)*taus_vec_ptr->at(1) +
    hsat_l_vec_ptr->at(f2_2)*taus_vec_ptr->at(2);
  out->at(f2_11) = 0;
  out->at(f2_12) = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
  out->at(f2_22) = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
            2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
            (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
}

void vf_up2(comp_enum comp, double ht, double pr, std::vector<double> *out){
  out->resize(6);
  if(pr >= param::Pc[comp]){ // consider supercritical fluid liquid
    out->at(0) = 0;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }
  if(pr < param::Pt[comp]){  // either vapor or ice, ignoring ice
    out->at(0) = 1;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  std::vector<double> *taus_vec_ptr, *hsat_l_vec_ptr, *hsat_v_vec_ptr;

  taus_vec_ptr = sat_tau(comp, pr);
  hsat_l_vec_ptr = memo2_internal_energy_liquid(comp, pr, taus_vec_ptr->at(0));
  hsat_v_vec_ptr = memo2_internal_energy_vapor(comp, pr, taus_vec_ptr->at(0));
  double hv = hsat_v_vec_ptr->at(0);
  double hl = hsat_l_vec_ptr->at(0);

  if(ht < hl){
    out->at(0) = 0;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  if(ht > hv){
    out->at(0) = 1;
    out->at(f2_1) = 0;
    out->at(f2_2) = 0;
    out->at(f2_11) = 0;
    out->at(f2_12) = 0;
    out->at(f2_22) = 0;
    return;
  }

  // if you're here then your in the saturated area
  double vf = (ht - hl)/(hv - hl);
  out->at(0) = vf;

  double dhvdp = hsat_v_vec_ptr->at(f2_1) + hsat_v_vec_ptr->at(f2_2)*taus_vec_ptr->at(1);
  double dhldp = hsat_l_vec_ptr->at(f2_1) + hsat_l_vec_ptr->at(f2_2)*taus_vec_ptr->at(1);

  out->at(f2_1) = 1.0/(hv - hl);
  out->at(f2_2) = -dhldp/(hv - hl) - vf/(hv-hl)*(dhvdp - dhldp);

  double d2hvdp2 = hsat_v_vec_ptr->at(f2_11) +
    2*hsat_v_vec_ptr->at(f2_12)*taus_vec_ptr->at(1) +
    hsat_v_vec_ptr->at(f2_22)*taus_vec_ptr->at(1)*taus_vec_ptr->at(1) +
    hsat_v_vec_ptr->at(f2_2)*taus_vec_ptr->at(2);
  double d2hldp2 = hsat_l_vec_ptr->at(f2_11) +
    2*hsat_l_vec_ptr->at(f2_12)*taus_vec_ptr->at(1) +
    hsat_l_vec_ptr->at(f2_22)*taus_vec_ptr->at(1)*taus_vec_ptr->at(1) +
    hsat_l_vec_ptr->at(f2_2)*taus_vec_ptr->at(2);
  out->at(f2_11) = 0;
  out->at(f2_12) = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
  out->at(f2_22) = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
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

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
 Calaculate components of dimensionless ideal and residual Helmholtz free energy
 and derivatives to fourth order (2 for thermo + 2 for optimization solver). The
 file also contains function calls for aux curves (approx sat densities).

 Author: John Eslick
 File: phi.cpp
--------------------------------------------------------------------------------*/

#include"config.h"
#include"phi.h"
#include"read_params.h"
#include<asl.h>
#include<iostream>
#include<stdlib.h>

prop_memo_table23 memo_table_phi_resi;
prop_memo_table23 memo_table_phi_ideal;
prop_memo_table24 memo_table_phi_resi4;
prop_memo_table24 memo_table_phi_ideal4;

prop_memo_table22 memo_table_viscosity;
prop_memo_table22 memo_table_thermal_conductivity;
prop_memo_table22 memo_table_surface_tension;


f23_struct phi_resi(uint comp, double delta, double tau){
  /*
    Calculate the dimensionless resi part of Helmholtz free energy (phi^r), and
    derivatives up to fourth order. The second order derivatives are need to
    thermodynamic relations for the properties.  Then the 3rd and 4th order
    derivatives are needed to calculate first and second order derivatives of
    the properties.

    The results of this function are cached, so they don't need to be repeated
    when calculating multiple properties.  After the first calculation the AD
    tape is also reused for subsequent calculations.
  */

  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-8){
    tau = 1e-8;
  }

  try{
    return memo_table_phi_resi.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }

  real x[2] = {(real)delta, (real)tau};
  f23_struct *res_ptr;
  ASL *asl;
  real G[2] = {0};
  fint err = 0;
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;
  if(memo_table_phi_resi.size() > MAX_MEMO_PHI) memo_table_phi_resi.clear();
  res_ptr = &memo_table_phi_resi[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);

  res_ptr->f = (double)objval(dat->expr_map[phir], x, &err);
  res_ptr->f_1 = (double)objval(dat->expr_map[phir_d], x, &err);
  res_ptr->f_11 = (double)objval(dat->expr_map[phir_dd], x, &err);
  res_ptr->f_2 = (double)objval(dat->expr_map[phir_t], x, &err);
  res_ptr->f_22 = (double)objval(dat->expr_map[phir_tt], x, &err);
  res_ptr->f_12 = (double)objval(dat->expr_map[phir_dt], x, &err);
  objgrd(dat->expr_map[phir_dd], x, G, &err);
  res_ptr->f_111 = (double)G[0];
  res_ptr->f_112 = (double)G[1];
  objgrd(dat->expr_map[phir_tt], x, G, &err);
  res_ptr->f_122 = (double)G[0];
  res_ptr->f_222 = (double)G[1];

  return *res_ptr;
}

f24_struct phi_resi4(uint comp, double delta, double tau){
  /*
    Calculate the dimensionless resi part of Helmholtz free energy (phi^r), and
    derivatives up to fourth order. The second order derivatives are need to
    thermodynamic relations for the properties.  Then the 3rd and 4th order
    derivatives are needed to calculate first and second order derivatives of
    the properties.

    The results of this function are cached, so they don't need to be repeated
    when calculating multiple properties.  After the first calculation the AD
    tape is also reused for subsequent calculations.
  */
  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-8){
    tau = 1e-8;
  }

  try{
    return memo_table_phi_resi4.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }

  real x[2] = {(real)delta, (real)tau};
  f24_struct *res_ptr;
  ASL *asl;
  real G[2] = {0};  // There may be two or three vars so accomadate 3
  real H[3] = {0};
  fint err = 0;
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;
  if(memo_table_phi_resi4.size() > MAX_MEMO_PHI) memo_table_phi_resi4.clear();
  res_ptr = &memo_table_phi_resi4[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);

  res_ptr->f = (double)objval(dat->expr_map[phir], x, &err);
  res_ptr->f_1 = (double)objval(dat->expr_map[phir_d], x, &err);
  res_ptr->f_11 = (double)objval(dat->expr_map[phir_dd], x, &err);
  res_ptr->f_2 = (double)objval(dat->expr_map[phir_t], x, &err);
  res_ptr->f_22 = (double)objval(dat->expr_map[phir_tt], x, &err);
  res_ptr->f_12 = (double)objval(dat->expr_map[phir_dt], x, &err);
  
  objgrd(dat->expr_map[phir_dd], x, G, &err);

  res_ptr->f_111 = (double)G[0];
  res_ptr->f_112 = (double)G[1];
  objgrd(dat->expr_map[phir_tt], x, G, &err);

  res_ptr->f_122 = (double)G[0];
  res_ptr->f_222 = (double)G[1];

  duthes(H, dat->expr_map[phir_tt], nullptr, nullptr);

  res_ptr->f_1222 = (double)H[1];
  res_ptr->f_2222 = (double)H[2];

  duthes(H, dat->expr_map[phir_dd], nullptr, nullptr);

  res_ptr->f_1111 = (double)H[0];
  res_ptr->f_1112 = (double)H[1];
  res_ptr->f_1122 = (double)H[2];

  return *res_ptr;
}

f23_struct phi_ideal(uint comp, double delta, double tau){
  /*
    Calculate the dimensionless ideal part of Helmholtz free energy (phi^0), and
    derivatives up to fourth order. The second order derivatives are need to
    thermodynamic relations for the properties.  Then the 3rd and 4th order
    derivatives are needed to calculate first and second order derivatives of
    the properties.

    The results of this function are cached, so they don't need to be repeated
    when calculating multiple properties.  After the first calculation the AD
    tape is also reused for subsequent calculations.
  */
  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-2){
    tau = 1e-2;
  }
  try{
    return memo_table_phi_ideal.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }
  real x[2] = {(real)delta, (real)tau};
  f23_struct *res_ptr;
  ASL *asl;
  fint err = 0;
  real G[2] = {0};
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;

  if(memo_table_phi_ideal.size() > MAX_MEMO_PHI) memo_table_phi_ideal.clear();
  res_ptr = &memo_table_phi_ideal[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);
  res_ptr->f = (double)objval(dat->expr_map[phii], x, &err);
  res_ptr->f_1 = (double)objval(dat->expr_map[phii_d], x, &err);
  res_ptr->f_11 = (double)objval(dat->expr_map[phii_dd], x, &err);
  res_ptr->f_2 = (double)objval(dat->expr_map[phii_t], x, &err);
  res_ptr->f_22 = (double)objval(dat->expr_map[phii_tt], x, &err);
  res_ptr->f_12 = (double)objval(dat->expr_map[phii_dt], x, &err);

  objgrd(dat->expr_map[phii_dd], x, G, &err);
  res_ptr->f_111 = (double)G[0];
  res_ptr->f_112 = (double)G[1];
  objgrd(dat->expr_map[phii_tt], x, G, &err);
  res_ptr->f_122 = (double)G[0];
  res_ptr->f_222 = (double)G[1];
  return *res_ptr;
}

f24_struct phi_ideal4(uint comp, double delta, double tau){
  /*
    Calculate the dimensionless ideal part of Helmholtz free energy (phi^0), and
    derivatives up to fourth order. The second order derivatives are need to
    thermodynamic relations for the properties.  Then the 3rd and 4th order
    derivatives are needed to calculate first and second order derivatives of
    the properties.

    The results of this function are cached, so they don't need to be repeated
    when calculating multiple properties.  After the first calculation the AD
    tape is also reused for subsequent calculations.
  */
  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-2){
    tau = 1e-2;
  }
  try{
    return memo_table_phi_ideal4.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }
  real x[2] = {(real)delta, (real)tau};
  f24_struct *res_ptr;
  ASL *asl;
  fint err = 0;
  real G[2] = {0};
  real H[3] = {0};
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;

  if(memo_table_phi_ideal4.size() > MAX_MEMO_PHI) memo_table_phi_ideal4.clear();
  res_ptr = &memo_table_phi_ideal4[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);
  res_ptr->f = (double)objval(dat->expr_map[phii], x, &err);
  res_ptr->f_1 = (double)objval(dat->expr_map[phii_d], x, &err);
  res_ptr->f_11 = (double)objval(dat->expr_map[phii_dd], x, &err);
  res_ptr->f_2 = (double)objval(dat->expr_map[phii_t], x, &err);
  res_ptr->f_22 = (double)objval(dat->expr_map[phii_tt], x, &err);
  res_ptr->f_12 = (double)objval(dat->expr_map[phii_dt], x, &err);

  objgrd(dat->expr_map[phii_dd], x, G, &err);
  res_ptr->f_111 = (double)G[0];
  res_ptr->f_112 = (double)G[1];
  objgrd(dat->expr_map[phii_tt], x, G, &err);
  res_ptr->f_122 = (double)G[0];
  res_ptr->f_222 = (double)G[1];

  duthes(H, dat->expr_map[phii_dd], nullptr, nullptr);
  res_ptr->f_1111 = (double)H[0];
  res_ptr->f_1112 = (double)H[1];
  res_ptr->f_1122 = (double)H[2];

  duthes(H, dat->expr_map[phii_tt], nullptr, nullptr);
  res_ptr->f_1222 = (double)H[1];
  res_ptr->f_2222 = (double)H[2];
  return *res_ptr;
}

f12_struct phi_resi_for_sat(uint comp, double delta, double tau){
  /*
    Calculate the dimensionless resi part of Helmholtz free energy (phi^r), and
    derivatives second order with respect to delta. After the first calculation
    the AD tape is reused for subsequent calculations. This is not cached, since
    it is used specifically in solving for the saturation curve calculations.
  */
  real x[2] = {(real)delta, (real)tau};
  ASL *asl;
  fint err = 0;
  f12_struct res;
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;
  asl = (ASL*)(dat->asl);
  res.f = (double)objval(dat->expr_map[phir], x, &err);
  res.f_1 = (double)objval(dat->expr_map[phir_d], x, &err);
  res.f_11 = (double)objval(dat->expr_map[phir_dd], x, &err);
  return res;
}

double sat_delta_v_approx(uint comp, double tau){
  real x[2] = {(real)1.0, (real)tau}; // delta doesn't matter here
  ASL *asl;
  fint err = 0;
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;
  asl = (ASL*)(dat->asl);
  return (double)objval(dat->expr_map[delta_v_sat_approx], x, &err);
}

double sat_delta_l_approx(uint comp, double tau){
  real x[2] = {(real)1.0, (real)tau}; // delta doesn't matter here
  ASL *asl;
  fint err = 0;
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;
  asl = (ASL*)(dat->asl);
  return (double)objval(dat->expr_map[delta_l_sat_approx], x, &err);
}

f22_struct memo2_viscosity(uint comp, double delta, double tau){
  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-2){
    tau = 1e-2;
  }
  try{
    return memo_table_viscosity.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }
  real x[2] = {(real)delta, (real)tau};
  f22_struct *res_ptr;
  ASL *asl;
  fint err = 0;
  real G[2] = {0};
  real H[3] = {0};
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;

  if(memo_table_viscosity.size() > MAX_MEMO_PHI) memo_table_viscosity.clear();
  res_ptr = &memo_table_viscosity[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);
  if (dat->expr_map[viscosity_idx]!=1000){
    objgrd(dat->expr_map[viscosity_idx], x, G, &err);
    duthes(H, dat->expr_map[viscosity_idx], nullptr, nullptr);
    res_ptr->f = (double)objval(dat->expr_map[viscosity_idx], x, &err);
    res_ptr->f_1 = (double)G[0];
    res_ptr->f_2 = (double)G[1];
    res_ptr->f_11 = (double)H[0];
    res_ptr->f_12 = (double)H[1];
    res_ptr->f_22 = (double)H[2];
  }
  else{
    res_ptr->f = 0.0;
    res_ptr->f_1 = 0.0;
    res_ptr->f_2 = 0.0;
    res_ptr->f_11 = 0.0;
    res_ptr->f_12 = 0.0;
    res_ptr->f_22 = 0.0;
  }
  return *res_ptr;
}

f22_struct memo2_thermal_conductivity(uint comp, double delta, double tau){
  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-2){
    tau = 1e-2;
  }
  try{
    return memo_table_thermal_conductivity.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }
  real x[2] = {(real)delta, (real)tau};
  f22_struct *res_ptr;
  ASL *asl;
  fint err = 0;
  real G[2] = {0};
  real H[3] = {0};
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;

  if(memo_table_thermal_conductivity.size() > MAX_MEMO_PHI) memo_table_thermal_conductivity.clear();
  res_ptr = &memo_table_thermal_conductivity[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);
  if (dat->expr_map[thermal_conductivity_idx]!=1000){
    objgrd(dat->expr_map[thermal_conductivity_idx], x, G, &err);
    duthes(H, dat->expr_map[thermal_conductivity_idx], nullptr, nullptr);
    res_ptr->f = (double)objval(dat->expr_map[thermal_conductivity_idx], x, &err);
    res_ptr->f_1 = (double)G[0];
    res_ptr->f_2 = (double)G[1];
    res_ptr->f_11 = (double)H[0];
    res_ptr->f_12 = (double)H[1];
    res_ptr->f_22 = (double)H[2];
  }
  else{
    res_ptr->f = 0.0;
    res_ptr->f_1 = 0.0;
    res_ptr->f_2 = 0.0;
    res_ptr->f_11 = 0.0;
    res_ptr->f_12 = 0.0;
    res_ptr->f_22 = 0.0;
  }
  return *res_ptr;
}

f22_struct memo2_surface_tension(uint comp, double delta, double tau){
  if (delta < 1e-14){
    delta = 1e-14;
  }
  if (tau < 1e-2){
    tau = 1e-2;
  }
  try{
    return memo_table_surface_tension.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range const&){
  }
  real x[2] = {(real)delta, (real)tau};
  f22_struct *res_ptr;
  ASL *asl;
  fint err = 0;
  real G[2] = {0};
  real H[3] = {0};
  parameters_struct *dat = &cdata[comp];
  using namespace expr_idx;

  if(memo_table_surface_tension.size() > MAX_MEMO_PHI) memo_table_surface_tension.clear();
  res_ptr = &memo_table_surface_tension[std::make_tuple(comp, delta, tau)];
  asl = (ASL*)(dat->asl);
  if (dat->expr_map[surface_tension_idx]!=1000){
    objgrd(dat->expr_map[surface_tension_idx], x, G, &err);
    duthes(H, dat->expr_map[surface_tension_idx], nullptr, nullptr);
    res_ptr->f = (double)objval(dat->expr_map[surface_tension_idx], x, &err);
    res_ptr->f_1 = (double)G[0];
    res_ptr->f_2 = (double)G[1];
    res_ptr->f_11 = (double)H[0];
    res_ptr->f_12 = (double)H[1];
    res_ptr->f_22 = (double)H[2];
  }
  else{
    res_ptr->f = 0.0;
    res_ptr->f_1 = 0.0;
    res_ptr->f_2 = 0.0;
    res_ptr->f_11 = 0.0;
    res_ptr->f_12 = 0.0;
    res_ptr->f_22 = 0.0;
  }
  return *res_ptr;
}
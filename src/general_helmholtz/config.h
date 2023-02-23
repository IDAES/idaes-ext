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
 This file provides some configuration parameters.

 Author: John Eslick
 File: config.h
--------------------------------------------------------------------------------*/

#include<vector>
#include<unordered_map>
#include<boost/functional/hash.hpp>
#include<climits>
#include<math.h>

#ifndef _INCLUDE_CONFIG_H_
#define _INCLUDE_CONFIG_H_

#define MAX_MEMO_PHI 500000
#define MAX_MEMO_PROP 1000000

typedef unsigned int uint;
typedef unsigned char uchar;

namespace expr_idx{
  const long phii = 0;
  const long phii_d = 1;
  const long phii_dd = 2;
  const long phii_t = 3;
  const long phii_tt = 4;
  const long phii_dt = 5;
  const long phir = 6;
  const long phir_d = 7;
  const long phir_dd = 8;
  const long phir_t = 9;
  const long phir_tt = 10;
  const long phir_dt = 11;
  const long delta_v_sat_approx = 12;
  const long delta_l_sat_approx = 13;
  const long viscosity_idx = 14;
  const long thermal_conductivity_idx = 15;
  const long surface_tension_idx = 16;
}
const long expr_map_size = 17;

// Structure for unary function value and 1st order derivatives
struct f11_struct {
  double f;
  double f_1;
};

// Structure for unary function value and 2nd order derivatives
struct f12_struct {
  double f;
  double f_1;
  double f_11;
};

// Structure for binary function value and 1st order derivatives
struct f21_struct {
  double f;
  double f_1;
  double f_2;
};

// Structure for binary function value and 2nd order derivatives
struct f22_struct {
  double f;
  double f_1;
  double f_11;
  double f_2;
  double f_12;
  double f_22;
};

// Structure for binary function value and 3rd order derivatives
struct f23_struct {
  double f;
  double f_1;
  double f_11;
  double f_111;
  double f_2;
  double f_12;
  double f_112;
  double f_22;
  double f_122;
  double f_222;
};

// Structure for binary function value and 4th order derivatives
struct f24_struct {
  double f;
  double f_1;
  double f_11;
  double f_111;
  double f_1111;
  double f_2;
  double f_12;
  double f_112;
  double f_1112;
  double f_22;
  double f_122;
  double f_1122;
  double f_222;
  double f_1222;
  double f_2222;
};

// Structure for parameters and expressions for specific component
struct parameters_struct {
  void *asl = nullptr;// expressions in ASL (cast to ASL*)
  long expr_map[expr_map_size];  // index of expressions in NL file
  long var_map[3];    // index of variables in NL file
  double R;           // specific ideal gas constant [kJ/kg/K]
  double MW;          // molecular weight [g/mol]
  double T_star;      // for calculating tau = T_star/T [K]
  double rho_star;    // for calculating delta = rho_star/rho [kg/m3]
  double Tc;          // critical temperature [K]
  double rhoc;        // critical density [kg/m3]
  double Pc;          // critical pressure [kPa]
  double Tt;          // triple point temperature [K]
  double Pt;          // triple point pressure [K]
  double rhot_l;      // liquid triple point density [kg/m3]
  double rhot_v;      // vapor triple point density [kg/m3]
  double P_min;       // minimum pressure [kPa]
  double P_max;       // maximum pressure [kPa]
  double rho_max;     // maximum density [kg/m3]
  double T_min;       // minimum temperature [kPa]
  double T_max;       // maximum temperature [kPa]
};

// AMPL models for a component
const uint MISSING_DATA = UINT_MAX;
extern std::unordered_map<std::string, uint> cindex;
extern std::vector<parameters_struct> cdata;

inline double tau_c(uint comp){
  return cdata[comp].T_star/cdata[comp].Tc;
}

inline double delta_c(uint comp){
  return cdata[comp].rhoc/cdata[comp].rho_star;
}

typedef std::unordered_map<
  std::tuple<uint, double>,
  f11_struct,
  boost::hash<std::tuple<uint, double>>
> prop_memo_table11;

typedef std::unordered_map<
  std::tuple<uint, double>,
  f12_struct,
  boost::hash<std::tuple<uint, double>>
> prop_memo_table12;

typedef std::unordered_map<
  std::tuple<uint, double, double>,
  f21_struct,
  boost::hash<std::tuple<uint, double, double>>
> prop_memo_table21;

typedef std::unordered_map<
  std::tuple<uint, double, double>,
  f22_struct,
  boost::hash<std::tuple<uint, double, double>>
> prop_memo_table22;

typedef std::unordered_map<
  std::tuple<uint, double, double>,
  f23_struct,
  boost::hash<std::tuple<uint, double, double>>
> prop_memo_table23;

typedef std::unordered_map<
  std::tuple<uint, double, double>,
  f24_struct,
  boost::hash<std::tuple<uint, double, double>>
> prop_memo_table24;

#define MEMO2_FUNCTION(new_func, calc_func, table) \
f22_struct new_func(uint comp, double delta, double tau){ \
  if(std::isnan(delta) || std::isnan(tau)){ \
    f22_struct res; \
    res.f = nan(""); \
    res.f_1 = nan(""); \
    res.f_2 = nan(""); \
    res.f_11 = nan(""); \
    res.f_22 = nan(""); \
    res.f_12 = nan(""); \
    return res; \
  } \
  try{ \
    return table.at(std::make_tuple(comp, delta, tau)); \
  } \
  catch(std::out_of_range const&){ \
  } \
  if(table.size() > MAX_MEMO_PROP) table.clear(); \
  f22_struct *res_ptr; \
  res_ptr = &table[std::make_tuple(comp, delta, tau)]; \
  calc_func(comp, delta, tau, res_ptr); \
  return *res_ptr; \
}

#endif

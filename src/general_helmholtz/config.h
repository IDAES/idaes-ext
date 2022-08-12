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
#include<boost/functional/hash.hpp>
#include"components/param.h"
#include"components/function_pointers.h"

#ifndef _INCLUDE_CONFIG_H_
#define _INCLUDE_CONFIG_H_

#define MAX_MEMO_PHI 500000
#define MAX_MEMO_PROP 1000000

typedef unsigned int uint;
typedef unsigned char uchar;

enum deriv1_enum {
  f1 = 0,
  f1_1 = 1,
  f1_2 = 2,
};

enum deriv2_enum {
  f2 = 0,
  f2_1 = 1,
  f2_11 = 2,
  f2_2 = 3,
  f2_12 = 4,
  f2_22 = 5,
};

enum deriv4_enum {
  f4 = 0,
  f4_1 = 1,
  f4_11 = 2,
  f4_111 = 3,
  f4_1111 = 4,
  f4_2 = 5,
  f4_12 = 6,
  f4_112 = 7,
  f4_1112 = 8,
  f4_22 = 9,
  f4_122 = 10,
  f4_1122 = 11,
  f4_222 = 12,
  f4_1222 = 13,
  f4_2222 = 14,
};

#define PHI_IDEAL_TAPE(c) c
#define PHI_RESI_TAPE(c) NCOMPS + c

extern unsigned int taped_ideal[NCOMPS];
extern unsigned int taped_resi[NCOMPS];

static std::vector<double> nan_vec3 = {
  nan(""),
  nan(""),
  nan("")
};

static std::vector<double> zero_vec3 = {
  0.0,
  0.0,
  0.0
};

static std::vector<double> one_vec3 = {
  1.0,
  1.0,
  1.0
};

static std::vector<double> nan_vec6 = {
  nan(""),
  nan(""),
  nan(""),
  nan(""),
  nan(""),
  nan("")
};

static std::vector<double> zero_vec6 = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
};

static std::vector<double> one_vec6 = {
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
};

typedef std::unordered_map<
  std::tuple<comp_enum, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double>>
> prop_memo_table1;

typedef std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> prop_memo_table2;

#define MEMO2_FUNCTION(new_func, calc_func, table) \
std::vector<double> *new_func(comp_enum comp, double delta, double tau){ \
  if(std::isnan(delta) || std::isnan(tau)) return &nan_vec6; \
  try{ \
    return &table.at(std::make_tuple(comp, delta, tau)); \
  } \
  catch(std::out_of_range const&){ \
  } \
  std::vector<double> *yvec_ptr; \
  if(table.size() > MAX_MEMO_PROP) table.clear(); \
  yvec_ptr = &table[std::make_tuple(comp, delta, tau)]; \
  calc_func(comp, delta, tau, yvec_ptr); \
  return yvec_ptr; \
}

#endif

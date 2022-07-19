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

#include<unordered_map>
#include<string>

#ifndef _INCLUDE_CONFIG_H_
#define _INCLUDE_CONFIG_H_

#define MAX_MEMO_PHI 500000
#define MAX_MEMO_PROP 1000000

typedef unsigned int uint;
typedef unsigned char uchar;

#define NCOMPS 4

// if adding components update NCOMPS above should be 1 more than last index
enum comp_enum{
  h2o = 1,
  co2 = 2,
  r134a = 3,
};

static std::unordered_map<std::string, comp_enum> const comp_string_table = {
  {"h2o", comp_enum::h2o},
  {"co2", comp_enum::co2},
  {"r134a", comp_enum::r134a},
};

enum class deriv1_enum {
  f = 0,
  f_d = 1,
  f_t = 2,
};

enum class deriv2_enum {
  f = 0,
  f_d = 1,
  f_dd = 2,
  f_t = 3,
  f_dt = 4,
  f_tt = 5,
};

enum class deriv4_enum {
  f = 0,
  f_d = 1,
  f_dd = 2,
  f_ddd = 3,
  f_dddd = 4,
  f_t = 5,
  f_dt = 6,
  f_ddt = 7,
  f_dddt = 8,
  f_tt = 9,
  f_dtt = 10,
  f_ddtt = 11,
  f_ttt = 12,
  f_dttt = 13,
  f_tttt = 14,
};

#define PHI_IDEAL_TAPE_H2O 1
#define PHI_RESI_TAPE_H2O 2

extern unsigned int taped_ideal[NCOMPS];
extern unsigned int taped_resi[NCOMPS];

#endif

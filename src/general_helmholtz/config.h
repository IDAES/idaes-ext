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

#define PHI_IDEAL_TAPE_H2O 1
#define PHI_RESI_TAPE_H2O 2

extern unsigned int taped_ideal[NCOMPS];
extern unsigned int taped_resi[NCOMPS];

#endif

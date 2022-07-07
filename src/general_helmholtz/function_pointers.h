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

#include"config.h"
#include"h2o.h"


#ifndef _INCLUDE_FUNCTION_POINTERS_H_
#define _INCLUDE_FUNCTION_POINTERS_H_

typedef double (*uni_double_function_type)(double);
typedef double (*bin_double_function_type)(double);
typedef void (*zero_void_function_type)();

/*
In the function pointer arrays the indexing needs to match the comp_enum
defined in config.h
*/
const uni_double_function_type delta_l_sat_guess_func[] = {
  (uni_double_function_type)NULL,  // 0 - not used
  delta_sat_l_approx_h2o,          // 1 - h2o
};

const uni_double_function_type delta_v_sat_guess_func[] = {
  (uni_double_function_type)NULL,  // 0 - not used
  delta_sat_v_approx_h2o,          // 1 - h2o
};

const zero_void_function_type phi_resi_tape_func[] = {
  (zero_void_function_type)NULL,   // 0 - not used
  phi_h2o_resi_tape,               // 1 - h2o
};

const zero_void_function_type phi_ideal_tape_func[] = {
  (zero_void_function_type)NULL,   // 0 - not used
  phi_h2o_ideal_tape,              // 1 - h2o
};

#endif

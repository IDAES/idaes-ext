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
 Pointers to component specific functions

 Author: John Eslick
 File: function_pointers.h
--------------------------------------------------------------------------------*/

#include<cstdlib>

#include"h2o.h"
#include"r1234ze.h"
#include"co2.h"
#include"r134a.h"

#ifndef _INCLUDE_FUNCTION_POINTERS_H_
#define _INCLUDE_FUNCTION_POINTERS_H_

typedef double (*uni_double_function_type)(double);
typedef double (*bin_double_function_type)(double);
typedef void (*zero_void_function_type)();

/*
In the function pointer arrays the indexing needs to match the comp_enum
defined in config.h
*/

// This provides melting point temperature as a function of pressure
// these functions are only to asssit the solver and do not need to be
// highly accurate
const uni_double_function_type melting_tau_func[] = {
  (uni_double_function_type)NULL,  // 0 - not used
  melting_tau_h2o,         // 1 - h2o
  melting_tau_co2,         // 2 - co2
  melting_tau_r1234ze,     // 3 - r1234ze
  melting_tau_r134a,       // 4 - r134a
};

// This provides melting liquid density as a function of pressure
// these functions are only to asssit the solver and do not need to be
// highly accurate
const uni_double_function_type melting_liquid_delta_func[] = {
  (uni_double_function_type)NULL,  // 0 - not used
  melting_liquid_delta_h2o,      // 1 - h2o
  melting_liquid_delta_co2,      // 2 - co2
  melting_liquid_delta_r1234ze,  // 3 - r1234ze
  melting_liquid_delta_r134a,       // 4 - r134a
};

// This provides a guess for saturated liquid density as a function of tau
const uni_double_function_type delta_l_sat_guess_func[] = {
  (uni_double_function_type)NULL,  // 0 - not used
  delta_sat_l_approx_h2o,          // 1 - h2o
  delta_sat_l_approx_co2,          // 2 - co2
  delta_sat_l_approx_r1234ze,      // 3 - r1234ze
  delta_sat_l_approx_r134a,        // 4 - r134a

};

// This provides a guess for saturated liquid density as a function of tau
const uni_double_function_type delta_v_sat_guess_func[] = {
  (uni_double_function_type)NULL,  // 0 - not used
  delta_sat_v_approx_h2o,          // 1 - h2o
  delta_sat_v_approx_co2,          // 2 - co2
  delta_sat_v_approx_r1234ze,      // 3 - r1234ze
  delta_sat_v_approx_r134a,        // 4 - r134a
};

// Make AD/calculation tape for residual dimensionless Helmholtz free energy
const zero_void_function_type phi_resi_tape_func[] = {
  (zero_void_function_type)NULL,   // 0 - not used
  phi_h2o_resi_tape,               // 1 - h2o
  phi_co2_resi_tape,               // 2 - co2
  phi_r1234ze_resi_tape,           // 3 - r1234ze
  phi_r134a_resi_tape,             // 4 - r134a

};

// Make AD/calculation tape for ideal dimensionless Helmholtz free energy
const zero_void_function_type phi_ideal_tape_func[] = {
  (zero_void_function_type)NULL,   // 0 - not used
  phi_h2o_ideal_tape,              // 1 - h2o
  phi_co2_ideal_tape,              // 2 - co2
  phi_r1234ze_ideal_tape,          // 3 - r1234ze
  phi_r134a_ideal_tape,            // 4 - r134a
};

#endif

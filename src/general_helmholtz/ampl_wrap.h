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
 Provide the AMPL user-defined function wrapper for property functions

 Author: John Eslick
 File: ampl_wrap.h
--------------------------------------------------------------------------------*/

#ifndef _AMPL_WRAP_H_
#define _AMPL_WRAP_H_

#include <funcadd.h>

#undef getenv
#undef printf
#undef fprintf
#undef strtod
#undef vsnprintf
#undef vsprintf
#undef vfprintf
#undef sprintf
#undef snprintf

#define ASL_WRAP_FUNC_2ARG(new_func, calc_func) \
double new_func(arglist *al){ \
  std::vector<double> *out; \
  out = calc_func(comp_string_table[al->sa[0]], al->ra[0], al->ra[1]); \
  if(al->derivs != NULL){ \
    al->derivs[0] = out->at(f2_1); \
    al->derivs[1] = out->at(f2_2); \
  } \
  if(al->hes != NULL){ \
    al->hes[0] = out->at(f2_11); \
    al->hes[1] = out->at(f2_12); \
    al->hes[2] = out->at(f2_22); \
  } \
  return out->at(0); \
}

#define ASL_WRAP_FUNC_1ARG(new_func, calc_func) \
double new_func(arglist *al){ \
  std::vector<double> *out; \
  out = calc_func(comp_string_table[al->sa[0]], al->ra[0]); \
  if(al->derivs != NULL){ \
    al->derivs[0] = out->at(1); \
  } \
  if(al->hes != NULL){ \
    al->hes[0] = out->at(2); \
  } \
  return out->at(0); \
}

#define ASL_WRAP_FUNC_0ARG(new_func, parameter) \
double new_func(arglist *al){ \
  return parameter[comp_string_table[al->sa[0]]]; \
}

#endif

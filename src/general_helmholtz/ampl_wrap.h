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

#include "read_params.h"
#include "config.h"
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
  f22_struct out; \
  std::string data_path(""); \
  if (al->n - al->nr > 1){ \
    data_path.assign(al->sa[1]); \
  } \
  uint comp = read_params(al->sa[0], data_path); \
  if (comp==MISSING_DATA){ \
    out.f = nan("missing data"); \
    out.f_1 = nan("missing data"); \
     out.f_2 = nan("missing data"); \
    out.f_11 = nan("missing data"); \
    out.f_12 = nan("missing data"); \
    out.f_22 = nan("missing data"); \
  } \
  else{ \
    out = calc_func(comp, al->ra[0], al->ra[1]); \
  } \
  if(al->derivs != NULL){ \
    al->derivs[0] = out.f_1; \
    al->derivs[1] = out.f_2; \
  } \
  if(al->hes != NULL){ \
    al->hes[0] = out.f_11; \
    al->hes[1] = out.f_12; \
    al->hes[2] = out.f_22; \
  } \
  return out.f; \
}

#define ASL_WRAP_FUNC_1ARG(new_func, calc_func) \
double new_func(arglist *al){ \
  f12_struct out; \
  std::string data_path(""); \
  if (al->n - al->nr > 1){ \
    data_path.assign(al->sa[1]); \
  } \
  uint comp = read_params(al->sa[0], data_path); \
  if (comp==MISSING_DATA){ \
    out.f = nan("missing data"); \
    out.f_1 = nan("missing data"); \
    out.f_11 = nan("missing data"); \
  } \
  else{ \
    out = calc_func(comp, al->ra[0]); \
  } \
  if(al->derivs != NULL){ \
    al->derivs[0] = out.f_1; \
  } \
  if(al->hes != NULL){ \
    al->hes[0] = out.f_11; \
  } \
  return out.f; \
}

#define ASL_WRAP_FUNC_0ARG(new_func, parameter) \
double new_func(arglist *al){ \
  std::string data_path(""); \
  if (al->n - al->nr > 1){ \
    data_path.assign(al->sa[1]); \
  } \
  uint comp = read_params(al->sa[0], data_path); \
  if (comp==MISSING_DATA){ \
    return nan("missing data"); \
  } \
  return cdata[comp].parameter; \
}

#endif

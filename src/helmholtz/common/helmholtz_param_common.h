#ifndef _INCLUDE_HELMHOLTZ_PARAM_COMMON_H_
#define _INCLUDE_HELMHOLTZ_PARAM_COMMON_H_

#include"helmholtz_param_common.h"

extern const s_real *n0;
extern const s_real *gamma0;
extern const s_real *c = param;
extern const s_real *d = param;
extern const s_real *t = param;
extern const s_real *n = param;
extern const s_real *alpha = param;
extern const s_real *theta = param;
extern const s_real *eps = param;
extern const s_real *beta = param;

// Functions to convert to dimensionless state variables from T and Rho.
// temperature and density
inline s_real f_tau(s_real T){return T_c/T;}
inline s_real f_delta(s_real rho){return rho/rho_c;}

#endif

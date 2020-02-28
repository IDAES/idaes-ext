#include"helmholtz_param.h"
#include"helmholtz_param_common.h"

const s_real *n0 = param + n0_offset;
const s_real *gamma0 = param + gamma0_offset;
const s_real *c = param + c_offset;
const s_real *d = param + d_offset;
const s_real *t = param + t_offset;
const s_real *n = param + n_offset;
const s_real *alpha = param + alpha_offset;
const s_real *theta = param + theta_offset;
const s_real *eps = param + eps_offset;
const s_real *beta = param + beta_offset;

inline s_real f_tau(s_real T){return T_c/T;}
inline s_real f_delta(s_real rho){return rho/rho_c;}

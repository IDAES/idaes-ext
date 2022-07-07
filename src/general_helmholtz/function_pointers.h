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

const zero_void_function_type phi_real_tape_func[] = {
  (zero_void_function_type)NULL,   // 0 - not used
  phi_h2o_real_tape,               // 1 - h2o
};

const zero_void_function_type phi_ideal_tape_func[] = {
  (zero_void_function_type)NULL,   // 0 - not used
  phi_h2o_ideal_tape,              // 1 - h2o
};

#endif

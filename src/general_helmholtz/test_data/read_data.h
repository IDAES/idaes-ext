#include "../config.h"

#ifndef _INCLUDE_READ_DATA_H_
#define _INCLUDE_READ_DATA_H_


namespace test_data {
  enum col_enum{
    T_col = 0,
    P_col = 1,
    rho_col = 2,
    v_col = 3,
    u_col = 4,
    h_col = 5,
    s_col = 6,
    cv_col = 7,
    cp_col = 8,
    w_col = 9,
    last_col = 13,
  };

  enum data_set_enum{
    vapor_set = 0,
    liquid_set = 1,
    saturated_set = 2,
    supercritical_set = 3,
    mixed_set = 4, // data for phases mixed in the same file
  };

}

std::vector< std::vector<double> > read_data(comp_enum comp, test_data::data_set_enum data_set);



#endif

#include<string>
#include<vector>
#include"config.h"

#ifndef _INCLUDE_READ_PARAMS_H_
#define _INCLUDE_READ_PARAMS_H_


struct tests_struct{
    std::string comp_str;
    double u_off = 0;
    double h_off = 0;
    double s_off = 0;
    std::string test_set;
};

uint read_params(std::string comp, std::string data_path="");
std::vector<tests_struct> read_run_tests(void);

#endif
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
#include"read_data.h"

#ifndef _INCLUDE_TEST_UTIL_H_
#define _INCLUDE_TEST_UTIL_H_

typedef f12_struct (*test_fptr1)(uint comp, double x);
typedef f22_struct (*test_fptr2)(uint comp, double x1, double x2);

int fd1(test_fptr1 func, uint comp, double x, double h, double tv, double tol, bool dbg);
int fd2(test_fptr2 func, uint comp, double x1, double x2, double h1, double h2, double tv, double tol, bool dbg);

uint test_basic_properties(uint comp, std::string comp_str, test_data::data_set_enum data_set, double u_off=0, double h_off=0, double s_off=0);
uint test_sat_curve(uint comp, std::string comp_str, double u_off=0, double h_off=0, double s_off=0);
uint test_delta_function(uint comp, std::string comp_str, test_data::data_set_enum data_set, double u_off=0, double h_off=0, double s_off=0);
uint test_state(uint comp, std::string comp_str, test_data::data_set_enum data_set, double u_off=0, double h_off=0, double s_off=0);
uint test_sat_curve_more(uint comp, std::string comp_str, double u_off=0, double h_off=0, double s_off=0);

uint run_set_all(uint comp, std::string comp_str, double u_off=0, double h_off=0, double s_off=0);
uint run_set_mixed(uint comp, std::string comp_str, double u_off=0, double h_off=0, double s_off=0);



#endif
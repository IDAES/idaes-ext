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

#define TEST_FUNCTION_OF_DELTA_TAU(PROP, FUNC, DAT, TOL, IGNORE_ERROR) \
start = std::chrono::high_resolution_clock::now(); \
std::cout << "    " << PROP << "(" << comp_str << ", delta, tau) "; \
for(i=0; i<dat.size(); ++i){ \
    tau = pdat->T_star/dat[i][test_data::T_col]; \
    delta = dat[i][test_data::rho_col]/pdat->rho_star; \
    err = fd2(FUNC, comp, delta, tau, 1e-4, 1e-4, DAT, TOL, 0); \
    if(err){ \
        std::cout << std::endl; \
        std::cout << "-------------------------ERROR--" << err << "------------------------" << std::endl; \
        std::cout << "Property: " << PROP << " Comp: " << comp_str << std::endl; \
        std::cout << "density " << dat[i][test_data::rho_col] << " delta: " << delta; \
        std::cout << " tau: " << tau << " pressure: " << dat[i][test_data::P_col]*1000; \
        std::cout << ", T= " << dat[i][test_data::T_col] <<  std::endl; \
        err = fd2(FUNC, comp, delta, tau, 1e-8, 1e-6, DAT, TOL, 1); \
        std::cout << "---" << std::endl; \
        if (!IGNORE_ERROR){ \
            return err; \
        } \
    } \
} \
stop = std::chrono::high_resolution_clock::now(); \
duration =  stop - start; \
std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

#define TEST_FUNCTION_OF_P_TAU(PROP, FUNC, DAT, TOL) \
start = std::chrono::high_resolution_clock::now(); \
std::cout << "    " << PROP << "(" << comp_str << ", P, tau) "; \
for(i=0; i<dat.size(); ++i){ \
    tau = pdat->T_star/dat[i][test_data::T_col]; \
    P = dat[i][test_data::P_col]*1000.0; \
    Psat = sat_p(comp, tau).f; \
    if((P >= Psat && FUNC == memo2_delta_liquid) || (P <= Psat && FUNC == memo2_delta_vapor)){ \
        if(fabs(P - pdat->Pc) < 0.1 && fabs(tau - tau_c(comp)) < 0.001) continue; \
        err = fd2(FUNC, comp, P, tau, 1e-3, 1e-8, DAT, TOL, 0); \
        if(err){ \
            std::cout << std::endl; \
            std::cout << "-------------------------ERROR--" << err << "------------------------" << std::endl; \
            std::cout << "Property: " << PROP << " Comp: " << comp_str << std::endl; \
            std::cout << "density " << dat[i][test_data::rho_col]; \
            std::cout << " tau: " << tau << " pressure: " << dat[i][test_data::P_col]*1000; \
            std::cout << ", T= " << dat[i][test_data::T_col] <<  std::endl; \
            fd2(FUNC, comp, P, tau, 1e-3, 1e-8, DAT, TOL, 1); \
            std::cout << "---" << std::endl; \
            return err; \
        } \
    } \
} \
stop = std::chrono::high_resolution_clock::now(); \
duration =  stop - start; \
std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

#define TEST_FUNCTION_OF_STATE_VARS(PROP, STATE1, FUNC, DATP, DATS1, TOL) \
start = std::chrono::high_resolution_clock::now(); \
std::cout << "    " << PROP << "(" << comp_str << ", " << STATE1 << ", P) "; \
for(i=0; i<dat.size(); ++i){ \
    P = dat[i][test_data::P_col]*1000.0; \
    if (fabs(P - pdat->Pc) < 1000){ \
        continue; \
    } \
    err = fd2(FUNC, comp, DATS1, P, 1e-3, 1e-8, DATP, TOL, 0); \
    if(err){ \
        std::cout << std::endl; \
        std::cout << "-------------------------ERROR--" << err << "------------------------" << std::endl; \
        std::cout << "Property: " << PROP << " Comp: " << comp_str << std::endl; \
        std::cout << "density " << dat[i][test_data::rho_col] << " delta: " << delta; \
        std::cout << " tau: " << tau << " pressure: " << dat[i][test_data::P_col]*1000; \
        std::cout << ", T= " << dat[i][test_data::T_col] <<  std::endl; \
        fd2(FUNC, comp, DATS1, P, 1e-3, 1e-8, DATP, TOL, 1); \
        std::cout << "---" << std::endl; \
        if (0){ \
            return err; \
        } \
    } \
} \
stop = std::chrono::high_resolution_clock::now(); \
duration =  stop - start; \
std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;

#define TEST_FUNCTION_OF_TP(PROP, STATE1, FUNC, DATP, DATS1, TOL) \
start = std::chrono::high_resolution_clock::now(); \
std::cout << "    " << PROP << "(" << comp_str << ", " << STATE1 << ", P) "; \
for(i=0; i<dat.size(); ++i){ \
    P = dat[i][test_data::P_col]*1000.0; \
    if (fabs(P - pdat->Pc) < 1000){ \
        continue; \
    } \
    err = fd2(FUNC, comp, DATS1, P, 1e-2, 1e-2, DATP, TOL, 0); \
    if(err){ \
        std::cout << std::endl; \
        std::cout << "-------------------------ERROR--" << err << "------------------------" << std::endl; \
        std::cout << "Property: " << PROP << " Comp: " << comp_str << std::endl; \
        std::cout << "density " << dat[i][test_data::rho_col] << " delta: " << delta; \
        std::cout << " tau: " << tau << " pressure: " << dat[i][test_data::P_col]*1000; \
        std::cout << ", T= " << dat[i][test_data::T_col] <<  std::endl; \
        fd2(FUNC, comp, DATS1, P, 1e-3, 1e-2, DATP, TOL, 1); \
        std::cout << "---" << std::endl; \
        if (0){ \
            return err; \
        } \
    } \
} \
stop = std::chrono::high_resolution_clock::now(); \
duration =  stop - start; \
std::cout << "Passed " << 5*dat.size() << " points in " << duration.count() << "s" << std::endl;


#define TEST_FUNCTION_SAT(PROP, STATE, FUNC, DATP, DATS, TOL, HSTEP) \
start = std::chrono::high_resolution_clock::now(); \
std::cout << "    " << PROP << "(" << comp_str << ", " << STATE << ") "; \
for(i=0; i<sat_liq_data.size(); ++i){ \
    tau = pdat->T_star/sat_liq_data[i][test_data::T_col]; \
    pressure = sat_liq_data[i][test_data::P_col]*1000; \
    err = fd1(FUNC, comp, DATS, HSTEP, DATP, TOL, 0); \
    if(err){ \
        std::cout << err; \
        fd1(FUNC, comp, DATS, HSTEP, DATP, TOL, 1); \
        return err; \
    } \
} \
stop = std::chrono::high_resolution_clock::now(); \
duration =  stop - start; \
std::cout << "Passed " << 3*sat_liq_data.size() << " points in " << duration.count() << "s" << std::endl;

#endif
#include"test_util.h"
#include"config.h"
#include"phi.h"
#include"props.h"
#include"props_hp.h"
#include"sat.h"
#include"delta.h"
#include"read_params.h"
#include"read_data.h"
#include"state.h"
#include<iostream>
#include<chrono>
#include<string>


int main(){
    uint comp, err;
    f22_struct res;
    f12_struct res1;
    std::vector<tests_struct> tests;
    tests = read_run_tests();

    // Verify some points
    comp = read_params("r32");
    parameters_struct *pdat = &cdata[comp];
    double rho=1055.258, T=273.15;
    std::cout << "r32 properties for rho = " << rho << " T = " << T << std::endl;
    res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " p = " << res.f << std::endl;
    res = memo2_enthalpy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " h = " << res.f << std::endl;
    res = memo2_entropy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_isochoric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cv = " << res.f << std::endl;
    res = memo2_isobaric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cp = " << res.f << std::endl;
    rho=22.09, T=273.15;
    std::cout << "r32 properties for rho = " << rho << " T = " << T << std::endl;
    res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " p = " << res.f << std::endl;
    res = memo2_enthalpy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " h = " << res.f << std::endl;
    res = memo2_entropy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_isochoric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cv = " << res.f << std::endl;
    res = memo2_isobaric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cp = " << res.f << std::endl; 

    comp = read_params("r125");
    pdat = &cdata[comp];
    rho=1319.818, T=273.15;
    std::cout << "r125 properties for rho = " << rho << " T = " << T << std::endl;
    res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " p = " << res.f << std::endl;
    res = memo2_enthalpy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " h = " << res.f << std::endl;
    res = memo2_entropy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_isochoric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cv = " << res.f << std::endl;
    res = memo2_isobaric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cp = " << res.f << std::endl;
    rho=42.070, T=273.15;
    std::cout << "r125 properties for rho = " << rho << " T = " << T << std::endl;
    res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " p = " << res.f << std::endl;
    res = memo2_enthalpy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " h = " << res.f << std::endl;
    res = memo2_entropy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_isochoric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cv = " << res.f << std::endl;
    res = memo2_isobaric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cp = " << res.f << std::endl;

    /*
    comp = read_params("r410a");
    pdat = &cdata[comp];
    rho=26.384, T=273.15+82.647;
    std::cout << "r410a properties for rho = " << rho << " T = " << T << std::endl;
    res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " p = " << res.f << std::endl;
    res = memo2_enthalpy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " h = " << res.f << std::endl;
    res = memo2_entropy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_isochoric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cv = " << res.f << std::endl;
    res = memo2_isobaric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cp = " << res.f << std::endl;
    */

    comp = read_params("r134a");
    pdat = &cdata[comp];
    //double rho=915.15, T=360.00; // P = pc
    //rho=1516.8, T=200.00; // P = pc
    rho=1417.7, T=233.15; // P = pc
    set_reference_state_offset(comp, 9.763426918083415, -4.858535799600055);
    //double rho=143.64, T=440.00; // P = pc
    // double rho=33.034, T=400.00; // P = 1 MPa
    // double rho=1512.1, T=200.00;
    // double rho=458.75, T=400.00; // P = 6 MPa
    // double rho=1519.8, T=200.00;
    // double rho=1022.3, T=350.00;
    std::cout << "r134a properties" << std::endl;
    res = memo2_entropy_hp(comp, 221.21, 4200);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_temperature_hp(comp, 221.21, 4200);
    std::cout << " s = " << res.f << std::endl;
    //return 0;
    std::cout << "r134a properties for rho = " << rho << " T = " << T << std::endl;
    res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " p = " << res.f << std::endl;
    res = memo2_enthalpy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " h = " << res.f << std::endl;
    res = memo2_entropy(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " s = " << res.f << std::endl;
    res = memo2_isobaric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cp = " << res.f << std::endl;
    res = memo2_isochoric_heat_capacity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " cv = " << res.f << std::endl;
    res = memo2_speed_of_sound(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " w = " << res.f << std::endl;
    res = memo2_isothermal_compressibility(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " isothermal compressibility = " << res.f << std::endl;
    res = memo2_viscosity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " mu = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " lambda = " << res.f << std::endl;
    res = memo2_surface_tension(comp, rho/pdat->rho_star, pdat->T_star/T);
    std::cout << " surface_tension = " << res.f << std::endl;
    set_reference_state_offset(comp, 0.0, 0.0);

    T=273.15 - 40.0;
    res1 = sat_delta_l(comp, pdat->T_star/T);
    std::cout.precision(12);
    std::cout << "r134a delta and tau for sat liquid at -40 C " << res1.f << " " << pdat->T_star/T << std::endl;

    comp = read_params("r1234ze");
    pdat = &cdata[comp];
    T=273.15 - 40.0;
    res1 = sat_delta_l(comp, pdat->T_star/T);
    std::cout << "r1234ze(e) delta and tau for sat liquid at -40 C " << res1.f << " " << pdat->T_star/T << std::endl;


    comp = read_params("co2");
    pdat = &cdata[comp];
    res = memo2_viscosity(comp, 2.440/pdat->rho_star, pdat->T_star/220.0);
    std::cout << " co2 mu = " << res.f << std::endl;
    res = memo2_viscosity(comp, 1.773/pdat->rho_star, pdat->T_star/300.0);
    std::cout << " co2 mu = " << res.f << std::endl;
    res = memo2_viscosity(comp, 1.773/pdat->rho_star, pdat->T_star/300.0);
    std::cout << " co2 mu = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 1.773/pdat->rho_star, pdat->T_star/300.0);
    std::cout << " co2 lambda = " << res.f << std::endl;
    res = memo2_surface_tension(comp, 1.773/pdat->rho_star, pdat->T_star/300.0);
    std::cout << " co2 sigma = " << res.f << std::endl;

    comp = read_params("h2o");
    pdat = &cdata[comp];
    res = memo2_thermal_conductivity(comp, 1/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 122/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 272/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 322/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 372/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 422/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_thermal_conductivity(comp, 750/pdat->rho_star, pdat->T_star/647.35);
    std::cout << " h2o lambda = " << res.f << std::endl;
    res = memo2_isothermal_compressibility(comp, 995.6/pdat->rho_star, pdat->T_star/303.15);
    std::cout << " h2o beta_T= " << res.f << std::endl;

    for (auto t = tests.begin(); t != tests.end(); ++t){
        std::cout << t->comp_str << std::endl;
        comp = read_params(t->comp_str);
        std::cout << "Calculated critical pressure for " << t->comp_str << ": " << pressure(comp, delta_c(comp), tau_c(comp)) << std::endl;
        if (t->test_set == "all"){
            err = run_set_all(comp, t->comp_str);
        }
        if (t->test_set == "mixed"){
            err = run_set_mixed(comp, t->comp_str);
        }
        if (err){
            exit(err);
        }
    }
    exit(0);
}
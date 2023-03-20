#include"test_util.h"
#include"config.h"
#include"phi.h"
#include"props.h"
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
    std::vector<tests_struct> tests;
    tests = read_run_tests();

    comp = read_params("r134a");
    parameters_struct *pdat = &cdata[comp];
    //double rho=915.15, T=360.00; // P = pc
    double rho=1516.8, T=200.00; // P = pc
    //double rho=143.64, T=440.00; // P = pc
    // double rho=33.034, T=400.00; // P = 1 MPa
    // double rho=1512.1, T=200.00;
    // double rho=458.75, T=400.00; // P = 6 MPa
    // double rho=1519.8, T=200.00;
    // double rho=1022.3, T=350.00;
    std::cout << "r134a properties for rho = " << rho << " T = " << T << std::endl;
    f22_struct res = memo2_pressure(comp, rho/pdat->rho_star, pdat->T_star/T);
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

    //return 0;
    for (auto t = tests.begin(); t != tests.end(); ++t){
        std::cout << t->comp_str << std::endl;
        comp = read_params(t->comp_str);
        std::cout << "Calculated critical pressure for " << t->comp_str << ": " << pressure(comp, delta_c(comp), tau_c(comp)) << std::endl;
        if (t->test_set == "all"){
            err = run_set_all(comp, t->comp_str, t->u_off, t->h_off, t->s_off);
        }
        if (t->test_set == "mixed"){
            err = run_set_mixed(comp, t->comp_str, t->u_off, t->h_off, t->s_off);
        }
        if (err){
            exit(err);
        }
    }
    exit(0);
}
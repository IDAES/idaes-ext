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
    std::string comp_str("co2");
    uint comp = read_params(comp_str);
    int err;

    
    double u_off = 506.778;
    double h_off = 506.778;
    double s_off = 2.738255753;
    
    /*
    double u_off = 0;
    double h_off = 0;
    double s_off = 0;
    */

    std::cout << std::endl << std::endl;

    std::cout << std::endl;
    std::cout << "Test basic comp liquid properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::liquid_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test basic comp vapor properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::vapor_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test basic comp supercritical properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(comp, comp_str, test_data::supercritical_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp sat. curve" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_sat_curve(comp, comp_str, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp liquid delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(comp, comp_str, test_data::liquid_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp vapor delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(comp, comp_str, test_data::vapor_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp supercritical delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(comp, comp_str, test_data::supercritical_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp functions for state var change on liquid data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(comp, comp_str, test_data::liquid_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp functions for state var change on vapor data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(comp, comp_str, test_data::vapor_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp functions for state var change on sc data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(comp, comp_str, test_data::supercritical_set, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test comp two phase" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_sat_curve_more(comp, comp_str, u_off, h_off, s_off);
    if(err){
        exit(err);
    }

    std::cout << std::endl << std::endl;
    exit(0);
}
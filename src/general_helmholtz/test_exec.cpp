#include"test_util.h"
#include"config.h"
#include"phi.h"
#include"props.h"
#include"sat.h"
#include"delta.h"
#include"read_params.h"
#include"read_data.h"
#include<iostream>
#include<chrono>

int main(){
    uint h2o = read_params("H2O");
    int err;

    std::cout << std::endl << std::endl;

    std::cout << std::endl;
    std::cout << "Test basic h2o liquid properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(h2o, "h2o", test_data::liquid_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test basic h2o supercritical properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(h2o, "h2o", test_data::supercritical_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test basic h2o vapor properties" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_basic_properties(h2o, "h2o", test_data::vapor_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o sat. curve" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_sat_curve(h2o, "h2o");
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o liquid delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(h2o, "h2o", test_data::liquid_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o vapor delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(h2o, "h2o", test_data::vapor_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o supercritical delta" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_delta_function(h2o, "h2o", test_data::supercritical_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o functions for state var change on liquid data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(h2o, "h2o", test_data::liquid_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o functions for state var change on vapor data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(h2o, "h2o", test_data::vapor_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o functions for state var change on sc data" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_state(h2o, "h2o", test_data::supercritical_set);
    if(err){
        exit(err);
    }

    std::cout << std::endl;
    std::cout << "Test h2o two phase" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    err = test_sat_curve_more(h2o, "h2o");
    if(err){
        exit(err);
    }


    std::cout << std::endl << std::endl;
    exit(0);
}
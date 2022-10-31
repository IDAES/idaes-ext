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

    for (auto t = tests.begin(); t != tests.end(); ++t){
        std::cout << t->comp_str << std::endl;
        uint comp = read_params(t->comp_str);
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
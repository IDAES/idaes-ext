#include"config.h"
#include<vector>

#ifndef _INCLUDE_PHI_H_
#define _INCLUDE_PHI_H_

std::vector<double> *phi_ideal(comp_enum comp, double delta, double tau);
std::vector<double> *phi_real(comp_enum comp, double delta, double tau);
void phi_real_for_sat(comp_enum comp, double delta, double tau, std::vector<double> *yvec_ptr);

#endif

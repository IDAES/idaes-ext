#include"config.h"
#include<vector>

#ifndef _INCLUDE_PROPS_H_
#define _INCLUDE_PROPS_H_

double pressure(comp_enum comp, double delta, double tau);
void pressure1(comp_enum comp, double delta, double tau, std::vector<double> *out);
void pressure2(comp_enum comp, double delta, double tau, std::vector<double> *out);

#endif

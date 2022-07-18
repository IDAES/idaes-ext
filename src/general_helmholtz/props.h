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
#include<vector>

#ifndef _INCLUDE_PROPS_H_
#define _INCLUDE_PROPS_H_

double pressure(comp_enum comp, double delta, double tau);
double entropy(comp_enum comp, double delta, double tau);
double enthalpy(comp_enum comp, double delta, double tau);
double internal_energy(comp_enum comp, double delta, double tau);

void pressure1(comp_enum comp, double delta, double tau, std::vector<double> *out);

void pressure2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void internal_energy2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void entropy2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void enthalpy2(comp_enum comp, double delta, double tau, std::vector<double> *out);


std::vector<double> *memo2_pressure(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_internal_energy(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_entropy(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_enthalpy(comp_enum comp, double delta, double tau);


#endif

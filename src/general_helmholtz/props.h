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

/*------------------------------------------------------------------------------
Author: John Eslick
File props.h

 This file contains basic property calculations as a function of delta and tau
 where delta = rho/rho_c and tau = T/T_c. The memo2_{prop} functions also
 memoize the property calculations with first and second derivatives.
------------------------------------------------------------------------------*/

#include"config.h"
#include<vector>
#include<unordered_map>

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
void gibbs2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void helmholtz2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void isochoric_heat_capacity2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void isobaric_heat_capacity2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void speed_of_sound2(comp_enum comp, double delta, double tau, std::vector<double> *out);

// Rather than calculate all the properties here, provide the terms needed to caluclate
void phi_ideal2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_ideal_d2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_ideal_t2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_ideal_dd2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_ideal_dt2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_ideal_tt2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_resi2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_resi_d2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_resi_t2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_resi_dd2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_resi_dt2(comp_enum comp, double delta, double tau, std::vector<double> *out);
void phi_resi_tt2(comp_enum comp, double delta, double tau, std::vector<double> *out);

std::vector<double> *memo2_pressure(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_internal_energy(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_entropy(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_enthalpy(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_gibbs(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_helmholtz(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_isochoric_heat_capacity(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_isobaric_heat_capacity(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_speed_of_sound(comp_enum comp, double delta, double tau);

// Rather than calculate all the properties here, provide the terms needed to caluclate
std::vector<double> *memo2_phi_ideal(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_ideal_d(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_ideal_t(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_ideal_dd(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_ideal_dt(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_ideal_tt(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_resi(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_resi_d(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_resi_t(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_resi_dd(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_resi_dt(comp_enum comp, double delta, double tau);
std::vector<double> *memo2_phi_resi_tt(comp_enum comp, double delta, double tau);

#endif

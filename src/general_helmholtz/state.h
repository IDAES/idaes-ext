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
#include "config.h"

#ifndef _INCLUDE_STATE_H_
#define _INCLUDE_STATE_H_

void enthalpy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out);
void entropy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out);
void internal_energy_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out);
void enthalpy_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out);
void entropy_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out);
void internal_energy_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out);

std::vector<double> *memo2_enthalpy_vapor(comp_enum comp, double pr, double tau);
std::vector<double> *memo2_entropy_vapor(comp_enum comp, double pr, double tau);
std::vector<double> *memo2_internal_energy_vapor(comp_enum comp, double pr, double tau);
std::vector<double> *memo2_enthalpy_liquid(comp_enum comp, double pr, double tau);
std::vector<double> *memo2_entropy_liquid(comp_enum comp, double pr, double tau);
std::vector<double> *memo2_internal_energy_liquid(comp_enum comp, double pr, double tau);

void tau_hp2(comp_enum comp, double ht, double pr, std::vector<double> *out);
void tau_sp2(comp_enum comp, double st, double pr, std::vector<double> *out);
void tau_up2(comp_enum comp, double ut, double pr, std::vector<double> *out);

std::vector<double> *memo2_tau_hp(comp_enum comp, double ht, double pr);
std::vector<double> *memo2_tau_sp(comp_enum comp, double st, double pr);
std::vector<double> *memo2_tau_up(comp_enum comp, double ut, double pr);

void vf_hp2(comp_enum comp, double ht, double pr, std::vector<double> *out);
void vf_sp2(comp_enum comp, double ht, double pr, std::vector<double> *out);
void vf_up2(comp_enum comp, double ht, double pr, std::vector<double> *out);

std::vector<double> *memo2_vf_hp(comp_enum comp, double ht, double pr);
std::vector<double> *memo2_vf_sp(comp_enum comp, double ht, double pr);
std::vector<double> *memo2_vf_up(comp_enum comp, double ht, double pr);

#endif

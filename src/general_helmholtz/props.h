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

double pressure(uint comp, double delta, double tau);
double entropy(uint comp, double delta, double tau);
double enthalpy(uint comp, double delta, double tau);
double internal_energy(uint comp, double delta, double tau);

void pressure1(uint comp, double delta, double tau, f21_struct *out);

void pressure2(uint comp, double delta, double tau, f22_struct *out);
void internal_energy2(uint comp, double delta, double tau, f22_struct *out);
void entropy2(uint comp, double delta, double tau, f22_struct *out);
void enthalpy2(uint comp, double delta, double tau, f22_struct *out);
void gibbs2(uint comp, double delta, double tau, f22_struct *out);
void helmholtz2(uint comp, double delta, double tau, f22_struct *out);
void isochoric_heat_capacity2(uint comp, double delta, double tau, f22_struct *out);
void isobaric_heat_capacity2(uint comp, double delta, double tau, f22_struct *out);
void speed_of_sound2(uint comp, double delta, double tau, f22_struct *out);
void specific_volume2(uint comp, double delta, double tau, f22_struct *out);
void isothermal_compressibility2(uint comp, double delta, double tau, f22_struct *out);

// Rather than calculate all the properties here, provide the terms needed to caluclate
void phi_ideal2(uint comp, double delta, double tau, f22_struct *out);
void phi_ideal_d2(uint comp, double delta, double tau, f22_struct *out);
void phi_ideal_t2(uint comp, double delta, double tau, f22_struct *out);
void phi_ideal_dd2(uint comp, double delta, double tau, f22_struct *out);
void phi_ideal_dt2(uint comp, double delta, double tau, f22_struct *out);
void phi_ideal_tt2(uint comp, double delta, double tau, f22_struct *out);
void phi_resi2(uint comp, double delta, double tau, f22_struct *out);
void phi_resi_d2(uint comp, double delta, double tau, f22_struct *out);
void phi_resi_t2(uint comp, double delta, double tau, f22_struct *out);
void phi_resi_dd2(uint comp, double delta, double tau, f22_struct *out);
void phi_resi_dt2(uint comp, double delta, double tau, f22_struct *out);
void phi_resi_tt2(uint comp, double delta, double tau, f22_struct *out);

f22_struct memo2_pressure(uint comp, double delta, double tau);
f22_struct memo2_internal_energy(uint comp, double delta, double tau);
f22_struct memo2_entropy(uint comp, double delta, double tau);
f22_struct memo2_enthalpy(uint comp, double delta, double tau);
f22_struct memo2_gibbs(uint comp, double delta, double tau);
f22_struct memo2_helmholtz(uint comp, double delta, double tau);
f22_struct memo2_isochoric_heat_capacity(uint comp, double delta, double tau);
f22_struct memo2_isobaric_heat_capacity(uint comp, double delta, double tau);
f22_struct memo2_speed_of_sound(uint comp, double delta, double tau);
f22_struct memo2_specific_volume(uint comp, double delta, double tau);
f22_struct memo2_isothermal_compressibility(uint comp, double delta, double tau);

// Rather than calculate all the properties here, provide the terms needed to calculate
f22_struct memo2_phi_ideal(uint comp, double delta, double tau);
f22_struct memo2_phi_ideal_d(uint comp, double delta, double tau);
f22_struct memo2_phi_ideal_t(uint comp, double delta, double tau);
f22_struct memo2_phi_ideal_dd(uint comp, double delta, double tau);
f22_struct memo2_phi_ideal_dt(uint comp, double delta, double tau);
f22_struct memo2_phi_ideal_tt(uint comp, double delta, double tau);
f22_struct memo2_phi_resi(uint comp, double delta, double tau);
f22_struct memo2_phi_resi_d(uint comp, double delta, double tau);
f22_struct memo2_phi_resi_t(uint comp, double delta, double tau);
f22_struct memo2_phi_resi_dd(uint comp, double delta, double tau);
f22_struct memo2_phi_resi_dt(uint comp, double delta, double tau);
f22_struct memo2_phi_resi_tt(uint comp, double delta, double tau);

#endif

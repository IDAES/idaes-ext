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
File state.cpp

Functions to enable change of state variables. The end goal of this section is
given a set of state varible, calculate T, P, and vapor fraction.  From there
you can calculate delta for each phase and then all the rest of the properties
can be calculated.

The general method for changing state varaibles is (first step is here):
  1) tau = tau(v1, v2), vf=vf(v1, v2), p(v1, v2); so far p is always a state
     variable but it doesn't need to be that case we can add on to support more
  2) delta_v = delta_v(p, tau), delta_l = delta_l(p, tau)
  3) property_v = f(delta_v, tau), property_l = f(delta_l, tau)
  4) calculate mixed phase properies
------------------------------------------------------------------------------*/


#include "config.h"

#ifndef _INCLUDE_STATE_H_
#define _INCLUDE_STATE_H_

void enthalpy_vapor2(uint comp, double pr, double tau, f22_struct *out);
void entropy_vapor2(uint comp, double pr, double tau, f22_struct *out);
void internal_energy_vapor2(uint comp, double pr, double tau, f22_struct *out);
void enthalpy_liquid2(uint comp, double pr, double tau, f22_struct *out);
void entropy_liquid2(uint comp, double pr, double tau, f22_struct *out);
void internal_energy_liquid2(uint comp, double pr, double tau, f22_struct *out);

f22_struct memo2_enthalpy_vapor(uint comp, double pr, double tau);
f22_struct memo2_entropy_vapor(uint comp, double pr, double tau);
f22_struct memo2_internal_energy_vapor(uint comp, double pr, double tau);
f22_struct memo2_enthalpy_liquid(uint comp, double pr, double tau);
f22_struct memo2_entropy_liquid(uint comp, double pr, double tau);
f22_struct memo2_internal_energy_liquid(uint comp, double pr, double tau);

void tau_hp2(uint comp, double ht, double pr, f22_struct *out);
void tau_sp2(uint comp, double st, double pr, f22_struct *out);
void tau_up2(uint comp, double ut, double pr, f22_struct *out);

f22_struct memo2_tau_hp(uint comp, double ht, double pr);
f22_struct memo2_tau_sp(uint comp, double st, double pr);
f22_struct memo2_tau_up(uint comp, double ut, double pr);

void vf_hp2(uint comp, double ht, double pr, f22_struct *out);
void vf_sp2(uint comp, double ht, double pr, f22_struct *out);
void vf_up2(uint comp, double ht, double pr, f22_struct *out);

f22_struct memo2_vf_hp(uint comp, double ht, double pr);
f22_struct memo2_vf_sp(uint comp, double ht, double pr);
f22_struct memo2_vf_up(uint comp, double ht, double pr);

#endif

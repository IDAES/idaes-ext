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

/*--------------------------------------------------------------------------------
 Calaculate components of dimensionless ideal and residual Helmholtz free energy
 and derivatives to fourth order (2 for thermo + 2 for optimization solver). The
 file also contains function calls for aux curves (approx sat densities).

 Author: John Eslick
 File: phi.h
--------------------------------------------------------------------------------*/
#include"config.h"

#ifndef _INCLUDE_PHI_H_
#define _INCLUDE_PHI_H_

// Three derivative versions (all we need for usual properties
// p, u, s, h, f, g)
f23_struct phi_ideal(uint comp, double delta, double tau);
f23_struct phi_resi(uint comp, double delta, double tau);

// Four derivative versions (cp, cv, speed of sound ...) only use when needed.
f24_struct phi_ideal4(uint comp, double delta, double tau);
f24_struct phi_resi4(uint comp, double delta, double tau);

// For sat functions only need derivatives wrt delta and no cache
f12_struct phi_resi_for_sat(uint comp, double delta, double tau);

// Approx sat densities
double sat_delta_v_approx(uint comp, double tau);
double sat_delta_l_approx(uint comp, double tau);

#endif

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
 and derivatives to fourth order (2 for thermo + 2 for optimization solver).

 Author: John Eslick
 File: phi.h
--------------------------------------------------------------------------------*/

#include"config.h"
#include<vector>

#ifndef _INCLUDE_PHI_H_
#define _INCLUDE_PHI_H_

std::vector<double> *phi_ideal(comp_enum comp, double delta, double tau);
std::vector<double> *phi_resi(comp_enum comp, double delta, double tau);
void phi_resi_for_sat(comp_enum comp, double delta, double tau, std::vector<double> *yvec_ptr);

#endif

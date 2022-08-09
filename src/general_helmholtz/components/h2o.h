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
 Specific water functions from:

 International Association for the Properties of Water and Steam (2016).
     IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
     the Properties of Ordinary Water Substance for General Scientific Use,"
     URL: http://iapws.org/relguide/IAPWS95-2016.pdf
 Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
     Thermodynamic Properties of Ordinary Water Substance for General and
     Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.
 Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
     Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
     and Power, 122, 150-182.

 Author: John Eslick
 File: h2o.h
--------------------------------------------------------------------------------*/

#ifndef _INCLUDE_H2O_H_
#define _INCLUDE_H2O_H_

double delta_sat_l_approx_h2o(double tau);
double delta_sat_v_approx_h2o(double tau);
double melting_temperature_h2o(double pr);
double melting_liquid_density_h2o(double pr);

void phi_h2o_resi_tape();
void phi_h2o_ideal_tape();

#endif

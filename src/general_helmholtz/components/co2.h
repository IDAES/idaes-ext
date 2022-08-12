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

 Span, R., and W. Wanger (1996). "A New Equation of State for Carbon Dioxide
     Covering the Fluid Region from the Triple-Point Temperature to 1100 K as
     Pressures up to 800 MPa." Journal of Physical and Chemical Reference Data,
     25, 1509.

 Author: John Eslick
 File: co2.h
--------------------------------------------------------------------------------*/

#ifndef _INCLUDE_CO2_H_
#define _INCLUDE_CO2_H_

double delta_sat_l_approx_co2(double tau);
double delta_sat_v_approx_co2(double tau);
double melting_temperature_co2(double pr);
double melting_liquid_density_co2(double pr);

void phi_co2_resi_tape();
void phi_co2_ideal_tape();

#endif

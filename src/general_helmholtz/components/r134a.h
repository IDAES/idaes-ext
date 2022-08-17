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
 Specific r134a (r1234ze(e)) functions from:

Monika Thol and Eric W. Lemmon. "Equation of State for the Thermodynamic
  Properties of trans-1,3,3,3-Tetrafluoropropene [R-1234ze(E)]." Int. J.
  Thermophys, 37(3):1–16, 2016. doi:10.1007/s10765-016-2040-6.

 Author: John Eslick
 File: r134a.h
--------------------------------------------------------------------------------*/

#ifndef _INCLUDE_R134A_H_
#define _INCLUDE_R134A_H_

double delta_sat_l_approx_r134a(double tau);
double delta_sat_v_approx_r134a(double tau);
double melting_tau_r134a(double pr);
double melting_liquid_delta_r134a(double pr);

void phi_r134a_resi_tape();
void phi_r134a_ideal_tape();

#endif

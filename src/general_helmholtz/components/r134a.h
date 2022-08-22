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
Specific r134a functions from:

Perkins, R.A.; Laesecke, A.; Howley, J.; Ramires, M.L.V.; Gurova, A.N.; Cusco, L.,
    Experimental thermal conductivity values for the IUPAC round-robin sample of
    1,1,1,2-tetrafluoroethane (R134a), NIST Interagency/Internal Report (NISTIR)
    - 6605, 2000, https://doi.org/10.6028/NIST.IR.6605.

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

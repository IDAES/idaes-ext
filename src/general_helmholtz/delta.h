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
 Provide functions to solve for liquid and vapor density from T and P

 Author: John Eslick
 File: delta.h
--------------------------------------------------------------------------------*/

#ifndef _INCLUDE_DELTA_H_
#define _INCLUDE_DELTA_H_

double delta_liquid(uint comp, double pr, double tau);
double delta_vapor(uint comp, double pr, double tau);

void delta_liquid2(uint comp, double pr, double tau, f22_struct *out);
void delta_vapor2(uint comp, double pr, double tau, f22_struct *out);

f22_struct memo2_delta_liquid(uint comp, double pr, double tau);
f22_struct memo2_delta_vapor(uint comp, double pr, double tau);

#endif

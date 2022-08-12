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

struct pressure_wrap_state {
  comp_enum comp;
  double p;
  double tau;
};

double delta_liquid(comp_enum comp, double pr, double tau);
double delta_vapor(comp_enum comp, double pr, double tau);

void delta_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out);
void delta_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out);


std::vector<double> *memo2_delta_liquid(comp_enum comp, double pr, double tau);
std::vector<double> *memo2_delta_vapor(comp_enum comp, double pr, double tau);

#endif

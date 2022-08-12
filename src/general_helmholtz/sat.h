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
File sat.h

 This file contains the functions to solve the saturation curve.  The method
 is described in the following paper:

 Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
     State from Helmholtz Energy Equations of State." Journal of Thermal
     Science and Technology, 3(3), 442-451.
------------------------------------------------------------------------------*/

#ifndef _INCLUDE_SAT_H_
#define _INCLUDE_SAT_H_

std::vector<double> *sat_tau(comp_enum comp, double pr);
std::vector<double> *sat_p(comp_enum comp, double tau);
std::vector<double> *sat_delta_v(comp_enum comp, double tau);
std::vector<double> *sat_delta_l(comp_enum comp, double tau);

#endif

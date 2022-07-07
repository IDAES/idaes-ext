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

#include<vector>
#include"config.h"

#ifndef _INCLUDE_TESTING_H_
#define _INCLUDE_TESTING_H_

typedef std::vector<double>* (*test_fptr1)(comp_enum comp, double x);
typedef std::vector<double>* (*test_fptr2)(comp_enum comp, double x, double y);

#endif

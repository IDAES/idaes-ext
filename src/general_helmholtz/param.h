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

#ifndef _INCLUDE_PARAM_H_
#define _INCLUDE_PARAM_H_

// Specific gas constant (kJ/kg/K)
const double R[] = {
  0,          // not used
  0.46151805, // H2O
};

// Critical temperature (K)
const double Tc[] = {
  0,        // not used
  647.096,  // H2O
};

// Critical density (kg/m3)
const double rhoc[] = {
  0,        // not used
  322.0     // H2O
};

// Critical Pressure
const double Pc[] = {
  0,          // not used
  22064.0     // H2O
};

// Triple point temperature (K)
const double Tt[] = {
  0,        // not used
  273.16    // H2O
};

// Triple point pressure (kPa)
const double Pt[] = {
  0,          // not used
  0.611655    // H2O
};

// Triple point liquid density (kPa)
const double rhot_l[] = {
  0,          // not used
  999.793     // H2O
};

// Triple point vapor density (kPa)
const double rhot_v[] = {
  0,          // not used
  0.00485458  // H2O
};

// upper pressue bound (kPa)
const double P_max[] = {
  0,       // not used
  1e6,     // h2o
};

// upper density bound (kg/m3)
const double rho_max[] = {
  0,       // not used
  1200.0,  // h2o
};

// upper pressue bound (kPa)
const double T_min[] = {
  0,       // not used
  235,     // h2o
};

#endif

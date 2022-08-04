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

namespace param {
  // MW g/mol
  const double mw[] = {
    0,           // not used
    18.015268,   // H2O
    44.0098,     // CO2
    114.0416     // r1234ze
  };


  // Specific gas constant (kJ/kg/K)
  const double R[] = {
    0,           // not used
    0.46151805,  // H2O
    0.1889241,   // CO2
    0.07290727,  //r1234ze
  };

  // Critical temperature (K)
  const double Tc[] = {
    0,        // not used
    647.096,  // H2O
    304.1282, // CO2
    382.513,  // r1234ze
  };

  // Critical density (kg/m3)
  const double rhoc[] = {
    0,        // not used
    322.0,    // H2O
    467.6,    // CO2
    489.238,  // r1234ze
  };

  // Critical Pressure (kPa)
  const double Pc[] = {
    0,          // not used
    22064.0,    // H2O
    7377.3,     // CO2
    3634.9,     // r1234ze
  };

  // Triple point temperature (K)
  const double Tt[] = {
    0,        // not used
    273.16,   // H2O
    216.592,  // CO2
    169.0,    // r1234ze
  };

  // Triple point pressure (kPa)
  const double Pt[] = {
    0,          // not used
    0.611655,   // H2O
    517.95,     // CO2
    0.228564,   // r1234ze (estimated from Tt and sat P)
  };

  // Triple point liquid density (kPa)
  const double rhot_l[] = {
    0,          // not used
    999.793,    // H2O
    1178.46,    // CO2
    1510.76,    // r1234ze
  };

  // Triple point vapor density (kPa)
  const double rhot_v[] = {
    0,          // not used
    0.00485458, // H2O
    13.761,     // CO2
    0.0185566,  // r1234ze
  };

  // upper pressue bound (kPa)
  const double P_min[] = {
    0,       // not used
    1e-9,    // h2o
    1e-9,    // CO2
    1e-9,    // r1234ze
  };

  // upper pressue bound (kPa)
  const double P_max[] = {
    0,       // not used
    1e6,     // h2o
    1e6,     // co2
    5e5,     // r1234ze
  };

  // upper density bound (kg/m3)
  const double rho_max[] = {
    0,       // not used
    1250.0,  // h2o
    1500.0,  // CO2
    2500.0,  //r1234ze

  };

  // upper pressue bound (kPa)
  const double T_min[] = {
    0,       // not used
    235,     // h2o
    200,     // CO2
    150,     // r1234ze
  };

  // upper bound on T (K)
  const double T_max[] = {
    0,       // not used
    1300,    // h2o
    1300,    // CO2
    1300,    // r1234ze
  };

}

#endif

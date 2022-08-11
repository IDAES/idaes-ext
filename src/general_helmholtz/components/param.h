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

 Component parameters, these arrays need to be in order so that the indexes line
 up with the component enum

 Author: John Eslick
 File: param.h
--------------------------------------------------------------------------------*/

#include<unordered_map>
#include<string>

#ifndef _INCLUDE_PARAM_H_
#define _INCLUDE_PARAM_H_

#define NCOMPS 4  // make this one more than the last component index

// if adding components update NCOMPS above should be 1 more than last index
enum comp_enum{
  h2o = 1,
  co2 = 2,
  r1234ze = 3,
};

static std::unordered_map<std::string, comp_enum> comp_string_table = {
  {"h2o", comp_enum::h2o},
  {"co2", comp_enum::co2},
  {"r1234ze", comp_enum::r1234ze},
};

static std::unordered_map<uint, std::string> comp_enum_table = {
  {comp_enum::h2o, "h2o"},
  {comp_enum::co2, "co2"},
  {comp_enum::r1234ze, "r1234ze"},
};

namespace param {
  // MW g/mol
  const double mw[] = {
    0,           // not used
    18.015268,   // h2o
    44.0098,     // co2
    114.0416     // r1234ze
  };


  // Specific gas constant (kJ/kg/K)
  const double R[] = {
    0,           // not used
    0.46151805,  // h2o
    0.1889241,   // co2
    0.07290727,  //r1234ze
  };

  // Critical temperature (K)
  const double Tc[] = {
    0,        // not used
    647.096,  // h2o
    304.1282, // co2
    382.513,  // r1234ze
  };

  // Critical density (kg/m3)
  const double rhoc[] = {
    0,        // not used
    322.0,    // h2o
    467.6,    // co2
    489.238,  // r1234ze
  };

  // Critical Pressure (kPa)
  const double Pc[] = {
    0,          // not used
    22064.0,    // h2o
    7377.3,     // co2
    3634.9,     // r1234ze
  };

  // Triple point temperature (K)
  const double Tt[] = {
    0,        // not used
    273.16,   // h2o
    216.592,  // co2
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
    999.793,    // h2o
    1178.46,    // co2
    1510.76,    // r1234ze
  };

  // Triple point vapor density (kPa)
  const double rhot_v[] = {
    0,          // not used
    0.00485458, // h2o
    13.761,     // co2
    0.0185566,  // r1234ze
  };

  // upper pressue bound (kPa)
  const double P_min[] = {
    0,       // not used
    1e-9,    // h2o
    1e-9,    // co2
    1e-9,    // r1234ze
  };

  // upper pressue bound (kPa)
  const double P_max[] = {
    0,       // not used
    1e6,     // h2o
    5e5,     // co2
    5e5,     // r1234ze
  };

  // upper density bound (kg/m3)
  const double rho_max[] = {
    0,       // not used
    1250.0,  // h2o
    1500.0,  // co2
    2500.0,  //r1234ze

  };

  // upper pressue bound (kPa)
  const double T_min[] = {
    0,       // not used
    235,     // h2o
    200,     // co2
    150,     // r1234ze
  };

  // upper bound on T (K)
  const double T_max[] = {
    0,       // not used
    1300,    // h2o
    1300,    // co2
    1300,    // r1234ze
  };

}

#endif

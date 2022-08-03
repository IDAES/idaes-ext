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
#include <math.h>
#include <adolc/adolc.h>
#include "config.h"
#include "param.h"

double melting_temperature_r1234ze(double pr){
  /*
    Estimate the melting temperature at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid tange of
    temperatures at a given pressure (kPa).
  */
  return 200.0;
}

double melting_liquid_density_r1234ze(double pr){
  /*
    Estimate the melting liquid density at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid range of
    temperatures at a given pressure (kPa).
  */
  return 2000;
}

double delta_sat_v_approx_r1234ze(double tau){
/*
  Approximate saturated vapor density
  This equation is from the original IAPWS-95 paper
*/
  double XX = 1 - 1.0/tau;
  return exp(
    -1.0308*pow(XX, 0.24) +
    -5.0422*pow(XX, 0.72) +
    -11.5*pow(XX, 2.1) +
    -37.499*pow(XX, 4.8) +
    -77.945*pow(XX, 9.5)
  );
}

double delta_sat_l_approx_r1234ze(double tau){
/*
  Approximate saturated vapor liquid
  This equation is from the original IAPWS-95 paper
*/
  double XX = 1 - 1.0/tau;
  return
    1.00
    + 1.1913*pow(XX, 0.27)
    + 2.2456*pow(XX, 0.7)
    - 1.7747*pow(XX, 1.25)
    + 1.3096*pow(XX, 1.9);
}


void phi_r1234ze_ideal_tape(){
// Create a ADOL-C tape for the ideal part of phi for R1234ze
  double out;

  taped_ideal[comp_enum::r1234ze] = PHI_IDEAL_TAPE(comp_enum::r1234ze);
  trace_on(PHI_IDEAL_TAPE(comp_enum::r1234ze));
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  // ideal part
  y[0] = log(x[0])
    -12.558347537 +
    8.7912297624*x[1] +
    (4.0 - 1.0)*log(x[1]) +
    9.3575*log(1 - exp(-513.0*x[1]/param::Tc[comp_enum::r1234ze])) +
    10.717*log(1 - exp(-1972.0*x[1]/param::Tc[comp_enum::r1234ze]));
  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}


void phi_r1234ze_resi_tape(){
// Create a ADOL-C tape for the ideal part of phi for R1234ze

  // parameters are only need once locally to create the tape
  double N[] = {
    0,            // 0 (not used)
    0.03982797,   // 1
    1.812227,     // 2
   -2.537512,     // 3
   -0.5333254,    // 4
    0.1677031,    // 5
   -1.323801,     // 6
   -0.6694654,    // 7
    0.8072718,    // 8
   -0.7740229,    // 9
   -0.01843846,   // 10
    1.407916,     // 11
   -0.4237082,    // 12
   -0.2270068,    // 13
   -0.805213,     // 14
    0.00994318,   // 15
   -0.008798793,  // 16
  };

  double d[] = {
    0,   // 0 (not used)
    4,   // 1
    1,   // 2
    1,   // 3
    2,   // 4
    3,   // 5
    1,   // 6
    3,   // 7
    2,   // 8
    2,   // 9
    7,   // 10
    1,   // 11
    1,   // 12
    3,   // 13
    3,   // 14
    2,   // 15
    1,   // 16
  };

  double t[] = {
    0.0,    // 0 (not used)
    1.0,    // 1
    0.223,  // 2
    0.755,  // 3
    1.24,   // 4
    0.44,   // 5
    2.0,    // 6
    2.2,    // 7
    1.2,    // 8
    1.5,    // 9
    0.9,    // 10
    1.33,   // 11
    1.75,   // 12
    2.11,   // 13
    1.0,    // 14
    1.5,    // 15
    1.0,    // 16
  };

  double l[] = {
      0,  // 0 (not used)
      0,  // 1 (not used)
      0,  // 2 (not used)
      0,  // 3 (not used)
      0,  // 4 (not used)
      0,  // 5 (not used)
      2,  // 6
      2,  // 7
      1,  // 8
      2,  // 9
      1,  // 10
      0,  // 11 (not used)
      0,  // 12 (not used)
      0,  // 13 (not used)
      0,  // 14 (not used)
      0,  // 15 (not used)
      0,  // 16 (not used)
  };

  double eta[] = {
    0,    // 0 (not used)
    0,    // 1 (not used)
    0,    // 2 (not used)
    0,    // 3 (not used)
    0,    // 4 (not used)
    0,    // 5 (not used)
    0,    // 6 (not used)
    0,    // 7 (not used)
    0,    // 8 (not used)
    0,    // 9 (not used)
    0,    // 10 (not used)
    1.0,  // 11
    1.61, // 12
    1.24, // 13
    9.34, // 14
    5.78, // 15
    3.08, // 16
  };

  double b[] = {
    0,    // 0 (not used)
    0,    // 1 (not used)
    0,    // 2 (not used)
    0,    // 3 (not used)
    0,    // 4 (not used)
    0,    // 5 (not used)
    0,    // 6 (not used)
    0,    // 7 (not used)
    0,    // 8 (not used)
    0,    // 9 (not used)
    0,    // 10 (not used)
    1.21, // 11
    1.37, // 12
    0.98, // 13
  171.0,  // 14
   47.4,  // 15
   15.4,  // 16
  };

  double g[] = {
    0,      // 0 (not used)
    0,      // 1 (not used)
    0,      // 2 (not used)
    0,      // 3 (not used)
    0,      // 4 (not used)
    0,      // 5 (not used)
    0,      // 6 (not used)
    0,      // 7 (not used)
    0,      // 8 (not used)
    0,      // 9 (not used)
    0,      // 10 (not used)
    0.943,  // 11 (not used)
    0.642,  // 12 (not used)
    0.59,   // 13 (not used)
    1.2,    // 14 (not used)
    1.33,   // 15 (not used)
    0.64,   // 16 (not used)
  };

  double e[] = {
    0,      // 0 (not used)
    0,      // 1 (not used)
    0,      // 2 (not used)
    0,      // 3 (not used)
    0,      // 4 (not used)
    0,      // 5 (not used)
    0,      // 6 (not used)
    0,      // 7 (not used)
    0,      // 8 (not used)
    0,      // 9 (not used)
    0,      // 10 (not used)
    0.728,  // 11
    0.87,   // 12
    0.855,  // 13
    0.79,   // 14
    1.3,    // 15
    0.71,   // 16
  };
  int i = 0;
  double out;

  taped_resi[comp_enum::r1234ze] = PHI_RESI_TAPE(comp_enum::r1234ze);
  trace_on(PHI_RESI_TAPE(comp_enum::r1234ze));
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  y[0] = 0;
  for(i=1; i<=5; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]);
  }
  for(i=6; i<=10; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-pow(x[0], l[i]));
  }
  for(i=11; i<=16; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]) *
      exp(
        -eta[i]*(x[0] - e[i])*(x[0] - e[i]) -
         b[i]*(x[1] - g[i])*(x[1] - g[i])
      );
  }

  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}

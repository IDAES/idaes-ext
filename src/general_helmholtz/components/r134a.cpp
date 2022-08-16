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

Tillner-Roth, R.; Baehr, H.D., An International Standard Formulation for the
    Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane (HFC-134a) for
    Temperatures from 170 K to 455 K and Pressures up to 70 MPa, J. Phys. Chem.
    Ref. Data, 1994, 23, 5, 657-729, https://doi.org/10.1063/1.555958

 Author: John Eslick
 File: r134a.cpp
--------------------------------------------------------------------------------*/

#include <math.h>
#include <adolc/adolc.h>
#include "../config.h"

double melting_tau_r134a(double pr){
  /*
    Estimate the melting temperature at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid tange of
    temperatures at a given pressure (kPa). If there is no good metling curve,
    data just supply a resonable upper limit on vapor temperature.
  */
  return param::T_star[r134a]/170.0;
}

double melting_liquid_delta_r134a(double pr){
  /*
    Estimate the melting liquid density at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid range of
    temperatures at a given pressure (kPa). If there is no good metling curve,
    data just supply a resonable upper limit on vapor density.
  */
  return 1500.0/param::rho_star[r134a];
}

double delta_sat_v_approx_r134a(double tau){
/*
  Approximate saturated vapor density
*/
  double XX = 1 - 1.0/tau;
  return exp(
    -2.837294*pow(XX, 1.0/3.0) +
    -7.875988*pow(XX, 2.0/3.0) +
    4.478586*pow(XX, 1.0/2.0) +
    -14.140125*pow(XX, 9.0/4.0) +
    -52.361297*pow(XX, 11.0/12.0)
  );
}

double delta_sat_l_approx_r134a(double tau){
/*
  Approximate saturated vapor liquid
*/
  double XX = 1 - 1.0/tau;
  return
    518.20
    + 884.13*pow(XX, 1.0/3.0)
    + 485.84*pow(XX, 2.0/3.0)
    + 193.29*pow(XX, 10.0/3.0);
}


void phi_r134a_ideal_tape(){
// Create a ADOL-C tape for the ideal part of phi for R1234ze
  double out;

  taped_ideal[comp_enum::r134a] = PHI_IDEAL_TAPE(comp_enum::r134a);
  trace_on(PHI_IDEAL_TAPE(comp_enum::r134a));
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  // ideal part
  y[0] = log(x[0])
    -1.019535 +
    9.047135*x[1] +
    -1.629789*log(x[1]) +
    -9.723916*pow(x[1], -1.0/2.0) +
    -3.927170*pow(x[1], -3.0/4.0);
  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}


void phi_r134a_resi_tape(){
// Create a ADOL-C tape for the residual part of phi for R1234ze

  // parameters are only need once locally to create the tape
  double N[] = {
    0,              // 0 (not used)
    0.5586817e-1,   // 1
    0.4982230e+0,   // 2
    0.2458698e-1,   // 3
    0.8570145e-3,   // 4
    0.4788584e-3,   // 5
   -0.1800808e+1,   // 6
    0.2671641e+0,   // 7
   -0.4781652e-1,   // 8
    0.1423987e-1,   // 9
    0.3324062e+0,   // 10
   -0.7485907e-2,   // 11
    0.1017263e-3,   // 12
   -0.5184567e+0,   // 13
   -0.8692288e-1,   // 14
    0.2057144e+0,   // 15
   -0.5000457e-2,   // 16
    0.4603262e-3,   // 17
   -0.3497836e-2,   // 18
    0.6995038e-2,  // 19
   -0.1452184e-1,  // 20
   -0.1285458e-3,  // 21
  };

  double d[] = {
    0,   // 0 (not used)
    2,   // 1
    1,   // 2
    3,   // 3
    6,   // 4
    6,   // 5
    1,   // 6
    1,   // 7
    2,   // 8
    5,   // 9
    2,   // 10
    2,   // 11
    4,   // 12
    1,   // 13
    4,   // 14
    1,   // 15
    2,   // 16
    4,   // 17
    1,   // 18
    5,   // 19
    3,   // 20
    10,  // 21
  };

  double t[] = {
    0.0,       // 0 (not used)
   -1.0/2.0,   // 1
    0.0,       // 2
    0.0,       // 3
    0.0,       // 4
    3.0/2.0,   // 5
    3.0/2.0,   // 6
    2.0,       // 7
    2.0,       // 8
    1.0,       // 9
    3.0,       // 10
    5.0,       // 11
    1.0,       // 12
    5.0,       // 13
    5.0,       // 14
    6.0,       // 15
    10.0,      // 16
    10.0,      // 17
    10.0,      // 18
    18.0,      // 19
    22.0,      // 20
    50.0,      // 21
  };

  int i = 0;
  double out;

  taped_resi[comp_enum::r134a] = PHI_RESI_TAPE(comp_enum::r134a);
  trace_on(PHI_RESI_TAPE(comp_enum::r134a));
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  y[0] = 0;
  for(i=1; i<=8; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]);
  }
  for(i=9; i<=11; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-x[0]);
  }
  for(i=12; i<=17; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-x[0]*x[0]);
  }
  for(i=18; i<=20; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-x[0]*x[0]*x[0]);
  }
  for(i=21; i<=21; ++i){
    y[0] += N[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-x[0]*x[0]*x[0]*x[0]);
  }

  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}

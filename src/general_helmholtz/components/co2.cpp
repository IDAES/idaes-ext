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

 Span, R., and W. Wanger (1996). "A New Equation of State for Carbon Dioxide
     Covering the Fluid Region from the Triple-Point Temperature to 1100 K as
     Pressures up to 800 MPa." Journal of Physical and Chemical Reference Data,
     25, 1509.

 Author: John Eslick
 File: co2.cpp
--------------------------------------------------------------------------------*/

#include <math.h>
#include <adolc/adolc.h>
#include "../config.h"

double melting_temperature_co2(double pr){
  /*
    Estimate the melting temperature at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid tange of
    temperatures at a given pressure (kPa). If there is no good metling curve,
    data just supply a resonable upper limit on vapor temperature.
  */
  if(pr < param::Pt[co2]){ // melting
    return -0.000000*pr*pr + 0.000201*pr + 216.687349;
  }
  // Sublimation
  return 150.56*pow(pr, 0.0559);
}

double melting_liquid_density_co2(double pr){
  /*
    Estimate the melting liquid density at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid range of
    temperatures at a given pressure (kPa). If there is no good metling curve,
    data just supply a resonable upper limit on vapor density.
  */
  return -3.776E-10*pr*pr + 6.860E-04*pr + 1.187E+03;
}

double delta_sat_v_approx_co2(double tau){
/*
  Approximate saturated vapor density
*/
double tt = 1 - 1/tau;
return exp(
  -1.7074879*pow(tt, 0.340) -
  0.82274670*pow(tt, 0.5) -
  4.6008549*pow(tt, 1.0) -
  10.111178*pow(tt, 7.0/3.0) -
  29.742252*pow(tt, 14.0/3.0)
);
}

double delta_sat_l_approx_co2(double tau){
/*
  Approximate saturated vapor liquid
*/
  double tt = 1 - 1/tau;
  return exp(
    1.9245108*pow(tt, 0.34) -
    0.62385555*pow(tt, 0.5) -
    0.32731127*pow(tt, 10.0/6.0) +
    0.39245142*pow(tt, 11.0/6.0)
  );
}

void phi_co2_ideal_tape(){
// Create a ADOL-C tape for the ideal part of phi for H2O
  double out;

  taped_ideal[comp_enum::h2o] = PHI_IDEAL_TAPE(comp_enum::h2o);
  trace_on(PHI_IDEAL_TAPE(comp_enum::h2o));
  adouble *x, *y;

  double a[] = {
    0,            // 0 not used
    8.37304456,   // 1
    -3.70454304,  // 2
    2.5,          // 3
    1.99427042,   // 4
    0.62105248,   // 5
    0.41195293,   // 6
    1.04028922,   // 7
    0.08327678,   // 8
  };
  double phi[] = {
    0,
    0,
    0,
    0,
    3.15163,   // 4
    6.1119,    // 5
    6.77708,   // 6
    11.32384,  // 7
    27.08792,  // 8
  };

  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  // ideal part
  y[0] = log(x[0]) + a[1] + a[2]*x[1] + a[3]*log(x[1]) +
    a[4]*log(1 - exp(phi[4]*x[1])) +
    a[5]*log(1 - exp(phi[5]*x[1])) +
    a[6]*log(1 - exp(phi[6]*x[1])) +
    a[7]*log(1 - exp(phi[7]*x[1])) +
    a[8]*log(1 - exp(phi[8]*x[1]));

  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}

void phi_co2_resi_tape(){
// Create a ADOL-C tape for the ideal part of phi for H2O

  // parameters are only need once locally to create the tape
  double c[] = {
    0,   // 0 (not used)
    0,   // 1 (not used)
    0,   // 2 (not used)
    0,   // 3 (not used)
    0,   // 4 (not used)
    0,   // 5 (not used)
    0,   // 6 (not used)
    0,   // 7 (not used)
    1,   // 8
    1,   // 9
    1,   // 10
    1,   // 11
    1,   // 12
    1,   // 13
    1,   // 14
    1,   // 15
    1,   // 16
    2,   // 17
    2,   // 18
    2,   // 19
    2,   // 20
    2,   // 21
    2,   // 22
    2,   // 23
    3,   // 24
    3,   // 25
    3,   // 26
    4,   // 27
    4,   // 28
    4,   // 29
    4,   // 30
    4,   // 31
    4,   // 32
    5,   // 33
    6,   // 34
  };

  double d[] = {
    6,   // 0 (not used)
    1,   // 1
    1,   // 2
    1,   // 3
    1,   // 4
    2,   // 5
    2,   // 6
    3,   // 7
    1,   // 8
    2,   // 9
    4,   // 10
    5,   // 11
    5,   // 12
    5,   // 13
    6,   // 14
    6,   // 15
    6,   // 16
    1,   // 17
    1,   // 18
    4,   // 19
    4,   // 20
    4,   // 21
    7,   // 22
    8,   // 23
    2,   // 24
    3,   // 25
    3,   // 26
    5,   // 27
    5,   // 28
    6,   // 29
    7,   // 30
    8,   // 31
    10,  // 32
    4,   // 33
    8,   // 34
    2,   // 35
    2,   // 36
    2,   // 37
    3,   // 38
    3,   // 39
  };

  double t[] = {
    0.00,   // 0 (not used)
    0.00,   // 1
    0.75,   // 2
    1.00,   // 3
    2.00,   // 4
    0.75,   // 5
    2.00,   // 6
    0.75,   // 7
    1.50,   // 8
    1.50,   // 9
    2.50,   // 10
    0.00,   // 11
    1.50,   // 12
    2.00,   // 13
    0.00,   // 14
    1.00,   // 15
    2.00,   // 16
    3.00,   // 17
    6.00,   // 18
    3.00,   // 19
    6.00,   // 20
    8.00,   // 21
    6.00,   // 22
    0.00,   // 23
    7.00,   // 24
    12.00,  // 25
    16.00,  // 26
    22.00,  // 27
    24.00,  // 28
    16.00,  // 29
    24.00,  // 30
    8.00,   // 31
    2.00,   // 32
    28.00,  // 33
    14.00,  // 34
    1.00,   // 35
    0.00,   // 36
    1.00,   // 37
    3.00,   // 38
    3.00,   // 39
  };

  double n[] = {
    0.0,                    // 0 (not used)
    3.88568232031610E-01,   // 1
    2.93854759427400E+00,   // 2
    -5.58671885349340E+00,  // 3
    -7.67531995924770E-01,  // 4
    3.17290055804160E-01,   // 5
    5.48033158977670E-01,   // 6
    1.22794112203350E-01,   // 7
    2.16589615432200E+00,   // 8
    1.58417351097240E+00,   // 9
    -2.31327054055030E-01,  // 10
    5.81169164314360E-02,   // 11
    -5.53691372053820E-01,  // 12
    4.89466159094220E-01,   // 13
    -2.42757398435010E-02,  // 14
    6.24947905016780E-02,   // 15
    -1.21758602252460E-01,  // 16
    -3.70556852700860E-01,  // 17
    -1.67758797004260E-02,  // 18
    -1.19607366379870E-01,  // 29
    -4.56193625087780E-02,  // 20
    3.56127892703460E-02,   // 21
    -7.44277271320520E-03,  // 22
    -1.73957049024320E-03,  // 23
    -2.18101212895270E-02,  // 24
    2.43321665592360E-02,   // 25
    -3.74401334234630E-02,  // 26
    1.43387157568780E-01,   // 27
    -1.34919690832860E-01,  // 28
    -2.31512250534800E-02,  // 29
    1.23631254929010E-02,   // 30
    2.10583219729400E-03,   // 31
    -3.39585190263680E-04,  // 32
    5.59936517715920E-03,   // 33
    -3.03351180556460E-04,  // 34
    -2.13654886883200E+02,  // 35
    2.66415691492720E+04,   // 36
    -2.40272122045570E+04,  // 37
    -2.83416034239990E+02,  // 38
    2.12472844001790E+02,   // 39
    -6.66422765407510E-01,  // 40
    7.26086323498970E-01,   // 41
    5.50686686128420E-02,   // 42
  };

  double a[] = {
    25, // 35
    25, // 36
    25, // 37
    15, // 38
    20, // 39
  };

  double b[] = {
    325, // 35
    300, // 36
    300, // 37
    275, // 38
    275, // 39
  };

  double g[] = {
    1.16, // 35
    1.19, // 36
    1.19, // 37
    1.25, // 38
    1.22, // 39
  };

  double e[] = {
    1, // 35
    1, // 36
    1, // 37
    1, // 38
    1, // 39
  };

  int i = 0;
  double out;

  taped_resi[comp_enum::h2o] = PHI_RESI_TAPE(comp_enum::h2o);
  trace_on(PHI_RESI_TAPE(comp_enum::h2o));
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  y[0] = 0;
  for(i=1; i<=7; ++i){
    y[0] += n[i] * pow(x[0], d[i]) * pow(x[1], t[i]);
  }
  for(i=8; i<=34; ++i){
    y[0] += n[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-pow(x[0], c[i]));
  }
  for(i=35; i<=39; ++i){
    y[0] += n[i] * pow(x[0], d[i]) * pow(x[1], t[i]) *
      exp(
        -a[i-35]*(x[0] - e[i-35])*(x[0] - e[i-35]) -
         b[i-35]*(x[1] - g[i-35])*(x[1] - g[i-35])
      );
  }

  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}

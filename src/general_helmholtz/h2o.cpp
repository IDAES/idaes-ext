#include <adolc/adolc.h>
#include"config.h"

void phi_h2o_ideal_tape(){
// Create a ADOL-C tape for the ideal part of phi for H2O
  double out;

  taped_ideal[comp_enum::h2o] = PHI_IDEAL_TAPE_H2O;
  trace_on(PHI_IDEAL_TAPE_H2O);
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  // ideal part
  y[0] = log(x[0]) -
    8.3204464837497 +
    6.6832105275932*x[1] +
    3.00632*log(x[1]) +
    0.012436*log(1 - exp(-1.28728967*x[1])) +
    0.97315*log(1 - exp(-3.53734222*x[1])) +
    1.27950*log(1 - exp(-7.74073708*x[1])) +
    0.96956*log(1 - exp(-9.24437796*x[1])) +
    0.24873*log(1 - exp(-27.5075105*x[1]));

  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}

void phi_h2o_real_tape(){
// Create a ADOL-C tape for the ideal part of phi for H2O

  // parameters are only need once locally to create the tape
  double c[] = {
    0,   // 0
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
    1,   // 17
    1,   // 18
    1,   // 19
    1,   // 20
    1,   // 21
    1,   // 22
    2,   // 23
    2,   // 24
    2,   // 25
    2,   // 26
    2,   // 27
    2,   // 28
    2,   // 29
    2,   // 30
    2,   // 31
    2,   // 32
    2,   // 33
    2,   // 34
    2,   // 35
    2,   // 36
    2,   // 37
    2,   // 38
    2,   // 39
    2,   // 40
    2,   // 41
    2,   // 42
    3,   // 43
    3,   // 44
    3,   // 45
    3,   // 46
    4,   // 47
    6,   // 48
    6,   // 49
    6,   // 50
    6,   // 51
  };

  double d[] = {
    6,   // 0 (not used)
    1,   // 1
    1,   // 2
    1,   // 3
    2,   // 4
    2,   // 5
    3,   // 6
    4,   // 7
    1,   // 8
    1,   // 9
    1,   // 10
    2,   // 11
    2,   // 12
    3,   // 13
    4,   // 14
    4,   // 15
    5,   // 16
    7,   // 17
    9,   // 18
    10,  // 19
    11,  // 20
    13,  // 21
    15,  // 22
    1,   // 23
    2,   // 24
    2,   // 25
    2,   // 26
    3,   // 27
    4,   // 28
    4,   // 29
    4,   // 30
    5,   // 31
    6,   // 32
    6,   // 33
    7,   // 34
    9,   // 35
    9,   // 36
    9,   // 37
    9,   // 38
    9,   // 39
    10,  // 40
    10,  // 41
    12,  // 42
    3,   // 43
    4,   // 44
    4,   // 45
    5,   // 46
    14,  // 47
    3,   // 48
    6,   // 49
    6,   // 50
    6,   // 51
    3,   // 52
    3,   // 53
    3,   // 54
  };

  double t[] = {
     0.0,    // 0 (not used)
    -0.5,    // 1
     0.875,  // 2
     1,      // 3
     0.5,    // 4
     0.75,   // 5
     0.375,  // 6
     1,      // 7
     4,      // 8
     6,      // 9
    12,      // 10
     1,      // 11
     5,      // 12
     4,      // 13
     2,      // 14
    13,      // 15
     9,      // 16
     3,      // 17
     4,      // 18
    11,      // 19
     4,      // 20
    13,      // 21
     1,      // 22
     7,      // 23
     1,      // 24
     9,      // 25
    10,      // 26
    10,      // 27
     3,      // 28
     7,      // 29
    10,      // 30
    10,      // 31
     6,      // 32
    10,      // 33
    10,      // 34
     1,      // 35
     2,      // 36
     3,      // 37
     4,      // 38
     8,      // 39
     6,      // 40
     9,      // 41
     8,      // 42
    16,      // 43
    22,      // 44
    23,      // 45
    23,      // 46
    10,      // 47
    50,      // 48
    44,      // 49
    46,      // 50
    50,      // 51
     0,      // 52
     1,      // 53
     4,      // 54
  };

  double n[] = {
      0,                    // 0 (not used)
      0.12533547935523e-1,  // 1
      0.78957634722828e1,   // 2
     -0.87803203303561e1,   // 3
      0.31802509345418,     // 4
     -0.26145533859358,     // 5
     -0.78199751687981e-2,  // 6
      0.88089493102134e-2,  // 7
     -0.66856572307965,     // 8
      0.20433810950965,     // 9
     -0.66212605039687e-4,  // 10
     -0.19232721156002,     // 11
     -0.25709043003438,     //
      0.16074868486251,     //
     -0.40092828925807e-1,  //
      0.39343422603254e-6,  //
     -0.75941377088144e-5,  //
      0.56250979351888e-3,  //
     -0.15608652257135e-4,  //
      0.11537996422951e-8,  //
      0.36582165144204e-6,  //
     -0.13251180074668e-11, //
     -0.62639586912454e-9,  //
     -0.10793600908932,     //
      0.17611491008752e-1,  //
      0.22132295167546,     //
     -0.40247669763528,     //
      0.58083399985759,     //
      0.49969146990806e-2,  //
     -0.31358700712549e-1,  //
     -0.74315929710341,     //
      0.47807329915480,     //
      0.20527940895948e-1,  //
     -0.13636435110343,     //
      0.14180634400617e-1,  //
      0.83326504880713e-2,  //
     -0.29052336009585e-1,  //
      0.38615085574206e-1,  //
     -0.20393486513704e-1,  //
     -0.16554050063743e-2,  //
      0.19955571979541e-2,  //
      0.15870308324157e-3,  //
     -0.16388568342530e-4,  //
      0.43613615723811e-1,  //
      0.34994005463765e-1,  //
     -0.76788197844621e-1,  //
      0.22446277332006e-1,  //
     -0.62689710414685e-4,  //
     -0.55711118565645e-9,  // 48
     -0.19905718354408,     // 49
      0.31777497330738,     // 50
     -0.11841182425981,     // 51
     -0.31306260323435e2,   // 52
      0.31546140237781e2,   // 53
     -0.25213154341695e4,   // 54
     -0.14874640856724,     // 55
      0.31806110878444,     // 56
  };

  double a[] = {
    20, // 52
    20, // 53
    20, // 54
  };

  double b[] = {
    150, // 52
    150, // 53
    250, // 54
  };

  double g[] = {
    1.21, // 52
    1.21, // 53
    1.25, // 54
  };

  double e[] = {
    1, // 52
    1, // 53
    1, // 54
  };
  int i = 0;
  double out;

  taped_real[comp_enum::h2o] = PHI_REAL_TAPE_H2O;
  trace_on(PHI_REAL_TAPE_H2O);
  adouble *x, *y;
  x = new adouble[2];
  y = new adouble[1];
  x[0] <<= 0.9;
  x[1] <<= 1.1;

  y[0] = 0;
  for(i=1; i<=7; ++i){
    y[0] += n[i] * pow(x[0], d[i]) * pow(x[1], t[i]);
  }
  for(i=8; i<=51; ++i){
    y[0] += n[i] * pow(x[0], d[i]) * pow(x[1], t[i]) * exp(-pow(x[0], c[i]));
  }
  for(i=52; i<=54; ++i){
    y[0] += n[i] * pow(x[0], d[i]) * pow(x[1], t[i]) *
      exp(
        -a[i-52]*(x[0] - e[i-52])*(x[0] - e[i-52]) -
         b[i-52]*(x[1] - g[i-52])*(x[1] - g[i-52])
      );
  }

  y[0] >>= out;
  delete[] y;
  delete[] x;
  trace_off();
}
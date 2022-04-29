/*
  General root calculator for cubic equations of state.
  Author: John Eslick
  Date: November 22, 2016
  Notes:
*/

#include "cubic_roots.h"

/***********************************************************************
 * CUBIC FORMULA FUNCTION
 * The root finding approach given in CRC Standard Mathematical
 * Tables and Formulae 32nd edition pg 68 was used to identify the
 * roots. There maybe a mistake in it depending on your printing,
 * if you want to check on it also check out the errata:
 * http://www.mathtable.com/smtf/
 **********************************************************************/

double cubic_root_low(double b, double c, double d){
 double F, G, H, I, J, K, M, N, P, R, S, T, U, zr[3];

 F = (3*c - b*b)/3.0;
 G = (2*b*b*b - 9*b*c + 27*d)/27.0;
 H = G*G/4.0 + F*F*F/27.0;
 P = -b/3.0;
 if(H <= 0.0){
     //Three roots, can include double or triple roots too
     I = sqrt(G*G/4.0 - H);
     J = cbrt(I);
     K = acos(-G/2.0/I);
     M = cos(K/3);
     N = sqrt_3*sin(K/3.0);
     zr[0] = P + 2*J*M;
     zr[1] = P - J*(M + N);
     zr[2] = P - J*(M - N);
     //Sort to get lowest (don't need full sort)
     if(zr[1] > zr[2]) zr[1] = zr[2];
     if(zr[0] > zr[1]) return zr[1];
     return zr[0];
 }
 R = -G/2.0 + sqrt(H);
 T = -G/2.0 - sqrt(H);
 S = cbrt(R);
 U = cbrt(T);
 return S + U + P;
}

double cubic_root_high(double b, double c, double d){
 double F, G, H, I, J, K, M, N, P, R, S, T, U, zr[3];

 F = (3*c - b*b)/3.0;
 G = (2*b*b*b - 9*b*c + 27*d)/27.0;
 H = G*G/4.0 + F*F*F/27.0;
 P = -b/3.0;
 if(H <= 0.0){
     //Three roots, can include double or triple roots too
     I = sqrt(G*G/4.0 - H);
     J = cbrt(I);
     K = acos(-G/2.0/I);
     M = cos(K/3);
     N = sqrt_3*sin(K/3.0);
     zr[0] = P + 2*J*M;
     zr[1] = P - J*(M + N);
     zr[2] = P - J*(M - N);
     //Sort to get highest (don't need full sort)
     if(zr[0] > zr[1]) zr[1] = zr[0];
     if(zr[1] > zr[2]) return zr[1];
     return zr[2];
 }
 R = -G/2.0 + sqrt(H);
 T = -G/2.0 - sqrt(H);
 S = cbrt(R);
 U = cbrt(T);
 return S + U + P;
}

double cubic_root_low_with_deriv(double b, double c, double d, double *derivs, double *hes){
  double z;
  z = cubic_root_low(b, c, d); // Liquid low compressibility
  if (derivs != NULL) {
    derivs[0] = 1.0/(-1.0 + c/z/z + 2*d/z/z/z); // dz/db
    derivs[1] = 1.0/(-2.0*z - b + d/z/z); // dz/dc
    derivs[2] = 1.0/(-3.0*z*z - 2*b*z - c); // dz/dd
  }
  if (hes != NULL){
    hes[0] = derivs[0]*derivs[0]*derivs[0]*(2*c/z/z/z + 6*d/z/z/z/z); // dz2/db2
    hes[1] = derivs[0]*derivs[0]*(2*c/z/z/z*derivs[1] + 6*d/z/z/z/z*derivs[1] - 1/z/z); // dz2/dbdc
    hes[2] = derivs[1]*derivs[1]*derivs[1]*(2 + 2*d/z/z/z); // dz2/dc2
    hes[3] = derivs[0]*derivs[0]*(2*c/z/z/z*derivs[2] + 6*d/z/z/z/z*derivs[2] - 2/z/z/z); // dz2/dbdd
    hes[4] = derivs[1]*derivs[1]*(2*derivs[2] + 2*d/z/z/z*derivs[2] - 1/z/z); // dz2/dcdd
    hes[5] = derivs[2]*derivs[2]*derivs[2]*(6*z + 2*b); // dz2/dd2
  }
  return z;
}

double cubic_root_high_with_deriv(double b, double c, double d, double *derivs, double *hes){
  double z;
  z = cubic_root_high(b, c, d); // Liquid low compressibility
  if (derivs != NULL) {
    derivs[0] = 1.0/(-1.0 + c/z/z + 2*d/z/z/z); // dz/db
    derivs[1] = 1.0/(-2.0*z - b + d/z/z); // dz/dc
    derivs[2] = 1.0/(-3.0*z*z - 2*b*z - c); // dz/dd
  }
  if (hes != NULL){
    hes[0] = derivs[0]*derivs[0]*derivs[0]*(2*c/z/z/z + 6*d/z/z/z/z); // dz2/db2
    hes[1] = derivs[0]*derivs[0]*(2*c/z/z/z*derivs[1] + 6*d/z/z/z/z*derivs[1] - 1/z/z); // dz2/dbdc
    hes[2] = derivs[1]*derivs[1]*derivs[1]*(2 + 2*d/z/z/z); // dz2/dc2
    hes[3] = derivs[0]*derivs[0]*(2*c/z/z/z*derivs[2] + 6*d/z/z/z/z*derivs[2] - 2/z/z/z); // dz2/dbdd
    hes[4] = derivs[1]*derivs[1]*(2*derivs[2] + 2*d/z/z/z*derivs[2] - 1/z/z); // dz2/dcdd
    hes[5] = derivs[2]*derivs[2]*derivs[2]*(6*z + 2*b); // dz2/dd2
  }
  return z;
}

double cubic_root_low_ext_with_deriv(double b, double c, double d, double *derivs, double *hes){
  double x, z, t1, t2, ft1, ft2, det, A;
  double dt1db, dt2db, dt1dc, dt2dc, dfdb_at_t1, dfdb_at_t2, dfdc_at_t1, dfdc_at_t2;
  double dxdb, dxdc, dxdd, dAdb, dAdc, dAdd;
  double d2t1db2, d2t2db2, d2t1dbdc, d2t2dbdc, d2t1dc2, d2t2dc2;
  double d2fdb2_at_t1, d2fdb2_at_t2, d2fdbdc_at_t1;
  double d2fdbdc_at_t2, d2fdc2_at_t1, d2fdc2_at_t2;

  det = 4*b*b - 12*c;
  t1 = (-2*b - sqrt(det))/6.0;
  t2 = (-2*b + sqrt(det))/6.0;

  ft1 = t1*t1*t1 + b*t1*t1 + c*t1 + d;
  ft2 = t2*t2*t2 + b*t2*t2 + c*t2 + d;
  z = cubic_root_high(b, c, d - ft2 + ft1) + t1 - t2;
  x = z - t1 + t2;
  A = 3*x*x + 2*b*x + c;
  dt1db = (-1 - 2*b/sqrt(4*b*b - 12*c))/3.0;
  dt2db = (-1 + 2*b/sqrt(4*b*b - 12*c))/3.0;
  dt1dc = 1.0/sqrt(4*b*b - 12*c);
  dt2dc = -1.0/sqrt(4*b*b - 12*c);
  dfdb_at_t1 = t1*t1 + 3*t1*t1*dt1db + 2*b*t1*dt1db + c*dt1db;
  dfdb_at_t2 = t2*t2 + 3*t2*t2*dt2db + 2*b*t2*dt2db + c*dt2db;
  dfdc_at_t1 = t1 + 3*t1*t1*dt1dc + 2*b*t1*dt1dc + c*dt1dc;
  dfdc_at_t2 = t2 + 3*t2*t2*dt2dc + 2*b*t2*dt2dc + c*dt2dc;

  derivs[0] = (-x*x + A*(dt1db - dt2db) - dfdb_at_t1 + dfdb_at_t2)/A;
  derivs[1] = (-x + A*(dt1dc - dt2dc) - dfdc_at_t1 + dfdc_at_t2)/A;
  derivs[2] = -1/A;

  dxdb = derivs[0] - dt1db + dt2db;
  dxdc = derivs[1] - dt1dc + dt2dc;
  dxdd = derivs[2];
  dAdb = 2*x + 6*x*dxdb + 2*b*dxdb;
  dAdc = 1 + 6*x*dxdc + 2*b*dxdc;
  dAdd = 6*x*dxdd + 2*b*dxdd;

  d2t1db2 = -2.0/3.0/sqrt(4*b*b - 12*c) + 8/3.0*b*b*pow(4*b*b - 12*c, -1.5);
  d2t2db2 = 2.0/3.0/sqrt(4*b*b - 12*c) - 8/3.0*b*b*pow(4*b*b - 12*c, -1.5);
  d2t1dbdc = -4*b*pow(4*b*b - 12*c, -1.5);
  d2t2dbdc = 4*b*pow(4*b*b - 12*c, -1.5);
  d2t1dc2 = 6*pow(4*b*b - 12*c, -1.5);
  d2t2dc2 = -6*pow(4*b*b - 12*c, -1.5);

  d2fdb2_at_t1 = 2*t1*dt1db + 6*t1*dt1db*dt1db + 3*t1*t1*d2t1db2 + 2*t1*dt1db + 2*b*dt1db*dt1db + 2*b*t1*d2t1db2 + c*d2t1db2;
  d2fdb2_at_t2 = 2*t2*dt2db + 6*t2*dt2db*dt2db + 3*t2*t2*d2t2db2 + 2*t2*dt2db + 2*b*dt2db*dt2db + 2*b*t2*d2t2db2 + c*d2t2db2;
  d2fdbdc_at_t1 = 2*t1*dt1dc + 6*t1*dt1db*dt1dc + 3*t1*t1*d2t1dbdc + 2*b*dt1db*dt1dc + 2*b*t1*d2t1dbdc + c*d2t1dbdc + dt1db;
  d2fdbdc_at_t2 = 2*t2*dt2dc + 6*t2*dt2db*dt2dc + 3*t2*t2*d2t2dbdc + 2*b*dt2db*dt2dc + 2*b*t2*d2t2dbdc + c*d2t2dbdc + dt2db;
  d2fdc2_at_t1 = dt1dc + 6*t1*dt1dc*dt1dc + 3*t1*t1*d2t1dc2 + 2*b*dt1dc*dt1dc + 2*b*t1*d2t1dc2 + c*d2t1dc2 + dt1dc;
  d2fdc2_at_t2 = dt2dc + 6*t2*dt2dc*dt2dc + 3*t2*t2*d2t2dc2 + 2*b*dt2dc*dt2dc + 2*b*t2*d2t2dc2 + c*d2t2dc2 + dt2dc;

  hes[0] = -1/A*derivs[0]*dAdb + 1/A*(-2*x*dxdb + dAdb*dt1db + A*d2t1db2 - dAdb*dt2db - A*d2t2db2 - d2fdb2_at_t1 + d2fdb2_at_t2);
  hes[1] = -1/A*derivs[0]*dAdc + 1/A*(-2*x*dxdc + dAdc*dt1db + A*d2t1dbdc - dAdc*dt2db - A*d2t2dbdc - d2fdbdc_at_t1 + d2fdbdc_at_t2);
  hes[2] = 1/A/A*dAdb;
  hes[3] = -1/A*derivs[1]*dAdc + 1/A*(-dxdc + dAdc*dt1dc + A*d2t1dc2 - dAdc*dt2dc - A*d2t2dc2 - d2fdc2_at_t1 + d2fdc2_at_t2);
  hes[4] = 1/A/A*dAdc;
  hes[5] = 1/A/A*dAdd;

  return z;
}

double cubic_root_high_ext_with_deriv(double b, double c, double d, double *derivs, double *hes){
  double x, z, t1, t2, ft1, ft2, det, A;
  double dt1db, dt2db, dt1dc, dt2dc, dfdb_at_t1, dfdb_at_t2, dfdc_at_t1, dfdc_at_t2;
  double dxdb, dxdc, dxdd, dAdb, dAdc, dAdd;
  double d2t1db2, d2t2db2, d2t1dbdc, d2t2dbdc, d2t1dc2, d2t2dc2;
  double d2fdb2_at_t1, d2fdb2_at_t2, d2fdbdc_at_t1;
  double d2fdbdc_at_t2, d2fdc2_at_t1, d2fdc2_at_t2;

  det = 4*b*b - 12*c;
  t1 = (-2*b - sqrt(det))/6.0;
  t2 = (-2*b + sqrt(det))/6.0;

  ft1 = t1*t1*t1 + b*t1*t1 + c*t1 + d;
  ft2 = t2*t2*t2+ b*t2*t2 + c*t2 + d;
  z = cubic_root_low(b, c, d + ft2 - ft1) - t1 + t2;
  x = z + t1 - t2;
  A = 3*x*x + 2*b*x + c;
  dt1db = (-1 - 2*b/sqrt(4*b*b - 12*c))/3.0;
  dt2db = (-1 + 2*b/sqrt(4*b*b - 12*c))/3.0;
  dt1dc = 1.0/sqrt(4*b*b - 12*c);
  dt2dc = -1.0/sqrt(4*b*b - 12*c);
  dfdb_at_t1 = t1*t1 + 3*t1*t1*dt1db + 2*b*t1*dt1db + c*dt1db;
  dfdb_at_t2 = t2*t2 + 3*t2*t2*dt2db + 2*b*t2*dt2db + c*dt2db;
  dfdc_at_t1 = t1 + 3*t1*t1*dt1dc + 2*b*t1*dt1dc + c*dt1dc;
  dfdc_at_t2 = t2 + 3*t2*t2*dt2dc + 2*b*t2*dt2dc + c*dt2dc;

  derivs[0] = (-x*x + A*(-dt1db + dt2db) + dfdb_at_t1 - dfdb_at_t2)/A;
  derivs[1] = (-x + A*(-dt1dc + dt2dc) + dfdc_at_t1 - dfdc_at_t2)/A;
  derivs[2] = -1/A;

  dxdb = derivs[0] + dt1db - dt2db;
  dxdc = derivs[1] + dt1dc - dt2dc;
  dxdd = derivs[2];
  dAdb = 2*x + 6*x*dxdb + 2*b*dxdb;
  dAdc = 1 + 6*x*dxdc + 2*b*dxdc;
  dAdd = 6*x*dxdd + 2*b*dxdd;

  d2t1db2 = -2.0/3.0/sqrt(4*b*b - 12*c) + 8/3.0*b*b*pow(4*b*b - 12*c, -1.5);
  d2t2db2 = 2.0/3.0/sqrt(4*b*b - 12*c) - 8/3.0*b*b*pow(4*b*b - 12*c, -1.5);
  d2t1dbdc = -4*b*pow(4*b*b - 12*c, -1.5);
  d2t2dbdc = 4*b*pow(4*b*b - 12*c, -1.5);
  d2t1dc2 = 6*pow(4*b*b - 12*c, -1.5);
  d2t2dc2 = -6*pow(4*b*b - 12*c, -1.5);

  d2fdb2_at_t1 = 2*t1*dt1db + 6*t1*dt1db*dt1db + 3*t1*t1*d2t1db2 + 2*t1*dt1db + 2*b*dt1db*dt1db + 2*b*t1*d2t1db2 + c*d2t1db2;
  d2fdb2_at_t2 = 2*t2*dt2db + 6*t2*dt2db*dt2db + 3*t2*t2*d2t2db2 + 2*t2*dt2db + 2*b*dt2db*dt2db + 2*b*t2*d2t2db2 + c*d2t2db2;
  d2fdbdc_at_t1 = 2*t1*dt1dc + 6*t1*dt1db*dt1dc + 3*t1*t1*d2t1dbdc + 2*b*dt1db*dt1dc + 2*b*t1*d2t1dbdc + c*d2t1dbdc + dt1db;
  d2fdbdc_at_t2 = 2*t2*dt2dc + 6*t2*dt2db*dt2dc + 3*t2*t2*d2t2dbdc + 2*b*dt2db*dt2dc + 2*b*t2*d2t2dbdc + c*d2t2dbdc + dt2db;
  d2fdc2_at_t1 = dt1dc + 6*t1*dt1dc*dt1dc + 3*t1*t1*d2t1dc2 + 2*b*dt1dc*dt1dc + 2*b*t1*d2t1dc2 + c*d2t1dc2 + dt1dc;
  d2fdc2_at_t2 = dt2dc + 6*t2*dt2dc*dt2dc + 3*t2*t2*d2t2dc2 + 2*b*dt2dc*dt2dc + 2*b*t2*d2t2dc2 + c*d2t2dc2 + dt2dc;

  hes[0] = -1/A*derivs[0]*dAdb + 1/A*(-2*x*dxdb - dAdb*dt1db - A*d2t1db2 + dAdb*dt2db + A*d2t2db2 + d2fdb2_at_t1 - d2fdb2_at_t2);
  hes[1] = -1/A*derivs[0]*dAdc + 1/A*(-2*x*dxdc - dAdc*dt1db - A*d2t1dbdc + dAdc*dt2db + A*d2t2dbdc + d2fdbdc_at_t1 - d2fdbdc_at_t2);
  hes[2] = 1/A/A*dAdb;
  hes[3] = -1/A*derivs[1]*dAdc + 1/A*(-dxdc - dAdc*dt1dc - A*d2t1dc2 + dAdc*dt2dc + A*d2t2dc2 + d2fdc2_at_t1 - d2fdc2_at_t2);
  hes[4] = 1/A/A*dAdc;
  hes[5] = 1/A/A*dAdd;

  return z;
}


/***********************************************************************
 * Functions to return the cubic roots
 **********************************************************************/

real cubic_root_l(arglist *al){
  return cubic_root_low_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
}

real cubic_root_h(arglist *al){
  return cubic_root_high_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
}

real cubic_root_l_nan(arglist *al){
  real z;
  z = cubic_root_low_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
  if (0 >= al->ra[al->at[0]]*al->ra[al->at[0]] - 3*al->ra[al->at[1]]) return z; // no turning points so use same
  if (z <= -al->ra[al->at[0]]/3.0) return z;
  return NAN;
}

real cubic_root_h_nan(arglist *al){
  real z;
  z = cubic_root_low_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
  if (0 >= al->ra[al->at[0]]*al->ra[al->at[0]] - 3*al->ra[al->at[1]]) return z; // no turning points so use same
  if (z >= -al->ra[al->at[0]]/3.0) return z;
  return NAN;
}

real cubic_root_l_ext(arglist *al){
  real z;
  z = cubic_root_low_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
  if (z <= -al->ra[al->at[0]]/3.0) return z;
  if (0 >= al->ra[al->at[0]]*al->ra[al->at[0]] - 3*al->ra[al->at[1]]) return z; // no turning points so use same
  return cubic_root_low_ext_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
}

real cubic_root_h_ext(arglist *al){
  real z;
  z = cubic_root_high_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
  if (z >= -al->ra[al->at[0]]/3.0) return z;
  if (0 >= al->ra[al->at[0]]*al->ra[al->at[0]] - 3*al->ra[al->at[1]]) return z; // no turning points so use same
  return cubic_root_high_ext_with_deriv(
      al->ra[al->at[0]],
      al->ra[al->at[1]],
      al->ra[al->at[2]],
      al->derivs, al->hes
  );
}

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info
     */
    //Real value function, and suppress anoying warings that, I think
    //happen when this gets called on an already loaded library.
    int t = FUNCADD_REAL_VALUED;
    addfunc("cubic_root_l", (rfunc)cubic_root_l, t, -1, NULL);
    addfunc("cubic_root_h", (rfunc)cubic_root_h, t, -1, NULL);
    addfunc("cubic_root_l_nan", (rfunc)cubic_root_l_nan, t, -1, NULL);
    addfunc("cubic_root_h_nan", (rfunc)cubic_root_h_nan, t, -1, NULL);
    addfunc("cubic_root_l_ext", (rfunc)cubic_root_l_ext, t, -1, NULL);
    addfunc("cubic_root_h_ext", (rfunc)cubic_root_h_ext, t, -1, NULL);
}

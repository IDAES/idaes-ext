/*
  General root calculator for cubic equations of state.

  Author: John Eslick
  Date: November 22, 2016

  Notes:
*/
#include <time.h>
#include <math.h>
#include "funcadd.h"
#include <stdio.h>

/***********************************************************************
 *
 * CONSTANTS
 *
 **********************************************************************/
static const double deriv_cap=1e10; // Largest derivative magnitude allowed
static const double sqrt_3=1.73205080756888; //square root of 3

double cubic_root_low(double b, double c, double d);
double cubic_root_low(double b, double c, double d);
double cubic_root_low_with_deriv(double b, double c, double d, double *derivs, double *hes);
double cubic_root_hign_with_deriv(double b, double c, double d, double *derivs, double *hes);

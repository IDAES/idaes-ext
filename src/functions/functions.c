#include "functions.h"
#include <math.h>

#undef printf

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info */
    addfunc("cbrt", (rfunc)scbrt, FUNCADD_REAL_VALUED, 1, NULL);
    addfunc("x_over_exp_x_minus_one", (rfunc)x_over_exp_x_minus_one, FUNCADD_REAL_VALUED, 1, NULL);
}

extern real scbrt(arglist *al){
    // al is the argument list data structure
    // al->ra is an array of real arguments
    // al->at is an array of of argument positions
    // The reason for using al->at for the position is
    // that there could also be integer or string
    // arguments
    real x = al->ra[al->at[0]];

    //al->derivs is a pointer to an array of first derivatives
    //of the function with respect to each real arg if
    //derivs is NULL solver isn't requesting derivatives.
    if(al->derivs!=NULL){
      if(fabs(x) < 6e-9) al->derivs[0] = 1e5;
      else al->derivs[0] = pow(cbrt(x), -2.0)/3.0;
      //al->hes is a pointer the Hessian matrix if NULL second
      //derivatives are not requested.  This function takes a
      //single argument, but the Hessian form in general is
      //the upper triangle in column-major form like
      // Args Index  0 1 2 3 ...
      //          0  0 1 3 6
      //          1    2 4 7
      //          2      5 8
      //          3        9
      //        ...
      if(al->hes!=NULL){
        al->hes[0] = -2.0*pow(cbrt(x), -5.0)/9.0;
      }
    }
    return cbrt(x);
}

extern real x_over_exp_x_minus_one(arglist *al){
    real x = al->ra[al->at[0]];

    if(al->derivs!=NULL){
      if(fabs(x) < 0.01){
        al->derivs[0] = -0.5 + x/6.0 + x*x*x/180.0 - x*x*x*x*x/5040.0;
      } 
      else {
        al->derivs[0] = 1.0/(exp(x) - 1) - x*exp(x)/(exp(x) - 1)/(exp(x) - 1);
      }

      if(al->hes!=NULL){
        if(fabs(x) < 0.01){
          al->hes[0] = 1.0/6.0 + x*x/60 - x*x*x*x/1008.0;
        }
        else{
          al->hes[0] = exp(x)*(exp(x)*(x - 2) + x + 2) /
            (exp(x) - 1)/(exp(x) - 1)/(exp(x) - 1);
        }
      }
    }
    if(fabs(x) < 0.01){
      return  1 - x/2.0 + x*x/12.0 + x*x*x*x/720.0 - x*x*x*x*x*x/30240.0;
    }
    return x/(exp(x) - 1);
}

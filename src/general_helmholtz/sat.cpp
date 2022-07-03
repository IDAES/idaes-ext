/*------------------------------------------------------------------------------
Author: John Eslick
File sat.cpp

 This file contains the functions to solve the saturation curve.  The method
 is described in the following paper:

 Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
     State from Helmholtz Energy Equations of State." Journal of Thermal
     Science and Technology, 3(3), 442-451.
------------------------------------------------------------------------------*/
#include<vector>
#include<math.h>
#include<iostream>
#include"phi.h"
#include"config.h"

inline double J(double delta, std::vector<double> *phi){
  // Term from Akasaka method for saturation state
  return delta*(1 + delta*phi->at(1));
}

inline double K(double delta, std::vector<double> *phi){
  // Term from Akasaka method for saturation state
  return delta*phi->at(1) + phi->at(0) + log(delta);
}

inline double J_delta(double delta, std::vector<double> *phi){
  // Derivative term from Akasaka method for saturation state
  return 1.0 + 2.0*delta*phi->at(1) + delta*delta*phi->at(2);
}

inline double K_delta(double delta, std::vector<double> *phi){
  // Derivative term from Akasaka method for saturation state
  return 2.0*phi->at(1) + delta*phi->at(2) + 1.0/delta;
}

int sat(comp_enum comp, double tau, double *delta_l, double *delta_v){
  double Kdiff=1, Jdiff=1, Jv, Jl, Kv, Kl, det, dJv, dJl, dKv, dKl;
  int n=0;
  std::vector<double> phir_v, phir_l;

  static const double tol_sat = 1e-11;
  static const double sat_gamma = 1.0;
  static const unsigned int max_iter = 100;

  if(tau - 1 < tol_sat){
    // So close to the critical point assume we are at the pritical point, if
    // over the critical point this will also use the critical density
    *delta_l = 1.0;
    *delta_v = 1.0;
  }
  else{
    // okay so you've decided to solve this thing
    *delta_l = 2.89; // DELTA_LIQ_SAT_GUESS;
    *delta_v = 0.0149; // DELTA_VAP_SAT_GUESS;
    while(n < max_iter){
      phi_real_for_sat(comp, *delta_v, tau, &phir_v);
      phi_real_for_sat(comp, *delta_l, tau, &phir_l);
      Jv = J(*delta_v, &phir_v);
      Jl = J(*delta_l, &phir_l);
      Kv = K(*delta_v, &phir_v);
      Kl = K(*delta_l, &phir_l);
      std::cout << Jv << " " << Jl << " " << Kv << " " << Kl << std::endl;
      Jdiff = Jv - Jl;
      Kdiff = Kv - Kl;
      if (fabs(Jdiff) < tol_sat && fabs(Kdiff) < tol_sat){
        break;
      }
      dJv = J_delta(*delta_v, &phir_v);
      dJl = J_delta(*delta_l, &phir_l);
      dKv = K_delta(*delta_v, &phir_v);
      dKl = K_delta(*delta_l, &phir_l);
      det = dJv*dKl - dJl*dKv;
      ++n; // Count iterations
      *delta_l += sat_gamma*(Kdiff*dJv - Jdiff*dKv)/det;
      *delta_v += sat_gamma*(Kdiff*dJl - Jdiff*dKl)/det;
    }
  }
  /*
  //Calculate grad and hes for and memoize
  double delta_l = *delta_l_sol, delta_v = *delta_v_sol;
  gradv[0] = LHM/LGM;
  gradl[0] = gradv[0]*LBV/LBL + (LCV - LCL)/LBL;
  hesv[0] = LdHdt(delta_l, delta_v, tau, gradl[0], gradv[0])/LGM
          - LHM/LGM/LGM*LdGdt(delta_l, delta_v, tau, gradl[0], gradv[0]);
  hesl[0] = hesv[0]*LBV*LFL + gradv[0]*(LBVt + LBVd*gradv[0])*LFL
           + gradv[0]*LBV*(LFLt + LFLd*gradl[0]) + (LFLt + LFLd*gradl[0])*(LCV - LCL)
           + LFL*(LCVt - LCLt + LCVd*gradv[0] - LCLd*gradl[0]);
  */
  return n;
}

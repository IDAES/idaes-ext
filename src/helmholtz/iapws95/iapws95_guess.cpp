#include"helmholtz_solve.h"
#include"iapws95_param.h"

#include<stdio.h>
#include<cmath>
#include<iostream>

s_real p_sat_iapws97(s_real tau){ //saturation pressure from tau IAPWS-97 eq.
  //the IAPWS-97 isn't as consistent as IAPWS-95, but this provides a good guess
  static const s_real n_psat[] = {
     0.11670521452767e4, //1
    -0.72421316703206e6, //2
    -0.17073846940092e2, //3
     0.12020824702470e5, //4
    -0.32325550322333e7, //5
     0.14915108613530e2, //6
    -0.48232657361591e4, //7
     0.40511340542057e6, //8
    -0.23855557567849,   //9
     0.65017534844798e3  //10
  };
  s_real T = T_c/tau;
  s_real tt = T + n_psat[8]/(T - n_psat[9]);
  s_real A = tt*tt + n_psat[0]*tt + n_psat[1];
  s_real B = n_psat[2]*tt*tt + n_psat[3]*tt + n_psat[4];
  s_real C = n_psat[5]*tt*tt + n_psat[6]*tt + n_psat[7];
  return 1000*pow(2*C/(-B + pow(B*B - 4*A*C, 0.5)), 4);
}

s_real delta_sat_v_approx_iapws95(s_real tau){ //approximate saturated vapor density
  // This equation is from the original IAPWS-95 paper
  s_real XX = 1 - 1.0/tau;
  s_real delta = exp(-2.03150240*pow(XX,2.0/6.0)
            - 2.68302940*pow(XX,4.0/6.0)
            - 5.38626492*pow(XX,8.0/6.0)
            - 17.2991605*pow(XX,18.0/6.0)
            - 44.7586581*pow(XX,37.0/6.0)
            - 63.9201063*pow(XX,71.0/6.0));
  return delta;
}

s_real delta_sat_l_approx_iapws95(s_real tau){ //approximate saturated vapor density
  // This equation is from the original IAPWS-95 paper.
  s_real XX = 1 - 1.0/tau;
  s_real delta = 1.001
           + 1.99274064*pow(XX,1.0/3.0)
           + 1.09965342*pow(XX,2.0/3.0)
           - 0.510839303*pow(XX,5.0/3.0)
           - 1.75493479*pow(XX,16.0/3.0)
           - 45.5170352*pow(XX,43.0/3.0)
           - 6.74694450e5*pow(XX,110.0/3.0);
  return delta;
}

s_real delta_p_tau_supercritical_guess_iapws95(s_real p, s_real tau){
  //since it's super critical there is not two phase region so there shouldn't
  //be multiple solutions to worry about.  There are some steep slopes in places
  //so I'm using a flase position backeting method to get close.
  return delta_p_tau_rf(p, tau, 1e-5/rho_c, 1100/rho_c);
}

s_real delta_p_tau_vap_guess_iapws95(s_real p, s_real tau){
  if (p >= P_c && tau <= 1){
    return delta_p_tau_supercritical_guess_iapws95(p, tau);
  }
  return delta_p_tau_rf(p, tau, 5e-6/rho_c, sat_delta_vap(tau));
}

s_real delta_p_tau_liq_guess_iapws95(s_real p, s_real tau){
  s_real T=T_c/tau;

  if (p >= P_c && tau <= 1){
    return delta_p_tau_supercritical_guess_iapws95(p, tau);
  }
  else if(T > 645.0 && p > 21500.0){
    return delta_p_tau_rf(p, tau, sat_delta_liq(tau), 450.0/rho_c);
  }
  return 1000.0/rho_c;
}

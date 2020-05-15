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

  /* This uses approximate isochors to get a reasonable initial guess for the
  density solution for a vapor at a particular temperature and pressure.  If
  there are tricky regions where a really good guess is needed, a bracketing
  method can be used to refine the guess.*/
  s_real T=T_c/tau;

  if (p - P_c < 200 && T - T_c < 5){ //close to critical stay near critical density
      return delta_p_tau_rf(p, tau, 290/rho_c, 390/rho_c);
  }
  else if(p < 7.29146E-09 + 4.61516E-06*T + 7.28228E-15*T*T){ //rho < 1e-5
    return 5e-6/rho_c;
  }
  else if(p < -1.50992E-07 + 4.61521E-05*T - 1.89168E-13*T*T){ //rho < 1e-4
    return 5e-5/rho_c;
  }
  else if(p < -1.38128E-05 + 0.000461547*T - 1.47273E-11*T*T){ //rho < 1e-3
    return 5e-4/rho_c;
  }
  else if(p < -0.001212888 + 0.004617609*T - 1.22486E-09*T*T){ //rho < 1e-2
    return 5e-3/rho_c;
  }
  else if(p < -0.090803154 + 0.046320106*T - 8.03954E-08*T*T){ //rho < 1e-1
    return 5e-2/rho_c;
  }
  else if(p < -6.934621112 + 0.473323909*T - 5.27629E-06*T*T){ //rho < 1
    return 5e-1/rho_c;
  }
  else if(p < -532.9921491 + 5.433894672*T - 0.0003339*T*T){ //rho < 10
    return delta_p_tau_rf(p, tau, 0.5/rho_c, 15/rho_c);
  }
  else if(p < -34931.47795 + 92.24848843*T - 0.015853248*T*T){ //rho < 100
    return delta_p_tau_rf(p, tau, 5/rho_c, 105/rho_c);
  }
  else if(p < -99354.97899 + 209.1578681*T - 0.032762406*T*T){ //rho < 200
    return delta_p_tau_rf(p, tau, 95/rho_c, 205/rho_c);
  }
  else if(p < -167178.3915 + 313.7993977*T - 0.033828951*T*T){ //rho < 300
    return delta_p_tau_rf(p, tau, 195/rho_c, 305/rho_c);
  }
  else if(p < -245383.5357 + 428.1316342*T - 0.027047284*T*T){ //rho < 400
    return delta_p_tau_rf(p, tau, 295/rho_c, 405/rho_c);
  }
  else if(p < -350195.9623 + 593.7310153*T - 0.02847895*T*T){ //rho < 500
    return delta_p_tau_rf(p, tau, 395/rho_c, 505/rho_c);
  }
  else if(p < -490776.283 + 849.812898*T - 0.056037941*T*T){ //rho < 600
    return delta_p_tau_rf(p, tau, 495/rho_c, 605/rho_c);
  }
  else if(p < -652011.1617 + 1201.907668*T - 0.113244604*T*T){ //rho < 700
    return delta_p_tau_rf(p, tau, 595/rho_c, 705/rho_c);
  }
  else if(p < -801476.7975 + 1633.589044*T - 0.19012386*T*T){ //rho < 800
    return delta_p_tau_rf(p, tau, 695/rho_c, 805/rho_c);
  }
  else if(p < -850674.915 + 2003.358868*T - 0.196401469*T*T){ //rho < 900
    return delta_p_tau_rf(p, tau, 795/rho_c, 905/rho_c);
  }
  else if(p < -417462.7784 + 1107.467347*T + 0.884115157*T*T){ //rho < 1000
    return delta_p_tau_rf(p, tau, 895/rho_c, 1005/rho_c);
  }
  return 1100/rho_c;
}

s_real delta_p_tau_vap_guess_iapws95(s_real p, s_real tau){

  /* This uses approximate isochors to get a reasonable initial guess for the
  density solution for a vapor at a particular temperature and pressure.  If
  there are tricky regions where a really good guess is needed, a bracketing
  method can be used to refine the guess.*/
  s_real T=T_c/tau;

  if (p >= P_c && tau <= 1){
    return delta_p_tau_supercritical_guess_iapws95(p, tau);
  }
  // for the rest just use a guess between isochors, while it's tempting to use
  // regula falsi between the isochors, don't do it.  There are multiple
  // solutions in there in some places.  It seems a low guess arrives at the
  // lower correct density for vapor
  else if(p < 7.29146E-09 + 4.61516E-06*T + 7.28228E-15*T*T){ //rho < 1e-5
      return 5e-6/rho_c;
  }
  else if(p < -1.50992E-07 + 4.61521E-05*T - 1.89168E-13*T*T){ //rho < 1e-4
    return 5e-5/rho_c;
  }
  else if(p < -1.38128E-05 + 0.000461547*T - 1.47273E-11*T*T){ //rho < 1e-3
    return 5e-4/rho_c;
  }
  else if(p < -0.001212888 + 0.004617609*T - 1.22486E-09*T*T){ //rho < 1e-2
    return 5e-3/rho_c;
  }
  else if(p < -0.090803154 + 0.046320106*T - 8.03954E-08*T*T){ //rho < 1e-1
    return 5e-2/rho_c;
  }
  else if(p < -6.934621112 + 0.473323909*T - 5.27629E-06*T*T){ //rho < 1
    return 5e-2/rho_c;
  }
  else if(p < -532.9921491 + 5.433894672*T - 0.0003339*T*T){ //rho < 10
    return 5e-1/rho_c;
  }
  else if(p < -34931.47795 + 92.24848843*T - 0.015853248*T*T){ //rho < 100
    return 9/rho_c;
  }
  else if(p < -99354.97899 + 209.1578681*T - 0.032762406*T*T){ //rho < 200
    return 90/rho_c;
  }
  else if(p < -167178.3915 + 313.7993977*T - 0.033828951*T*T){ //rho < 300
    return 190/rho_c;
  }
  return delta_p_tau_rf(p, tau, 298/rho_c, 1);
}

s_real delta_p_tau_liq_guess_iapws95(s_real p, s_real tau){
  /* This uses approximate isochors to get a reasonable initial guess for the
  density solution for a vapor at a particular temperature and pressure.  If
  there are tricky regions where a really good guess is needed, a bracketing
  method can be used to refine the guess.*/
  s_real T=T_c/tau;

  if (p >= P_c && tau <= 1){
    return delta_p_tau_supercritical_guess_iapws95(p, tau);
  }
  // for the rest just use a guess between isochors
  else if(T > 646.0 && p > 21900.0){
    return delta_p_tau_rf(p, tau, 322.0/rho_c, 390.0/rho_c);
  }
  else if(T > 645.0 && p > 21500.0){
    return delta_p_tau_rf(p, tau, 330.0/rho_c, 450.0/rho_c);
  }
  return 1000.0/rho_c;
}

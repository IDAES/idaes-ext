#include"helmholtz_solve.h"
#include"swco2_param.h"

s_real p_sat_approx(s_real tau){
  s_real tt = 1 - 1/tau;
  return P_c*exp(tau*(
      -7.0602087*tt +
      1.9391218*pow(tt, 1.5) -
      1.6463597*pow(tt, 2.0) -
      3.2995634*pow(tt, 4.0)
    )
  );
}

s_real delta_sat_v_approx(s_real tau){
  s_real tt = 1 - 1/tau;
  return exp(
    -1.7074879*pow(tt, 0.340) -
    0.82274670*pow(tt, 0.5) -
    4.6008549*pow(tt, 1.0) -
    10.111178*pow(tt, 7.0/3.0) -
    29.742252*pow(tt, 14.0/3.0)
  );
}

s_real delta_sat_l_approx(s_real tau){
  s_real tt = 1 - 1/tau;
  return exp(
    1.9245108*pow(tt, 0.34) -
    0.62385555*pow(tt, 0.5) -
    0.32731127*pow(tt, 10.0/6.0) +
    0.39245142*pow(tt, 11.0/6.0)
  );
}

s_real delta_p_tau_supercritical_guess_swco2(s_real p, s_real tau){

  /* This uses approximate isochors to get a reasonable initial guess for the
  density solution for a vapor at a particular temperature and pressure.  If
  there are tricky regions where a really good guess is needed, a bracketing
  method can be used to refine the guess.*/
  s_real T=T_c/tau;

  if(p < 1.38648E-09 + 1.88923E-06*T + 1.14747E-14*T*T){ //rho < 1e-5
    return 5e-6/rho_c;
  }
  else if(p < 1.38648E-08 + 1.88923E-05*T + 1.14747E-13*T*T){ //rho < 1e-4
    return 5e-5/rho_c;
  }
  else if(p < -5.29756E-07 + 0.000188925*T - 5.85444E-14*T*T){ //rho < 1e-3
    return 5e-4/rho_c;
  }
  else if(p < -3.38279E-05 + 0.001889309*T - 2.8997E-11*T*T){ //rho < 1e-2
    return 5e-3/rho_c;
  }
  else if(p < -0.0033456 + 0.018898965*T - 2.60833E-09*T*T){ //rho < 1e-1
    return 5e-2/rho_c;
  }
  else if(p < -0.342149621 + 0.189613095*T - 2.86939E-07*T*T){ //rho < 1
    return 5e-1/rho_c;
  }
  else if(p < -34.14334862 + 1.958265958*T - 2.87239E-05*T*T){ //rho < 10
    return delta_p_tau_rf(p, tau, 0.5/rho_c, 15/rho_c);
  }
  else if(p < -2820.065655 + 24.25970886*T - 0.001748359*T*T){ //rho < 100
    return delta_p_tau_rf(p, tau, 5/rho_c, 105/rho_c);
  }
  else if(p < -10141.29813 + 56.98420928*T - 0.005166199*T*T){ //rho < 200
    return delta_p_tau_rf(p, tau, 95/rho_c, 205/rho_c);
  }
  else if(p < -21125.61421 + 97.22905662*T - 0.008799283*T*T){ //rho < 300
    return delta_p_tau_rf(p, tau, 195/rho_c, 305/rho_c);
  }
  else if(p < -35876.37569 + 146.8404454*T - 0.012537284*T*T){ //rho < 400
    return delta_p_tau_rf(p, tau, 295/rho_c, 405/rho_c);
  }
  else if(p < -54964.94489 + 209.5598894*T - 0.017302282*T*T){ //rho < 500
    return delta_p_tau_rf(p, tau, 395/rho_c, 505/rho_c);
  }
  else if(p < -79898.47714 + 293.2504079*T - 0.026823627*T*T){ //rho < 600
    return delta_p_tau_rf(p, tau, 495/rho_c, 605/rho_c);
  }
  else if(p < -112368.0377 + 408.529043*T - 0.046634364*T*T){ //rho < 700
    return delta_p_tau_rf(p, tau, 595/rho_c, 705/rho_c);
  }
  else if(p < -153093.0031 + 566.7194992*T - 0.082889317*T*T){ //rho < 800
    return delta_p_tau_rf(p, tau, 695/rho_c, 805/rho_c);
  }
  else if(p < -202040.1487 + 782.056835*T - 0.144729326*T*T){ //rho < 900
    return delta_p_tau_rf(p, tau, 795/rho_c, 905/rho_c);
  }
  return 1180/rho_c;
}

s_real delta_p_tau_vap_guess_swco2(s_real p, s_real tau){

  /* This uses approximate isochors to get a reasonable initial guess for the
  density solution for a vapor at a particular temperature and pressure.  If
  there are tricky regions where a really good guess is needed, a bracketing
  method can be used to refine the guess.*/
  s_real T=T_c/tau;

  if (p >= P_c && tau <= 1){
    return delta_p_tau_supercritical_guess_swco2(p, tau);
  }
  // for the rest just use a guess between isochors
  else if(p < 1.38648E-09 + 1.88923E-06*T + 1.14747E-14*T*T){ //rho < 1e-5
    return 5e-6/rho_c;
  }
  else if(p < 1.38648E-08 + 1.88923E-05*T + 1.14747E-13*T*T){ //rho < 1e-4
    return 5e-5/rho_c;
  }
  else if(p < -5.29756E-07 + 0.000188925*T - 5.85444E-14*T*T){ //rho < 1e-3
    return 5e-4/rho_c;
  }
  else if(p < -3.38279E-05 + 0.001889309*T - 2.8997E-11*T*T){ //rho < 1e-2
    return 5e-3/rho_c;
  }
  else if(p < -0.0033456 + 0.018898965*T - 2.60833E-09*T*T){ //rho < 1e-1
    return 5e-2/rho_c;
  }
  else if(p < -0.342149621 + 0.189613095*T - 2.86939E-07*T*T){ //rho < 1
    return 5e-1/rho_c;
  }
  else if(p < -34.14334862 + 1.958265958*T - 2.87239E-05*T*T){ //rho < 10
    return delta_p_tau_rf(p, tau, 0.5/rho_c, 15/rho_c);
  }
  else if(p < -2820.065655 + 24.25970886*T - 0.001748359*T*T){ //rho < 100
    return delta_p_tau_rf(p, tau, 5/rho_c, 105/rho_c);
  }
  else if(p < -10141.29813 + 56.98420928*T - 0.005166199*T*T){ //rho < 200
    return delta_p_tau_rf(p, tau, 95/rho_c, 205/rho_c);
  }
  else if(p < -21125.61421 + 97.22905662*T - 0.008799283*T*T){ //rho < 300
    return delta_p_tau_rf(p, tau, 195/rho_c, 305/rho_c);
  }
  else if(p < -35876.37569 + 146.8404454*T - 0.012537284*T*T){ //rho < 400
    return delta_p_tau_rf(p, tau, 295/rho_c, 405/rho_c);
  }
  return delta_p_tau_rf(p, tau, 395/rho_c, 1);
}

s_real delta_p_tau_liq_guess_swco2(s_real p, s_real tau){
  /* This uses approximate isochors to get a reasonable initial guess for the
  density solution for a vapor at a particular temperature and pressure.  If
  there are tricky regions where a really good guess is needed, a bracketing
  method can be used to refine the guess.*/
  s_real T=T_c/tau;

  if (p >= P_c && tau <= 1){
    return delta_p_tau_supercritical_guess_swco2(p, tau);
  }
  // for the rest just use a guess between isochors
  else if(p < -54964.94489 + 209.5598894*T - 0.017302282*T*T){ //rho < 500
    return delta_p_tau_rf(p, tau, 1, 505/rho_c);
  }
  else if(p < -79898.47714 + 293.2504079*T - 0.026823627*T*T){ //rho < 600
    return delta_p_tau_rf(p, tau, 495/rho_c, 605/rho_c);
  }
  else if(p < -112368.0377 + 408.529043*T - 0.046634364*T*T){ //rho < 700
    return delta_p_tau_rf(p, tau, 595/rho_c, 705/rho_c);
  }
  else if(p < -153093.0031 + 566.7194992*T - 0.082889317*T*T){ //rho < 800
    return delta_p_tau_rf(p, tau, 695/rho_c, 805/rho_c);
  }
  else if(p < -202040.1487 + 782.056835*T - 0.144729326*T*T){ //rho < 900
    return delta_p_tau_rf(p, tau, 795/rho_c, 905/rho_c);
  }
  else if(p < -256602.6793 + 1066.954274*T - 0.241727489*T*T){ //rho < 1000
    return delta_p_tau_rf(p, tau, 895/rho_c, 1005/rho_c);
  }
  else if(p < -312602.0898 + 1436.311945*T - 0.38689172*T*T){ //rho < 1100
    return delta_p_tau_rf(p, tau, 995/rho_c, 1105/rho_c);
  }
  return 1180/rho_c;
}

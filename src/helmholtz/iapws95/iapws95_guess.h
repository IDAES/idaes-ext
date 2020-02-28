
#ifndef _INCLUDE_IAPWS95_GUESS_H_
#define _INCLUDE_IAPWS95_GUESS_H_

s_real delta_p_tau_supercritical_guess_iapws95(s_real p, s_real tau);
s_real delta_p_tau_vap_guess_iapws95(s_real p, s_real tau);
s_real delta_p_tau_liq_guess_iapws95(s_real p, s_real tau);
s_real p_sat_iapws97(s_real tau);
s_real delta_sat_v_approx_iapws95(s_real tau);
s_real delta_sat_l_approx_iapws95(s_real tau);

 #endif

#include "props_up.h"
#include<iostream>
#include<stdlib.h>

prop_memo_table22 memo_table_enthalpy_up;
prop_memo_table22 memo_table_entropy_up;
prop_memo_table22 memo_table_gibbs_up;
prop_memo_table22 memo_table_helmholtz_up;
prop_memo_table22 memo_table_isochoric_heat_capacity_up;
prop_memo_table22 memo_table_isobaric_heat_capacity_up;
prop_memo_table22 memo_table_speed_of_sound_up;
prop_memo_table22 memo_table_specific_volume_up;
prop_memo_table22 memo_table_temperature_up;
prop_memo_table22 memo_table_vapor_fraction_up;

prop_memo_table22 memo_table_internal_energy_liq_up;
prop_memo_table22 memo_table_entropy_liq_up;
prop_memo_table22 memo_table_gibbs_liq_up;
prop_memo_table22 memo_table_helmholtz_liq_up;
prop_memo_table22 memo_table_isochoric_heat_capacity_liq_up;
prop_memo_table22 memo_table_isobaric_heat_capacity_liq_up;
prop_memo_table22 memo_table_speed_of_sound_liq_up;
prop_memo_table22 memo_table_specific_volume_liq_up;

prop_memo_table22 memo_table_internal_energy_vap_up;
prop_memo_table22 memo_table_entropy_vap_up;
prop_memo_table22 memo_table_gibbs_vap_up;
prop_memo_table22 memo_table_helmholtz_vap_up;
prop_memo_table22 memo_table_isochoric_heat_capacity_vap_up;
prop_memo_table22 memo_table_isobaric_heat_capacity_vap_up;
prop_memo_table22 memo_table_speed_of_sound_vap_up;
prop_memo_table22 memo_table_specific_volume_vap_up;

//T
void temperature_up(uint comp, double u, double p, f22_struct *out){ 
    f22_struct tau_vec = memo2_tau_up(comp, u, p);
    out->f = cdata[comp].T_star/tau_vec.f;
    out->f_1 = -cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_1;
    out->f_2 = -cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_2;
    out->f_11 = 2*cdata[comp].T_star/tau_vec.f/tau_vec.f/tau_vec.f*tau_vec.f_1*tau_vec.f_1 - cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_11;
    out->f_12 = 2*cdata[comp].T_star/tau_vec.f/tau_vec.f/tau_vec.f*tau_vec.f_1*tau_vec.f_2 - cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_12;
    out->f_22 = 2*cdata[comp].T_star/tau_vec.f/tau_vec.f/tau_vec.f*tau_vec.f_2*tau_vec.f_2 - cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_22;
}
MEMO2_FUNCTION(memo2_temperature_up, temperature_up, memo_table_temperature_up)

//vapor_fraction
void vapor_fraction_up(uint comp, double u, double p, f22_struct *out){ 
    f22_struct vf_vec = memo2_vf_up(comp, u, p);
    out->f = vf_vec.f;
    out->f_1 = vf_vec.f_1;
    out->f_2 = vf_vec.f_2;
    out->f_11 = vf_vec.f_11;
    out->f_12 = vf_vec.f_12;
    out->f_22 = vf_vec.f_22;
}
MEMO2_FUNCTION(memo2_vapor_fraction_sp, vapor_fraction_sp, memo_table_vapor_fraction_sp)

//u
PROP_UP_SINGLE_PHASE(enthalpy_vap_up, memo2_enthalpy, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(enthalpy_liq_up, memo2_enthalpy, memo2_delta_liquid)
PROP_UP_GENERAL(enthalpy_up, enthalpy_liq_up, enthalpy_vap_up)
MEMO2_FUNCTION(memo2_enthalpy_up, enthalpy_up, memo_table_enthalpy_up)

//s
PROP_UP_SINGLE_PHASE(entropy_vap_up, memo2_entropy, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(entropy_liq_up, memo2_entropy, memo2_delta_liquid)
PROP_UP_GENERAL(entropy_up, entropy_liq_up, entropy_vap_up)
MEMO2_FUNCTION(memo2_entropy_up, entropy_up, memo_table_entropy_up)

//g
PROP_UP_SINGLE_PHASE(gibbs_vap_up, memo2_gibbs, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(gibbs_liq_up, memo2_gibbs, memo2_delta_liquid)
PROP_UP_GENERAL(gibbs_up, gibbs_liq_up, gibbs_vap_up)
MEMO2_FUNCTION(memo2_gibbs_up, gibbs_up, memo_table_gibbs_up)

//f
PROP_UP_SINGLE_PHASE(helmholtz_vap_up, memo2_helmholtz, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(helmholtz_liq_up, memo2_helmholtz, memo2_delta_liquid)
PROP_UP_GENERAL(helmholtz_up, helmholtz_liq_up, helmholtz_vap_up)
MEMO2_FUNCTION(memo2_helmholtz_up, helmholtz_up, memo_table_helmholtz_up)

//cv
PROP_UP_SINGLE_PHASE(isochoric_heat_capacity_vap_up, memo2_isochoric_heat_capacity, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(isochoric_heat_capacity_liq_up, memo2_isochoric_heat_capacity, memo2_delta_liquid)
PROP_UP_GENERAL(isochoric_heat_capacity_up, isochoric_heat_capacity_liq_up, isochoric_heat_capacity_vap_up)
MEMO2_FUNCTION(memo2_isochoric_heat_capacity_up, isochoric_heat_capacity_up, memo_table_isochoric_heat_capacity_up)

//cp
PROP_UP_SINGLE_PHASE(isobaric_heat_capacity_vap_up, memo2_isobaric_heat_capacity, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(isobaric_heat_capacity_liq_up, memo2_isobaric_heat_capacity, memo2_delta_liquid)
PROP_UP_GENERAL(isobaric_heat_capacity_up, isobaric_heat_capacity_liq_up, isobaric_heat_capacity_vap_up)
MEMO2_FUNCTION(memo2_isobaric_heat_capacity_up, isobaric_heat_capacity_up, memo_table_isobaric_heat_capacity_up)

//w (doesn't really mean much in the two phase region, so use with care)
PROP_UP_SINGLE_PHASE(speed_of_sound_vap_up, memo2_speed_of_sound, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(speed_of_sound_liq_up, memo2_speed_of_sound, memo2_delta_liquid)
PROP_UP_GENERAL(speed_of_sound_up, speed_of_sound_liq_up, speed_of_sound_vap_up)
MEMO2_FUNCTION(memo2_speed_of_sound_up, speed_of_sound_up, memo_table_speed_of_sound_up)

//v
PROP_UP_SINGLE_PHASE(specific_volume_vap_up, memo2_specific_volume, memo2_delta_vapor)
PROP_UP_SINGLE_PHASE(specific_volume_liq_up, memo2_specific_volume, memo2_delta_liquid)
PROP_UP_GENERAL(specific_volume_up, specific_volume_liq_up, specific_volume_vap_up)
MEMO2_FUNCTION(memo2_specific_volume_up, specific_volume_up, memo_table_specific_volume_up)


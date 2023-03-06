#include "props_hp.h"
#include<iostream>
#include<stdlib.h>

prop_memo_table22 memo_table_internal_energy_hp;
prop_memo_table22 memo_table_entropy_hp;
prop_memo_table22 memo_table_gibbs_hp;
prop_memo_table22 memo_table_helmholtz_hp;
prop_memo_table22 memo_table_isochoric_heat_capacity_hp;
prop_memo_table22 memo_table_isobaric_heat_capacity_hp;
prop_memo_table22 memo_table_speed_of_sound_hp;
prop_memo_table22 memo_table_specific_volume_hp;
prop_memo_table22 memo_table_temperature_hp;
prop_memo_table22 memo_table_vapor_fraction_hp;

prop_memo_table22 memo_table_internal_energy_liq_hp;
prop_memo_table22 memo_table_entropy_liq_hp;
prop_memo_table22 memo_table_gibbs_liq_hp;
prop_memo_table22 memo_table_helmholtz_liq_hp;
prop_memo_table22 memo_table_isochoric_heat_capacity_liq_hp;
prop_memo_table22 memo_table_isobaric_heat_capacity_liq_hp;
prop_memo_table22 memo_table_speed_of_sound_liq_hp;
prop_memo_table22 memo_table_specific_volume_liq_hp;

prop_memo_table22 memo_table_internal_energy_vap_hp;
prop_memo_table22 memo_table_entropy_vap_hp;
prop_memo_table22 memo_table_gibbs_vap_hp;
prop_memo_table22 memo_table_helmholtz_vap_hp;
prop_memo_table22 memo_table_isochoric_heat_capacity_vap_hp;
prop_memo_table22 memo_table_isobaric_heat_capacity_vap_hp;
prop_memo_table22 memo_table_speed_of_sound_vap_hp;
prop_memo_table22 memo_table_specific_volume_vap_hp;

//T
void temperature_hp(uint comp, double h, double p, f22_struct *out){ 
    f22_struct tau_vec = memo2_tau_hp(comp, h, p);
    out->f = cdata[comp].T_star/tau_vec.f;
    out->f_1 = -cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_1;
    out->f_2 = -cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_2;
    out->f_11 = 2*cdata[comp].T_star/tau_vec.f/tau_vec.f/tau_vec.f*tau_vec.f_1*tau_vec.f_1 - cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_11;
    out->f_12 = 2*cdata[comp].T_star/tau_vec.f/tau_vec.f/tau_vec.f*tau_vec.f_1*tau_vec.f_2 - cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_12;
    out->f_22 = 2*cdata[comp].T_star/tau_vec.f/tau_vec.f/tau_vec.f*tau_vec.f_2*tau_vec.f_2 - cdata[comp].T_star/tau_vec.f/tau_vec.f*tau_vec.f_22;
}
MEMO2_FUNCTION(memo2_temperature_hp, temperature_hp, memo_table_temperature_hp)

//vapor_fraction
void vapor_fraction_hp(uint comp, double h, double p, f22_struct *out){ 
    f22_struct vf_vec = memo2_vf_hp(comp, h, p);
    out->f = vf_vec.f;
    out->f_1 = vf_vec.f_1;
    out->f_2 = vf_vec.f_2;
    out->f_11 = vf_vec.f_11;
    out->f_12 = vf_vec.f_12;
    out->f_22 = vf_vec.f_22;
}
MEMO2_FUNCTION(memo2_vapor_fraction_hp, vapor_fraction_hp, memo_table_vapor_fraction_hp)

//u
PROP_HP_SINGLE_PHASE(internal_energy_vap_hp, memo2_internal_energy, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(internal_energy_liq_hp, memo2_internal_energy, memo2_delta_liquid)
PROP_HP_GENERAL(internal_energy_hp, internal_energy_liq_hp, internal_energy_vap_hp)
MEMO2_FUNCTION(memo2_internal_energy_hp, internal_energy_hp, memo_table_internal_energy_hp)
MEMO2_FUNCTION(memo2_internal_energy_vap_hp, internal_energy_vap_hp, memo_table_internal_energy_vap_hp)
MEMO2_FUNCTION(memo2_internal_energy_liq_hp, internal_energy_liq_hp, memo_table_internal_energy_liq_hp)


//s
PROP_HP_SINGLE_PHASE(entropy_vap_hp, memo2_entropy, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(entropy_liq_hp, memo2_entropy, memo2_delta_liquid)
PROP_HP_GENERAL(entropy_hp, entropy_liq_hp, entropy_vap_hp)
MEMO2_FUNCTION(memo2_entropy_hp, entropy_hp, memo_table_entropy_hp)
MEMO2_FUNCTION(memo2_entropy_vap_hp, entropy_vap_hp, memo_table_entropy_vap_hp)
MEMO2_FUNCTION(memo2_entropy_liq_hp, entropy_liq_hp, memo_table_entropy_liq_hp)

//g
PROP_HP_SINGLE_PHASE(gibbs_vap_hp, memo2_gibbs, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(gibbs_liq_hp, memo2_gibbs, memo2_delta_liquid)
PROP_HP_GENERAL(gibbs_hp, gibbs_liq_hp, gibbs_vap_hp)
MEMO2_FUNCTION(memo2_gibbs_hp, gibbs_hp, memo_table_gibbs_hp)
MEMO2_FUNCTION(memo2_gibbs_vap_hp, gibbs_vap_hp, memo_table_gibbs_vap_hp)
MEMO2_FUNCTION(memo2_gibbs_liq_hp, gibbs_liq_hp, memo_table_gibbs_liq_hp)


//f
PROP_HP_SINGLE_PHASE(helmholtz_vap_hp, memo2_helmholtz, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(helmholtz_liq_hp, memo2_helmholtz, memo2_delta_liquid)
PROP_HP_GENERAL(helmholtz_hp, helmholtz_liq_hp, helmholtz_vap_hp)
MEMO2_FUNCTION(memo2_helmholtz_hp, helmholtz_hp, memo_table_helmholtz_hp)
MEMO2_FUNCTION(memo2_helmholtz_vap_hp, helmholtz_vap_hp, memo_table_helmholtz_vap_hp)
MEMO2_FUNCTION(memo2_helmholtz_liq_hp, helmholtz_liq_hp, memo_table_helmholtz_liq_hp)

//cv
PROP_HP_SINGLE_PHASE(isochoric_heat_capacity_vap_hp, memo2_isochoric_heat_capacity, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(isochoric_heat_capacity_liq_hp, memo2_isochoric_heat_capacity, memo2_delta_liquid)
PROP_HP_GENERAL(isochoric_heat_capacity_hp, isochoric_heat_capacity_liq_hp, isochoric_heat_capacity_vap_hp)
MEMO2_FUNCTION(memo2_isochoric_heat_capacity_hp, isochoric_heat_capacity_hp, memo_table_isochoric_heat_capacity_hp)
MEMO2_FUNCTION(memo2_isochoric_heat_capacity_vap_hp, isochoric_heat_capacity_vap_hp, memo_table_isochoric_heat_capacity_vap_hp)
MEMO2_FUNCTION(memo2_isochoric_heat_capacity_liq_hp, isochoric_heat_capacity_liq_hp, memo_table_isochoric_heat_capacity_liq_hp)

//cp
PROP_HP_SINGLE_PHASE(isobaric_heat_capacity_vap_hp, memo2_isobaric_heat_capacity, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(isobaric_heat_capacity_liq_hp, memo2_isobaric_heat_capacity, memo2_delta_liquid)
PROP_HP_GENERAL(isobaric_heat_capacity_hp, isobaric_heat_capacity_liq_hp, isobaric_heat_capacity_vap_hp)
MEMO2_FUNCTION(memo2_isobaric_heat_capacity_hp, isobaric_heat_capacity_hp, memo_table_isobaric_heat_capacity_hp)
MEMO2_FUNCTION(memo2_isobaric_heat_capacity_vap_hp, isobaric_heat_capacity_vap_hp, memo_table_isobaric_heat_capacity_vap_hp)
MEMO2_FUNCTION(memo2_isobaric_heat_capacity_liq_hp, isobaric_heat_capacity_liq_hp, memo_table_isobaric_heat_capacity_liq_hp)

//w (doesn't really mean much in the two phase region, so use with care)
PROP_HP_SINGLE_PHASE(speed_of_sound_vap_hp, memo2_speed_of_sound, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(speed_of_sound_liq_hp, memo2_speed_of_sound, memo2_delta_liquid)
PROP_HP_GENERAL(speed_of_sound_hp, speed_of_sound_liq_hp, speed_of_sound_vap_hp)
MEMO2_FUNCTION(memo2_speed_of_sound_hp, speed_of_sound_hp, memo_table_speed_of_sound_hp)
MEMO2_FUNCTION(memo2_speed_of_sound_vap_hp, speed_of_sound_vap_hp, memo_table_speed_of_sound_vap_hp)
MEMO2_FUNCTION(memo2_speed_of_sound_liq_hp, speed_of_sound_liq_hp, memo_table_speed_of_sound_liq_hp)

//v
PROP_HP_SINGLE_PHASE(specific_volume_vap_hp, memo2_specific_volume, memo2_delta_vapor)
PROP_HP_SINGLE_PHASE(specific_volume_liq_hp, memo2_specific_volume, memo2_delta_liquid)
PROP_HP_GENERAL(specific_volume_hp, specific_volume_liq_hp, specific_volume_vap_hp)
MEMO2_FUNCTION(memo2_specific_volume_hp, specific_volume_hp, memo_table_specific_volume_hp)
MEMO2_FUNCTION(memo2_specific_volume_vap_hp, specific_volume_vap_hp, memo_table_specific_volume_vap_hp)
MEMO2_FUNCTION(memo2_specific_volume_liq_hp, specific_volume_liq_hp, memo_table_specific_volume_liq_hp)

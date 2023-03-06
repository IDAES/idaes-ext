#include "config.h"
#include "props.h"
#include "delta.h"
#include "state.h"
#include "phi.h"

f22_struct memo2_temperature_sp(uint comp, double s, double p);
f22_struct memo2_vapor_fraction_sp(uint comp, double s, double p);

f22_struct memo2_internal_energy_sp(uint comp, double s, double p);
f22_struct memo2_enthalpy_sp(uint comp, double s, double p);
f22_struct memo2_gibbs_sp(uint comp, double s, double p);
f22_struct memo2_helmholtz_sp(uint comp, double s, double p);
f22_struct memo2_isochoric_heat_capacity_sp(uint comp, double s, double p);
f22_struct memo2_isobaric_heat_capacity_sp(uint comp, double s, double p);
f22_struct memo2_speed_of_sound_sp(uint comp, double s, double p);
f22_struct memo2_specific_volume_sp(uint comp, double s, double p);
f22_struct memo2_viscosity_sp(uint comp, double s, double p);
f22_struct memo2_thermal_conductivity_sp(uint comp, double s, double p);
f22_struct memo2_surface_tension_sp(uint comp, double s, double p);

f22_struct memo2_internal_energy_vap_sp(uint comp, double s, double p);
f22_struct memo2_enthalpy_vap_sp(uint comp, double s, double p);
f22_struct memo2_gibbs_vap_sp(uint comp, double s, double p);
f22_struct memo2_helmholtz_sp(uint comp, double s, double p);
f22_struct memo2_isochoric_heat_capacity_vap_sp(uint comp, double s, double p);
f22_struct memo2_isobaric_heat_capacity_vap_sp(uint comp, double s, double p);
f22_struct memo2_speed_of_sound_vap_sp(uint comp, double s, double p);
f22_struct memo2_specific_volume_vap_sp(uint comp, double s, double p);
f22_struct memo2_viscosity_vap_sp(uint comp, double s, double p);
f22_struct memo2_thermal_conductivity_vap_sp(uint comp, double s, double p);
f22_struct memo2_surface_tension_vap_sp(uint comp, double s, double p);

f22_struct memo2_internal_energy_liq_sp(uint comp, double s, double p);
f22_struct memo2_enthalpy_liq_sp(uint comp, double s, double p);
f22_struct memo2_gibbs_liq_sp(uint comp, double s, double p);
f22_struct memo2_helmholtz_liq_sp(uint comp, double s, double p);
f22_struct memo2_isochoric_heat_capacity_liq_sp(uint comp, double s, double p);
f22_struct memo2_isobaric_heat_capacity_liq_sp(uint comp, double s, double p);
f22_struct memo2_speed_of_sound_liq_sp(uint comp, double s, double p);
f22_struct memo2_specific_volume_liq_sp(uint comp, double s, double p);
f22_struct memo2_viscosity_liq_sp(uint comp, double s, double p);
f22_struct memo2_thermal_conductivity_liq_sp(uint comp, double s, double p);
f22_struct memo2_surface_tension_liq_sp(uint comp, double s, double p);


#define PROP_SP_SINGLE_PHASE(new_func, prop_func, delta_func) \
void new_func(uint comp, double s, double p, f22_struct *out){ \
    f22_struct tau_vec, delta_vec, prop_vec; \
    tau_vec = memo2_tau_sp(comp, s, p); \
    delta_vec = delta_func(comp, p, tau_vec.f); \
    prop_vec = prop_func(comp, delta_vec.f, tau_vec.f); \
    double dddh = delta_vec.f_2 * tau_vec.f_1; \
    double dddp = delta_vec.f_1 + delta_vec.f_2 * tau_vec.f_2; \
    double d2ddh2 = delta_vec.f_22 * tau_vec.f_1 * tau_vec.f_1 + delta_vec.f_2 * tau_vec.f_11; \
    double d2ddp2 = delta_vec.f_11 + delta_vec.f_22 * tau_vec.f_2 * tau_vec.f_2 + delta_vec.f_2 * tau_vec.f_22 + 2 * delta_vec.f_12 * tau_vec.f_2; \
    double d2ddpdh = delta_vec.f_22 * tau_vec.f_1 * tau_vec.f_2 + delta_vec.f_2 * tau_vec.f_12 + delta_vec.f_12 * tau_vec.f_1; \
    out->f = prop_vec.f; \
    out->f_1 = prop_vec.f_1 * dddh \
        + prop_vec.f_2 * tau_vec.f_1; \
    out->f_2 = prop_vec.f_1 * dddp \
        + prop_vec.f_2 * tau_vec.f_2; \
    out->f_11 = \
        prop_vec.f_11 * dddh * dddh \
        + 2 * prop_vec.f_12 * dddh * tau_vec.f_1 \
        + prop_vec.f_1 * d2ddh2 \
        + prop_vec.f_22 * tau_vec.f_1 * tau_vec.f_1 \
        + prop_vec.f_2 * tau_vec.f_11; \
    out->f_22 = \
        prop_vec.f_11 * dddp * dddp \
        + 2 * prop_vec.f_12 * dddp * tau_vec.f_2 \
        + prop_vec.f_1 * d2ddp2 \
        + prop_vec.f_22 * tau_vec.f_2 * tau_vec.f_2 \
        + prop_vec.f_2 * tau_vec.f_22; \
    out->f_12 = \
        prop_vec.f_11 * dddp * dddh \
        + prop_vec.f_12 * dddh * tau_vec.f_2 \
        + prop_vec.f_12 * dddp * tau_vec.f_1 \
        + prop_vec.f_1 * d2ddpdh \
        + prop_vec.f_22 * tau_vec.f_1 * tau_vec.f_2 \
        + prop_vec.f_2 * tau_vec.f_12; \
} 

#define PROP_SP_GENERAL(new_func, prop_liq_func, prop_vap_func) \
void new_func(uint comp, double s, double p, f22_struct *out){ \
    f22_struct tau_vec, delta_vec, vf_vec, prop_vec; \
    vf_vec = memo2_vf_sp(comp, s, p); \
    if(vf_vec.f == 0.0){ \
        prop_liq_func(comp, s, p, out); \
    } \
    else if(vf_vec.f == 1.0){ \
        prop_vap_func(comp, s, p, out); \
    } \
    else{ \
        f22_struct out_vap, out_liq; \
        prop_vap_func(comp, s, p, &out_vap); \
        prop_liq_func(comp, s, p, &out_liq); \
        out->f = vf_vec.f * out_vap.f + (1 - vf_vec.f) * out_liq.f; \
        out->f_1 = vf_vec.f_1 * out_vap.f \
            + vf_vec.f * out_vap.f_1 \
            - vf_vec.f_1 * out_liq.f \
            + (1 - vf_vec.f) * out_liq.f_1; \
        out->f_2 = vf_vec.f_2 * out_vap.f \
            + vf_vec.f * out_vap.f_2 \
            - vf_vec.f_2 * out_liq.f \
            + (1 - vf_vec.f) * out_liq.f_2; \
        out->f_11 = vf_vec.f_11 * out_vap.f \
            + 2 * vf_vec.f_1 * out_vap.f_1 \
            + vf_vec.f * out_vap.f_11 \
            - vf_vec.f_11 * out_liq.f \
            - 2 * vf_vec.f_1 * out_liq.f_1 \
            - (1 - vf_vec.f) * out_liq.f_11; \
        out->f_12 = vf_vec.f_12 * out_vap.f \
            + vf_vec.f_1 * out_vap.f_2 \
            + vf_vec.f_2 * out_vap.f_1 \
            + vf_vec.f * out_vap.f_12 \
            - vf_vec.f_12 * out_liq.f \
            - vf_vec.f_1 * out_liq.f_2 \
            - vf_vec.f_2 * out_liq.f_1 \
            - (1 - vf_vec.f) * out_liq.f_12; \
        out->f_22 = vf_vec.f_22 * out_vap.f \
            + 2 * vf_vec.f_2 * out_vap.f_2 \
            + vf_vec.f * out_vap.f_22 \
            - vf_vec.f_22 * out_liq.f \
            - 2 * vf_vec.f_2 * out_liq.f_2 \
            + (1 - vf_vec.f) * out_liq.f_22; \
    } \
}
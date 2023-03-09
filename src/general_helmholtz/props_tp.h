#include "config.h"
#include "props.h"
#include "delta.h"
#include "state.h"
#include "phi.h"

void enthalpy_liq_tp(uint comp, double t, double p, f22_struct *out);
void entropy_liq_tp(uint comp, double t, double p, f22_struct *out);
void internal_energy_liq_tp(uint comp, double t, double p, f22_struct *out);
void gibbs_liq_tp(uint comp, double t, double p, f22_struct *out);
void helmholtz_liq_tp(uint comp, double t, double p, f22_struct *out);
void isochoric_heat_capacity_liq_tp(uint comp, double t, double p, f22_struct *out);
void isobaric_heat_capacity_liq_tp(uint comp, double t, double p, f22_struct *out);
void speed_of_sound_liq_tp(uint comp, double t, double p, f22_struct *out);
void specific_volume_liq_tp(uint comp, double t, double p, f22_struct *out);
void viscosity_liq_tp(uint comp, double t, double p, f22_struct *out);
void thermal_conductivity_liq_tp(uint comp, double t, double p, f22_struct *out);
void surface_tension_liq_tp(uint comp, double t, double p, f22_struct *out);

void enthalpy_vap_tp(uint comp, double t, double p, f22_struct *out);
void entropy_vap_tp(uint comp, double t, double p, f22_struct *out);
void internal_energy_vap_tp(uint comp, double t, double p, f22_struct *out);
void gibbs_vap_tp(uint comp, double t, double p, f22_struct *out);
void helmholtz_vap_tp(uint comp, double t, double p, f22_struct *out);
void isochoric_heat_capacity_vap_tp(uint comp, double t, double p, f22_struct *out);
void isobaric_heat_capacity_vap_tp(uint comp, double t, double p, f22_struct *out);
void speed_of_sound_vap_tp(uint comp, double t, double p, f22_struct *out);
void specific_volume_vap_tp(uint comp, double t, double p, f22_struct *out);
void viscosity_vap_tp(uint comp, double t, double p, f22_struct *out);
void thermal_conductivity_vap_tp(uint comp, double t, double p, f22_struct *out);
void surface_tension_vap_tp(uint comp, double t, double p, f22_struct *out);

f22_struct memo2_enthalpy_liq_tp(uint comp, double t, double p);
f22_struct memo2_internal_energy_liq_tp(uint comp, double t, double p);
f22_struct memo2_entropy_liq_tp(uint comp, double t, double p);
f22_struct memo2_gibbs_liq_tp(uint comp, double t, double p);
f22_struct memo2_helmholtz_liq_tp(uint comp, double t, double p);
f22_struct memo2_isochoric_heat_capacity_liq_tp(uint comp, double t, double p);
f22_struct memo2_isobaric_heat_capacity_liq_tp(uint comp, double t, double p);
f22_struct memo2_speed_of_sound_liq_tp(uint comp, double t, double p);
f22_struct memo2_specific_volume_liq_tp(uint comp, double t, double p);
f22_struct memo2_viscosity_liq_tp(uint comp, double t, double p);
f22_struct memo2_thermal_conductivity_liq_tp(uint comp, double t, double p);
f22_struct memo2_surface_tension_liq_tp(uint comp, double t, double p);

f22_struct memo2_enthalpy_vap_tp(uint comp, double t, double p);
f22_struct memo2_internal_energy_vap_tp(uint comp, double t, double p);
f22_struct memo2_entropy_vap_tp(uint comp, double t, double p);
f22_struct memo2_gibbs_vap_tp(uint comp, double t, double p);
f22_struct memo2_helmholtz_vap_tp(uint comp, double t, double p);
f22_struct memo2_isochoric_heat_capacity_vap_tp(uint comp, double t, double p);
f22_struct memo2_isobaric_heat_capacity_vap_tp(uint comp, double t, double p);
f22_struct memo2_speed_of_sound_vap_tp(uint comp, double t, double p);
f22_struct memo2_specific_volume_vap_tp(uint comp, double t, double p);
f22_struct memo2_viscosity_vap_tp(uint comp, double t, double p);
f22_struct memo2_thermal_conductivity_vap_tp(uint comp, double t, double p);
f22_struct memo2_surface_tension_vap_tp(uint comp, double t, double p);

#define PROP_TP_SINGLE_PHASE(new_func, prop_func, delta_func) \
void new_func(uint comp, double t, double p, f22_struct *out){ \
    f22_struct tau_vec, delta_vec, prop_vec; \
    tau_vec.f = cdata[comp].T_star/t; \
    tau_vec.f_1 = -cdata[comp].T_star/t/t; \
    tau_vec.f_11 = 2*cdata[comp].T_star/t/t/t; \
    tau_vec.f_2 = 0.0; \
    tau_vec.f_12 = 0.0; \
    tau_vec.f_22 = 0.0; \
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
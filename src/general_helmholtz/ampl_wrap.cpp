/*-------------------------------------------------------------------------------+
| The Institute for the Design of Advanced Energy Systems Integrated Platform    |
| Framework (IDAES IP) was produced under the DOE Institute for the              |
| Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021      |
| by the software owners: The Regents of the University of California, through   |
| Lawrence Berkeley National Laboratory,  National Technology & Engineering      |
| Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University |
| Research Corporation, et al.  All rights reserved.                             |
|                                                                                |
| Please see the files COPYRIGHT.md and LICENSE.md for full copyright and        |
| license information.                                                           |
+-------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------
 Provide the AMPL user-defined function wrapper for property functions

 Author: John Eslick
 File: ampl_wrap.cpp
--------------------------------------------------------------------------------*/

#include "ampl_wrap.h"
#include "config.h"
#include "props.h"
#include "props_hp.h"
#include "props_sp.h"
#include "props_up.h"
#include "state.h"
#include "delta.h"
#include "sat.h"

ASL_WRAP_FUNC_2ARG(p, memo2_pressure)                    // p(comp, delta, tau) [kPa]
ASL_WRAP_FUNC_2ARG(u, memo2_internal_energy)             // u(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(s, memo2_entropy)                     // s(comp, delta, tau) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(h, memo2_enthalpy)                    // h(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(g, memo2_gibbs)                       // g(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(f, memo2_helmholtz)                   // f(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(cv, memo2_isochoric_heat_capacity)    // cv(comp, delta, tau) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(cp, memo2_isobaric_heat_capacity)     // cp(comp, delta, tau) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(w, memo2_speed_of_sound)              // w(comp, delta, tau) [m/s]
ASL_WRAP_FUNC_2ARG(v, memo2_specific_volume)             // v(comp, delta, tau) [m3/kg]
ASL_WRAP_FUNC_2ARG(hvpt, memo2_enthalpy_vapor)           // hv(comp, pressure, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(hlpt, memo2_enthalpy_liquid)          // hl(comp, pressure, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(svpt, memo2_entropy_vapor)            // sv(comp, pressure, tau) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(slpt, memo2_entropy_liquid)           // sl(comp, pressure, tau) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(uvpt, memo2_internal_energy_vapor)    // uv(comp, pressure, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(ulpt, memo2_internal_energy_liquid)   // ul(comp, pressure, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(tau, memo2_tau_hp)                    // tau(comp, enthalpy, pressure) [none]
ASL_WRAP_FUNC_2ARG(taus, memo2_tau_sp)                   // tau(comp, entropy, pressure) [none]
ASL_WRAP_FUNC_2ARG(tauu, memo2_tau_up)                   // tau(comp, internal energy, pressure) [none]
ASL_WRAP_FUNC_2ARG(vf, memo2_vf_hp)                      // vf(comp, enthalpy, pressure) [none]
ASL_WRAP_FUNC_2ARG(vfs, memo2_vf_sp)                     // vf(comp, entropy, pressure) [none]
ASL_WRAP_FUNC_2ARG(vfu, memo2_vf_up)                     // vf(comp, internal energy, pressure) [none]
ASL_WRAP_FUNC_2ARG(delta_liq, memo2_delta_liquid)        // delta_liq(comp, pressure, tau) [none]
ASL_WRAP_FUNC_2ARG(delta_vap, memo2_delta_vapor)         // delta_vap(comp, pressure, tau) [none]

ASL_WRAP_FUNC_2ARG(phi0, memo2_phi_ideal)                // phi0(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phir, memo2_phi_resi)                 // phir(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phi0_d, memo2_phi_ideal_d)            // phi0_d(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phir_d, memo2_phi_resi_d)             // phir_d(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phi0_dd, memo2_phi_ideal_dd)          // phi0_dd(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phir_dd, memo2_phi_resi_dd)           // phir_dd(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phi0_t, memo2_phi_ideal_t)            // phi0_t(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phir_t, memo2_phi_resi_t)             // phir_t(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phi0_dt, memo2_phi_ideal_dt)          // phi0_dt(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phir_dt, memo2_phi_resi_dt)           // phir_dt(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phi0_tt, memo2_phi_ideal_tt)          // phi0_tt(comp, delta, tau) [none]
ASL_WRAP_FUNC_2ARG(phir_tt, memo2_phi_resi_tt)           // phir_tt(comp, delta, tau) [none]

// Saturation curve functions
ASL_WRAP_FUNC_1ARG(tau_sat, sat_tau)                     // tau_sat(comp, pressure) [none]
ASL_WRAP_FUNC_1ARG(p_sat, sat_p)                         // p_sat(comp, pressure) [kPa]
ASL_WRAP_FUNC_1ARG(delta_sat_v, sat_delta_v)             // delta_sat_v(comp, pressure) [none]
ASL_WRAP_FUNC_1ARG(delta_sat_l, sat_delta_l)             // delta_sat_l(comp, pressure) [none]

// General (liquid, vapor, two-phase) property functions of (h, p)
ASL_WRAP_FUNC_2ARG(tau_hp, memo2_tau_hp)                    // tau(comp, h, p) [none]
ASL_WRAP_FUNC_2ARG(vf_hp, memo2_vf_hp)                      // vf(comp, h, p) [none]
ASL_WRAP_FUNC_2ARG(u_hp, memo2_internal_energy_hp)          // u(comp, h, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(s_hp, memo2_entropy_hp)                  // s(comp, h, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(g_hp, memo2_gibbs_hp)                    // g(comp, h, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(f_hp, memo2_helmholtz_hp)                // f(comp, h, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(cv_hp, memo2_isochoric_heat_capacity_hp) // cv(comp, h, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(cp_hp, memo2_isobaric_heat_capacity_hp)  // cp(comp, h, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(w_hp, memo2_speed_of_sound_hp)           // w(comp, h, p) [m/s]
ASL_WRAP_FUNC_2ARG(v_hp, memo2_specific_volume_hp)          // v(comp, h, p) [m3/kg]

// General (liquid, vapor, two-phase) property functions of (s, p)
ASL_WRAP_FUNC_2ARG(tau_sp, memo2_tau_sp)                    // tau(comp, s, p) [none]
ASL_WRAP_FUNC_2ARG(vf_sp, memo2_vf_sp)                      // vf(comp, s, p) [none]
ASL_WRAP_FUNC_2ARG(u_sp, memo2_internal_energy_sp)          // u(comp, s, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(h_sp, memo2_enthalpy_sp)                 // h(comp, s, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(g_sp, memo2_gibbs_sp)                    // g(comp, s, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(f_sp, memo2_helmholtz_sp)                // f(comp, s, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(cv_sp, memo2_isochoric_heat_capacity_sp) // cv(comp, s, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(cp_sp, memo2_isobaric_heat_capacity_sp)  // cp(comp, s, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(w_sp, memo2_speed_of_sound_sp)           // w(comp, s, p) [m/s]
ASL_WRAP_FUNC_2ARG(v_sp, memo2_specific_volume_sp)          // v(comp, s, p) [m3/kg]

// General (liquid, vapor, two-phase) property functions of (u, p)
ASL_WRAP_FUNC_2ARG(tau_up, memo2_tau_up)                    // tau(comp, u, p) [none]
ASL_WRAP_FUNC_2ARG(vf_up, memo2_vf_up)                      // vf(comp, u, p) [none]
ASL_WRAP_FUNC_2ARG(s_up, memo2_entropy_up)                  // s(comp, u, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(h_up, memo2_enthalpy_up)                 // h(comp, u, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(g_up, memo2_gibbs_up)                    // g(comp, u, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(f_up, memo2_helmholtz_up)                // f(comp, u, p) [kJ/kg]
ASL_WRAP_FUNC_2ARG(cv_up, memo2_isochoric_heat_capacity_up) // cv(comp, u, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(cp_up, memo2_isobaric_heat_capacity_up)  // cp(comp, u, p) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(w_up, memo2_speed_of_sound_up)           // w(comp, u, p) [m/s]
ASL_WRAP_FUNC_2ARG(v_up, memo2_specific_volume_up)          // v(comp, u, p) [m3/kg]

// Some parameters to make it easier to sync Pyomo parameters with external functions
// these take the component name as an argument
ASL_WRAP_FUNC_0ARG(mw, MW)              // Critical Pressure     [g/mol]
ASL_WRAP_FUNC_0ARG(t_star, T_star)      // Temperature to calculate tau [K]
ASL_WRAP_FUNC_0ARG(rho_star, rho_star)  // Desity to calculate delta [kg/m^3]
ASL_WRAP_FUNC_0ARG(pc, Pc)              // Critical Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tc, Tc)              // Critical Temperature  [K]
ASL_WRAP_FUNC_0ARG(rhoc, rhoc)          // Critical Density      [kg/m^3]
ASL_WRAP_FUNC_0ARG(pt, Pt)              // Critical Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tt, Tt)              // Critical Temperature  [K]
ASL_WRAP_FUNC_0ARG(rhot_v, rhot_v)      // Critical Density      [kg/m^3]
ASL_WRAP_FUNC_0ARG(rhot_l, rhot_l)      // Critical Density      [kg/m^3]
ASL_WRAP_FUNC_0ARG(sgc, R)              // Specific gas constant [kJ/kg/K] or [kPa m^3/kg/K]
ASL_WRAP_FUNC_0ARG(pmin, P_min)         // Minimum Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tmin, T_min)         // Minumum Temperature  [K]
ASL_WRAP_FUNC_0ARG(pmax, P_max)         // Minimum Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tmax, T_max)         // Minumum Temperature  [K]

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info */
    int typ = FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS;
    addfunc("p", (rfunc)p, typ, -4, NULL);
    addfunc("u", (rfunc)u, typ, -4, NULL);
    addfunc("s", (rfunc)s, typ, -4, NULL);
    addfunc("h", (rfunc)h, typ, -4, NULL);
    addfunc("g", (rfunc)g, typ, -4, NULL);
    addfunc("f", (rfunc)f, typ, -4, NULL);
    addfunc("cv", (rfunc)cv, typ, -4, NULL);
    addfunc("cp", (rfunc)cp, typ, -4, NULL);
    addfunc("w", (rfunc)w, typ, -4, NULL);
    addfunc("v", (rfunc)v, typ, -4, NULL);

    addfunc("tau_hp", (rfunc)tau_hp, typ, -4, NULL);
    addfunc("vf_hp", (rfunc)vf_hp, typ, -4, NULL);
    addfunc("u_hp", (rfunc)u_hp, typ, -4, NULL);
    addfunc("s_hp", (rfunc)s_hp, typ, -4, NULL);
    addfunc("g_hp", (rfunc)g_hp, typ, -4, NULL);
    addfunc("f_hp", (rfunc)f_hp, typ, -4, NULL);
    addfunc("cv_hp", (rfunc)cv_hp, typ, -4, NULL);
    addfunc("cp_hp", (rfunc)cp_hp, typ, -4, NULL);
    addfunc("w_hp", (rfunc)w_hp, typ, -4, NULL);
    addfunc("v_hp", (rfunc)v_hp, typ, -4, NULL);

    addfunc("tau_sp", (rfunc)tau_sp, typ, -4, NULL);
    addfunc("vf_sp", (rfunc)vf_sp, typ, -4, NULL);
    addfunc("u_sp", (rfunc)u_sp, typ, -4, NULL);
    addfunc("h_sp", (rfunc)h_sp, typ, -4, NULL);
    addfunc("g_sp", (rfunc)g_sp, typ, -4, NULL);
    addfunc("f_sp", (rfunc)f_sp, typ, -4, NULL);
    addfunc("cv_sp", (rfunc)cv_sp, typ, -4, NULL);
    addfunc("cp_sp", (rfunc)cp_sp, typ, -4, NULL);
    addfunc("w_sp", (rfunc)w_sp, typ, -4, NULL);
    addfunc("v_sp", (rfunc)v_sp, typ, -4, NULL);

    addfunc("tau_up", (rfunc)tau_up, typ, -4, NULL);
    addfunc("vf_up", (rfunc)vf_up, typ, -4, NULL);
    addfunc("s_up", (rfunc)s_up, typ, -4, NULL);
    addfunc("h_up", (rfunc)h_up, typ, -4, NULL);
    addfunc("g_up", (rfunc)g_up, typ, -4, NULL);
    addfunc("f_up", (rfunc)f_up, typ, -4, NULL);
    addfunc("cv_up", (rfunc)cv_up, typ, -4, NULL);
    addfunc("cp_up", (rfunc)cp_up, typ, -4, NULL);
    addfunc("w_up", (rfunc)w_up, typ, -4, NULL);
    addfunc("v_up", (rfunc)v_up, typ, -4, NULL);

    addfunc("hvpt", (rfunc)hvpt, typ, -4, NULL);
    addfunc("hlpt", (rfunc)hlpt, typ, -4, NULL);
    addfunc("svpt", (rfunc)svpt, typ, -4, NULL);
    addfunc("slpt", (rfunc)slpt, typ, -4, NULL);
    addfunc("uvpt", (rfunc)uvpt, typ, -4, NULL);
    addfunc("ulpt", (rfunc)ulpt, typ, -4, NULL);
    addfunc("tau", (rfunc)tau, typ, -4, NULL);
    addfunc("taus", (rfunc)taus, typ, -4, NULL);
    addfunc("tauu", (rfunc)tauu, typ, -4, NULL);
    addfunc("vf", (rfunc)vf, typ, -4, NULL);
    addfunc("vfs", (rfunc)vfs, typ, -4, NULL);
    addfunc("vfu", (rfunc)vfu, typ, -4, NULL);
    addfunc("delta_liq", (rfunc)delta_liq, typ, -4, NULL);
    addfunc("delta_vap", (rfunc)delta_vap, typ, -4, NULL);
    // phi and derivatives for calculating more thermo properties.
    addfunc("phi0", (rfunc)phi0, typ, -4, NULL);
    addfunc("phir", (rfunc)phir, typ, -4, NULL);
    addfunc("phi0_d", (rfunc)phi0_d, typ, -4, NULL);
    addfunc("phir_d", (rfunc)phir_d, typ, -4, NULL);
    addfunc("phi0_t", (rfunc)phi0_t, typ, -4, NULL);
    addfunc("phir_t", (rfunc)phir_t, typ, -4, NULL);
    addfunc("phi0_dd", (rfunc)phi0_dd, typ, -4, NULL);
    addfunc("phir_dd", (rfunc)phir_dd, typ, -4, NULL);
    addfunc("phi0_dt", (rfunc)phi0_dt, typ, -4, NULL);
    addfunc("phir_dt", (rfunc)phir_dt, typ, -4, NULL);
    addfunc("phi0_tt", (rfunc)phi0_tt, typ, -4, NULL);
    addfunc("phir_tt", (rfunc)phir_tt, typ, -4, NULL);
    // Unary functions for sat curve
    addfunc("delta_sat_l", (rfunc)delta_sat_l, typ, -3, NULL);
    addfunc("delta_sat_v", (rfunc)delta_sat_v, typ, -3, NULL);
    addfunc("p_sat", (rfunc)p_sat, typ, -3, NULL);
    addfunc("tau_sat", (rfunc)tau_sat, typ, -3, NULL);
    // Parameters
    addfunc("mw", (rfunc)mw, typ, -2, NULL);
    addfunc("t_star", (rfunc)t_star, typ, -2, NULL);
    addfunc("rho_star", (rfunc)rho_star, typ, -2, NULL);
    addfunc("pc", (rfunc)pc, typ, -2, NULL);
    addfunc("tc", (rfunc)tc, typ, -2, NULL);
    addfunc("rhoc", (rfunc)rhoc, typ, -2, NULL);
    addfunc("pt", (rfunc)pt, typ, -2, NULL);
    addfunc("tt", (rfunc)tt, typ, -2, NULL);
    addfunc("rhot_v", (rfunc)rhot_l, typ, -2, NULL);
    addfunc("rhot_l", (rfunc)rhot_v, typ, -2, NULL);
    addfunc("sgc", (rfunc)sgc, typ, -2, NULL);
    addfunc("pmin", (rfunc)pmin, typ, -2, NULL);
    addfunc("tmin", (rfunc)tmin, typ, -2, NULL);
    addfunc("pmax", (rfunc)pmax, typ, -2, NULL);
    addfunc("tmax", (rfunc)tmax, typ, -2, NULL);
}

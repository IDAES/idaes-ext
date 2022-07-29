#import "ampl_wrap.h"
#import "config.h"
#import "props.h"
#import "state.h"
#import "delta.h"
#import "sat.h"
#import "param.h"

ASL_WRAP_FUNC_2ARG(p, memo2_pressure)                    // p(comp, delta, tau) [kPa]
ASL_WRAP_FUNC_2ARG(u, memo2_internal_energy)             // u(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(s, memo2_entropy)                     // s(comp, delta, tau) [kJ/kg/K]
ASL_WRAP_FUNC_2ARG(h, memo2_enthalpy)                    // h(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(g, memo2_gibbs)                       // g(comp, delta, tau) [kJ/kg]
ASL_WRAP_FUNC_2ARG(f, memo2_helmholtz)                   // f(comp, delta, tau) [kJ/kg]
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

// Some parameters to make it easier to sync Pyomo parameters with external functions
ASL_WRAP_FUNC_0ARG(mw, param::mw)          // Critical Pressure     [g/mol]
ASL_WRAP_FUNC_0ARG(pc, param::Pc)          // Critical Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tc, param::Tc)          // Critical Temperature  [K]
ASL_WRAP_FUNC_0ARG(rhoc, param::rhoc)      // Critical Density      [kg/m^3]
ASL_WRAP_FUNC_0ARG(pt, param::Pt)          // Critical Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tt, param::Tt)          // Critical Temperature  [K]
ASL_WRAP_FUNC_0ARG(rhot_v, param::rhot_v)  // Critical Density      [kg/m^3]
ASL_WRAP_FUNC_0ARG(rhot_l, param::rhot_l)  // Critical Density      [kg/m^3]
ASL_WRAP_FUNC_0ARG(sgc, param::R)          // Specific gas constant [kJ/kg/K] or [kPa m^3/kg/K]
ASL_WRAP_FUNC_0ARG(pmin, param::P_min)     // Minimum Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tmin, param::T_min)     // Minumum Temperature  [K]
ASL_WRAP_FUNC_0ARG(pmax, param::P_max)     // Minimum Pressure     [kPa]
ASL_WRAP_FUNC_0ARG(tmax, param::T_max)     // Minumum Temperature  [K]

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of real arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info */
    int typ = FUNCADD_REAL_VALUED;
    addfunc("p", (rfunc)p, typ, 2, NULL);
    addfunc("u", (rfunc)u, typ, 2, NULL);
    addfunc("s", (rfunc)s, typ, 2, NULL);
    addfunc("h", (rfunc)h, typ, 2, NULL);
    addfunc("g", (rfunc)g, typ, 2, NULL);
    addfunc("f", (rfunc)f, typ, 2, NULL);
    /*
    addfunc("cv", (rfunc)cv_EOS_TAG, typ, 2, NULL);
    addfunc("cp", (rfunc)cp_EOS_TAG, typ, 2, NULL);
    addfunc("w", (rfunc)w_EOS_TAG, typ, 2, NULL);
    */
    addfunc("hvpt", (rfunc)hvpt, typ, 2, NULL);
    addfunc("hlpt", (rfunc)hlpt, typ, 2, NULL);
    addfunc("svpt", (rfunc)svpt, typ, 2, NULL);
    addfunc("slpt", (rfunc)slpt, typ, 2, NULL);
    addfunc("uvpt", (rfunc)uvpt, typ, 2, NULL);
    addfunc("ulpt", (rfunc)ulpt, typ, 2, NULL);
    addfunc("tau", (rfunc)tau, typ, 2, NULL);
    addfunc("taus", (rfunc)taus, typ, 2, NULL);
    addfunc("tauu", (rfunc)tauu, typ, 2, NULL);
    addfunc("vf", (rfunc)vf, typ, 2, NULL);
    addfunc("vfs", (rfunc)vfs, typ, 2, NULL);
    addfunc("vfu", (rfunc)vfu, typ, 2, NULL);
    addfunc("delta_liq", (rfunc)delta_liq, typ, 2, NULL);
    addfunc("delta_vap", (rfunc)delta_vap, typ, 2, NULL);
    // phi and derivatives for calculating more thermo properties.
    addfunc("phi0", (rfunc)phi0, typ, 2, NULL);
    addfunc("phir", (rfunc)phir, typ, 2, NULL);
    addfunc("phi0_d", (rfunc)phi0_d, typ, 2, NULL);
    addfunc("phir_d", (rfunc)phir_d, typ, 2, NULL);
    addfunc("phi0_t", (rfunc)phi0_t, typ, 2, NULL);
    addfunc("phir_t", (rfunc)phir_t, typ, 2, NULL);
    addfunc("phi0_dd", (rfunc)phi0_dd, typ, 2, NULL);
    addfunc("phir_dd", (rfunc)phir_dd, typ, 2, NULL);
    addfunc("phi0_dt", (rfunc)phi0_dt, typ, 2, NULL);
    addfunc("phir_dt", (rfunc)phir_dt, typ, 2, NULL);
    addfunc("phi0_tt", (rfunc)phi0_tt, typ, 2, NULL);
    addfunc("phir_tt", (rfunc)phir_tt, typ, 2, NULL);
    // Unary functions for sat curve
    addfunc("delta_sat_l", (rfunc)delta_sat_l, typ, 1, NULL);
    addfunc("delta_sat_v", (rfunc)delta_sat_v, typ, 1, NULL);
    addfunc("p_sat", (rfunc)p_sat, typ, 1, NULL);
    addfunc("tau_sat", (rfunc)tau_sat, typ, 1, NULL);
    // Parameters
    addfunc("mw", (rfunc)mw, typ, 0, NULL);
    addfunc("pc", (rfunc)pc, typ, 0, NULL);
    addfunc("tc", (rfunc)tc, typ, 0, NULL);
    addfunc("rhoc", (rfunc)rhoc, typ, 0, NULL);
    addfunc("pt", (rfunc)pt, typ, 0, NULL);
    addfunc("tt", (rfunc)tt, typ, 0, NULL);
    addfunc("rhot_v", (rfunc)rhot_l, typ, 0, NULL);
    addfunc("rhot_l", (rfunc)rhot_v, typ, 0, NULL);
    addfunc("sgc", (rfunc)sgc, typ, 0, NULL);
    addfunc("pmin", (rfunc)pmin, typ, 0, NULL);
    addfunc("tmin", (rfunc)tmin, typ, 0, NULL);
    addfunc("pmax", (rfunc)pmax, typ, 0, NULL);
    addfunc("tmax", (rfunc)tmax, typ, 0, NULL);
}

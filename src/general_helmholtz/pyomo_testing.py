import pyomo.environ as pyo
from helmholtz_functions import (
    add_helmholtz_external_functions,
    HelmholtzParameterBlock,
    HelmholtzThermoExpressions,
    AmountBasis,
)
from helmholtz_state import (
    HelmholtzStateBlock,
)


def main():
    m = pyo.ConcreteModel()
    m.param_block = HelmholtzParameterBlock(pure_component="h2o")

    m.r1234ze_ph_mass_param_block = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MASS
    )
    m.r1234ze_ph_mass_state_block = HelmholtzStateBlock(
        parameters=m.r1234ze_ph_mass_param_block
    )

    m.r1234ze_ph_mole_param_block = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MOLE
    )
    m.r1234ze_ph_mole_state_block = HelmholtzStateBlock(
        parameters=m.r1234ze_ph_mole_param_block
    )

    m.r1234ze_ph_mole_state_block.enth_mol = pyo.value(
        m.r1234ze_ph_mole_param_block.htpx(T=300 * pyo.units.K, p=101 * pyo.units.kPa)
    )
    m.r1234ze_ph_mole_state_block.pressure = 101
    m.r1234ze_ph_mole_state_block.entr_mass.display()

    add_helmholtz_external_functions(m)

    pc = pyo.value(m.pc_func("h2o"))
    tc = pyo.value(m.tc_func("h2o"))
    rhoc = pyo.value(m.rhoc_func("h2o"))
    pt = pyo.value(m.pc_func("h2o"))
    tt = pyo.value(m.tc_func("h2o"))
    rhot_l = pyo.value(m.rhot_l_func("h2o"))
    rhot_v = pyo.value(m.rhot_v_func("h2o"))
    sgc = pyo.value(m.sgc_func("h2o"))

    m.density = pyo.Var(initialize=961.618, doc="Density")
    m.temperature = pyo.Var(initialize=375, doc="Temperature")

    p = pyo.value(m.p_func("h2o", m.density / rhoc, tc / m.temperature))
    u = pyo.value(m.u_func("h2o", m.density / rhoc, tc / m.temperature))
    s = pyo.value(m.s_func("h2o", m.density / rhoc, tc / m.temperature))
    h = pyo.value(m.h_func("h2o", m.density / rhoc, tc / m.temperature))
    g = pyo.value(m.g_func("h2o", m.density / rhoc, tc / m.temperature))
    f = pyo.value(m.f_func("h2o", m.density / rhoc, tc / m.temperature))
    hvpt = pyo.value(m.hvpt_func("h2o", p, tc / m.temperature))
    hlpt = pyo.value(m.hlpt_func("h2o", p, tc / m.temperature))
    svpt = pyo.value(m.svpt_func("h2o", p, tc / m.temperature))
    slpt = pyo.value(m.slpt_func("h2o", p, tc / m.temperature))
    uvpt = pyo.value(m.uvpt_func("h2o", p, tc / m.temperature))
    ulpt = pyo.value(m.ulpt_func("h2o", p, tc / m.temperature))
    delta_vap = pyo.value(m.delta_vap_func("h2o", p, tc / m.temperature))
    delta_liq = pyo.value(m.delta_liq_func("h2o", p, tc / m.temperature))
    tauh = pyo.value(m.tau_func("h2o", h, p))
    vfh = pyo.value(m.vf_func("h2o", h, p))
    taus = pyo.value(m.taus_func("h2o", s, p))
    vfs = pyo.value(m.vfs_func("h2o", s, p))
    tauu = pyo.value(m.tauu_func("h2o", u, p))
    vfu = pyo.value(m.vfu_func("h2o", u, p))

    tau_sat = pyo.value(m.tau_sat_func("h2o", p))
    p_sat = pyo.value(m.p_sat_func("h2o", tc / m.temperature))
    delta_sat_l = pyo.value(m.delta_sat_l_func("h2o", tc / m.temperature))
    delta_sat_v = pyo.value(m.delta_sat_v_func("h2o", tc / m.temperature))

    phi0 = pyo.value(m.phi0_func("h2o", m.density / rhoc, tc / m.temperature))
    phi0_d = pyo.value(m.phi0_d_func("h2o", m.density / rhoc, tc / m.temperature))
    phi0_dd = pyo.value(m.phi0_dd_func("h2o", m.density / rhoc, tc / m.temperature))
    phi0_t = pyo.value(m.phi0_t_func("h2o", m.density / rhoc, tc / m.temperature))
    phi0_dt = pyo.value(m.phi0_dt_func("h2o", m.density / rhoc, tc / m.temperature))
    phi0_tt = pyo.value(m.phi0_tt_func("h2o", m.density / rhoc, tc / m.temperature))

    phir = pyo.value(m.phir_func("h2o", m.density / rhoc, tc / m.temperature))
    phir_d = pyo.value(m.phir_d_func("h2o", m.density / rhoc, tc / m.temperature))
    phir_dd = pyo.value(m.phir_dd_func("h2o", m.density / rhoc, tc / m.temperature))
    phir_t = pyo.value(m.phir_t_func("h2o", m.density / rhoc, tc / m.temperature))
    phir_dt = pyo.value(m.phir_dt_func("h2o", m.density / rhoc, tc / m.temperature))
    phir_tt = pyo.value(m.phir_tt_func("h2o", m.density / rhoc, tc / m.temperature))

    print("")
    print(f"Parameters {'h2o'}")
    print("-------------------------------------------------------------------")
    print(f"critical pressure = {pc} [kPa]")
    print(f"critical temperature = {tc} [K]")
    print(f"critical density = {rhoc} [kg/m^3]")
    print(f"tripple point pressure = {pt} [kPa]")
    print(f"tripple point temperature = {tt} [K]")
    print(f"tripple point liquid density = {rhot_l} [kg/m^3]")
    print(f"tripple point liquid density = {rhot_v} [kg/m^3]")
    print(f"specific gas constant = {sgc} [kJ/kg/K] or [kPa m^3/kg/K]")
    print("-------------------------------------------------------------------")
    print("")
    print(f"Basic Properties {'h2o'} ")
    print(f"  density = {pyo.value(m.density)} kg/m^3")
    print(f"  temperature = {pyo.value(m.temperature)} K")
    print("-------------------------------------------------------------------")
    print(f"p = {p} [kPa]")
    print(f"u = {u} [kJ/kg]")
    print(f"s = {s} [kJ/kg/K]")
    print(f"h = {h} [kJ/kg]")
    print(f"g = {g} [kJ/kg]")
    print(f"f = {f} [kJ/kg]")
    print("-------------------------------------------------------------------")
    print("")
    print(f"State Var Change PT {'h2o'}")
    print(f"  pressure = {p} kPa")
    print(f"  temperature = {pyo.value(m.temperature)} K")
    print("-------------------------------------------------------------------")
    print(f"hvpt = {hvpt} [kJ/kg]")
    print(f"hlpt = {hlpt} [kJ/kg]")
    print(f"svpt = {svpt} [kJ/kg]")
    print(f"slpt = {slpt} [kJ/kg]")
    print(f"uvpt = {uvpt} [kJ/kg]")
    print(f"ulpt = {ulpt} [kJ/kg]")
    print(f"rho_vap = {delta_vap*rhoc} [kg/m^3]")
    print(f"rho_liq = {delta_liq*rhoc} [kg/m^3]")
    print("-------------------------------------------------------------------")
    print("")
    print(f"State Var Change HP {'h2o'}")
    print(f"  enthalpy = {h} kJ/kg")
    print(f"  pressure = {p} kPa")
    print("-------------------------------------------------------------------")
    print(f"T = {tc/tauh} [K]")
    print(f"vf = {vfh} [none]")
    print("-------------------------------------------------------------------")
    print("")
    print(f"State Var Change SP {'h2o'}")
    print(f"  entropy = {s} kJ/kg/K")
    print(f"  pressure = {p} kPa")
    print("-------------------------------------------------------------------")
    print(f"T = {tc/taus} [K]")
    print(f"vf = {vfs} [none]")
    print("-------------------------------------------------------------------")
    print("")
    print(f"State Var Change UP {'h2o'}")
    print(f"  internal energy = {u} kJ/kg")
    print(f"  pressure = {p} kPa")
    print("-------------------------------------------------------------------")
    print(f"T = {tc/tauu} [K]")
    print(f"vf = {vfu} [none]")
    print("-------------------------------------------------------------------")
    print("")
    print(f"Saturation Curve {'h2o'}")
    print(f"  temprature = {pyo.value(m.temperature)} K")
    print("-------------------------------------------------------------------")
    print(f"p_sat = {p_sat}")
    print(f"rho_sat_l = {delta_sat_l*rhoc}")
    print(f"rho_sat_v = {delta_sat_v*rhoc}")
    print("-------------------------------------------------------------------")
    print("")
    print(f"Saturation Curve {'h2o'}")
    print(f"  pressure = {p} kPa")
    print("-------------------------------------------------------------------")
    print(f"T_sat = {tc/tau_sat} K")
    print("-------------------------------------------------------------------")
    print("")
    print(f"Dimensionless Helmholtz Free Energy {'h2o'}")
    print(f"  density = {pyo.value(m.density)} kg/m^3")
    print(f"  temperature = {pyo.value(m.temperature)} K")
    print("Can be used to calculate more thermo properties (e.g. cv, cp, ...)")
    print("-------------------------------------------------------------------")
    print(f"phi0 = {phi0} [none]")
    print(f"phi0_d = {phi0_d} [none]")
    print(f"phi0_t = {phi0_t} [none]")
    print(f"phi0_dd = {phi0_dd} [none]")
    print(f"phi0_dt = {phi0_dt} [none]")
    print(f"phi0_tt = {phi0_tt} [none]")
    print(f"phir = {phir} [none]")
    print(f"phir_d = {phir_d} [none]")
    print(f"phir_t = {phir_t} [none]")
    print(f"phir_dd = {phir_dd} [none]")
    print(f"phir_dt = {phir_dt} [none]")
    print(f"phir_tt = {phir_tt} [none]")
    print("-------------------------------------------------------------------")

    m.param_block.ph_diagram()
    m.param_block.st_diagram()
    m.param_block.pt_diagram()

    m.param_block2 = HelmholtzParameterBlock(pure_component="r1234ze")
    h = m.param_block2.htpx(
        T=200 * pyo.units.K,
        p=101.325 * pyo.units.kPa,
        amount_basis=AmountBasis.MASS,
        with_units=True,
    )
    print(pyo.value(pyo.units.convert(h, pyo.units.kJ / pyo.units.kg)))

    m.param_block2.temperature_crit.display()
    m.param_block2.ph_diagram()
    m.param_block2.st_diagram()
    m.param_block2.pt_diagram()

    m.param_block3 = HelmholtzParameterBlock(pure_component="co2")
    m.param_block3.temperature_crit.display()
    m.param_block3.ph_diagram()
    m.param_block3.st_diagram()
    m.param_block3.pt_diagram()

if __name__ == "__main__":
    main()

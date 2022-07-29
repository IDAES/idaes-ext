import pyomo.environ as pyo
from pyomo.common.fileutils import find_library

_flib = find_library("general_helmholtz_external")


def available():
    """Return whether the shared library is available. If it is not this cannot
    be used."""
    return _flib is not None


_external_function_map = {
    # Thermo properties as a function of delta and tau
    "p_func": {  # pressure
        "fname": "p",
        "units": pyo.units.kPa,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "u_func": {  # internal energy
        "fname": "u",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "s_func": {  # entropy
        "fname": "s",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "h_func": {  # enthaply
        "fname": "h",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "g_func": {  # Gibbs free energy
        "fname": "g",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "f_func": {  # Helmholtz free energy
        "fname": "f",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    # Dimensionless Helmholtz energy to calculate other thermo properties
    "phi0_func": {  # ideal part
        "fname": "phi0",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phi0_d_func": {  # ideal part derivative wrt delta
        "fname": "phi0_d",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phi0_dd_func": {  # ideal part second derivative wrt delta
        "fname": "phi0_dd",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phi0_t_func": {  # ideal part derivative wrt tau
        "fname": "phi0_t",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phi0_dt_func": {  # ideal part second derivative wrt delta and tau
        "fname": "phi0_dt",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phi0_tt_func": {  # ideal part second derivative wrt tau
        "fname": "phi0_tt",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phir_func": {  # residual part
        "fname": "phir",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phir_d_func": {  # residual part derivative wrt delta
        "fname": "phir_d",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phir_dd_func": {  # residual part second derivative wrt delta
        "fname": "phir_dd",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phir_t_func": {  # residual part derivative wrt tau
        "fname": "phir_t",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phir_dt_func": {  # residual part second derivative wrt delta and tau
        "fname": "phir_dt",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    "phir_tt_func": {  # residual part second derivative wrt tau
        "fname": "phir_tt",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless, pyo.units.dimensionless],
    },
    # Phase specific functions of pressure and tau
    "hvpt_func": {  # vapor enthalpy
        "fname": "hvpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "hlpt_func": {  # liquid enthalpy
        "fname": "hlpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "svpt_func": {  # vapor entropy
        "fname": "svpt",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "slpt_func": {  # liquid entropy
        "fname": "slpt",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "uvpt_func": {  # vapor internal energy
        "fname": "uvpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "ulpt_func": {  # liquid internal energy
        "fname": "ulpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "delta_liq_func": {  # liquid density
        "fname": "delta_liq",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    "delta_vap_func": {  # vapor density
        "fname": "delta_vap",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kPa, pyo.units.dimensionless],
    },
    # state variable change functions
    "tau_func": {  # tau as a function of h, p
        "fname": "tau",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
    },
    "vf_func": {  # vapor fraction as a function of h, p
        "fname": "vf",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
    },
    "taus_func": {  # tau as a function of s, p
        "fname": "taus",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    "vfs_func": {  # vapor fraction as a function of s, p
        "fname": "vfs",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    "tauu_func": {  # tau as a function of u, p
        "fname": "tauu",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    "vfu_func": {  # vapor fraction as a function of u, p
        "fname": "vfu",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    # saturation curve as a function of tau
    "p_sat_func": {
        "fname": "p_sat",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless],
    },
    "delta_sat_v_func": {
        "fname": "delta_sat_v",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless],
    },
    "delta_sat_l_func": {
        "fname": "delta_sat_l",
        "units": pyo.units.dimensionless,
        "arg_units": [pyo.units.dimensionless],
    },
    # saturation curve as a function of p
    "tau_sat_func": {
        "fname": "tau_sat",
        "units": pyo.units.kPa,
        "arg_units": [pyo.units.dimensionless],
    },
    # Parameters (these functions take no arguments)
    "sgc_func": {
        "fname": "sgc",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [],
    },
    "pc_func": {
        "fname": "pc",
        "units": pyo.units.kPa,
        "arg_units": [],
    },
    "tc_func": {
        "fname": "tc",
        "units": pyo.units.K,
        "arg_units": [],
    },
    "rhoc_func": {
        "fname": "rhoc",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [],
    },
    "pt_func": {
        "fname": "pt",
        "units": pyo.units.kPa,
        "arg_units": [],
    },
    "tt_func": {
        "fname": "tt",
        "units": pyo.units.K,
        "arg_units": [],
    },
    "rhot_l_func": {
        "fname": "rhot_l",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [],
    },
    "rhot_v_func": {
        "fname": "rhot_v",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [],
    },
    "pmin_func": {
        "fname": "pmin",
        "units": pyo.units.kPa,
        "arg_units": [],
    },
    "tmin_func": {
        "fname": "tmin",
        "units": pyo.units.K,
        "arg_units": [],
    },
    "pmax_func": {
        "fname": "pmin",
        "units": pyo.units.kPa,
        "arg_units": [],
    },
    "tmax_func": {
        "fname": "tmin",
        "units": pyo.units.K,
        "arg_units": [],
    },
}


def add_helmholtz_external_functions(blk, names=None):
    """Add a Helmholtz EoS function to a Pyomo Block.

    Args:
        blk: block to add function to
        names: if None, add all functions, if a string add the single function,
            otherwise a list of function names.

    Returns:
        None
    """
    if names is None:
        names = _external_function_map.keys()
    if isinstance(names, str):
        names = [names]
    for name in names:
        if hasattr(blk, name):
            continue
        fdict = _external_function_map[name]
        setattr(
            blk,
            name,
            pyo.ExternalFunction(
                library=_flib,
                function=fdict["fname"],
                units=fdict["units"],
                arg_units=fdict["arg_units"],
            ),
        )

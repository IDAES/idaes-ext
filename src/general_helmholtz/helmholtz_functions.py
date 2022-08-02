from matplotlib import pyplot as plt
import numpy as np

import pyomo.environ as pyo
from pyomo.common.fileutils import find_library
from pyomo.common.config import ConfigValue, In
from idaes.core import declare_process_block_class
from idaes.core import (
    StateBlock,
    StateBlockData,
    PhysicalParameterBlock,
)

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
    "mw_func": {
        "fname": "mw",
        "units": pyo.units.g / pyo.units.mol,
        "arg_units": [],
    },
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


class HelmholtzThermoExpressions(object):
    """Class to write thermodynamic property expressions.  Take one of these
    possible sets of state variables: {h, p}, {u, p}, {s, p}, {s, T}, {T, x},
    {P, x}, or {T, P, x}, and return an expression for a thermo property.
    This works by converting the given state varaibles to temperature, density,
    and vapor fraction expressions then using those to write an expression for
    requested property. This writes expressions in a way that looks like a
    thermodynaic property function.
    """

    def __init__(self, blk, parameters):
        """Create a new thermodynamic property expression writer class.

        Args:
            blk: the block to attach the external functions to
            parameters: property parameter block

        Returns:
            HelmholtzThermoExpressions
        """
        self.param = parameters
        self.blk = blk

    @staticmethod
    def _sv_str(**kwargs):
        a = [x for x in kwargs if kwargs[x] is not None]
        return ", ".join(a)

    def add_funcs(self, names=None):
        add_helmholtz_external_functions(self.blk, names=names)

    def basic_calculations(self, h=None, s=None, p=None, T=None, u=None, x=None):
        """This function is called as the basis for most thermo expression
        writer functions.  It takes the given state variables and returns
        expressions for liquid density, vapor density, vapor fraction and
        temperature, which can be used to write an expression for any thermo
        quantity.
        """
        mw = self.param.mw
        c = self.param.pure_component
        # 1.) convert units to those expected by external functions
        if h is not None:
            h *= self.param.uc_J_per_mol_to_kJ_per_kg
        if u is not None:
            u *= self.param.uc_J_per_mol_to_kJ_per_kg
        if s is not None:
            s *= self.param.uc_J_per_molK_to_kJ_per_kgK
        if p is not None:
            p *= self.param.uc_Pa_to_kPa
        if T is not None:
            tau = self.param.temperature_crit / T
        # 2.) find the block with the external functions
        blk = self.blk

        # 3.) Take given state varaibles and convert to density, T, and x
        if h is not None and p is not None:
            # h, p
            self.add_funcs(names=["tau_func", "vf_func"])
            tau = blk.tau_func(c, h, p)
            if x is None:
                x = blk.vf_func(c, h, p)
        elif s is not None and p is not None:
            # s, p
            self.add_funcs(names=["taus_func", "vfs_func"])
            tau = blk.taus_func(c, s, p)
            if x is None:
                x = blk.vfs_func(c, s, p)
        elif u is not None and p is not None:
            # u, p
            self.add_funcs(names=["tauu_func", "vfu_func"])
            tau = blk.tauu_func(c, u, p)
            if x is None:
                x = blk.vfu_func(c, u, p)
        elif x is not None and T is not None and p is not None:
            # T, P, x (okay, but I hope you know what you're doing)
            pass
        elif x is not None and p is not None:
            # x, p
            self.add_funcs(names=["tau_sat_func"])
            tau = blk.tau_sat_func(c, p)
        elif x is not None and T is not None:
            # x, T
            self.add_funcs(names=["p_sat_func"])
            p = blk.p_sat_func(c, tau)
        else:
            m = "This choice of state variables ({}) is not yet supported.".format(
                self._sv_str(h=h, s=s, p=p, T=T, u=u, x=x)
            )
            _log.error(m)
            raise NotImplementedError(m)

        # 4.) Calculate density
        self.add_funcs(names=["delta_liq_func", "delta_vap_func"])
        delta_liq = blk.delta_liq_func(c, p, tau)
        delta_vap = blk.delta_vap_func(c, p, tau)

        # 5.) From here its straight forward to calculate any property
        return blk, delta_liq, delta_vap, tau, x, c

    def s(self, **kwargs):
        """Mixed phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_liq, tau) * (1 - x) + blk.s_func(c, delta_vap, tau) * x
        return s * self.param.uc_kJ_per_kgK_to_J_per_molK

    def s_liq(self, **kwargs):
        """Liquid phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_liq, tau)
        return s * self.param.uc_kJ_per_kgK_to_J_per_molK

    def s_vap(self, **kwargs):
        """Vapor phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_vap, tau)
        return s * self.param.uc_kJ_per_kgK_to_J_per_molK

    def h(self, **kwargs):
        """Mixed phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_liq, tau) * (1 - x) + blk.h_func(c, delta_vap, tau) * x
        return h * self.param.uc_kJ_per_kg_to_J_per_mol

    def h_liq(self, **kwargs):
        """Liquid phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_liq, tau)
        return h * self.param.uc_kJ_per_kg_to_J_per_mol

    def h_vap(self, **kwargs):
        """Vapor phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_vap, tau)
        return h * self.param.uc_kJ_per_kg_to_J_per_mol

    def u(self, **kwargs):
        """Mixed phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_liq, tau) * (1 - x) + blk.u_func(c, delta_vap, tau) * x
        return u * self.param.uc_kJ_per_kg_to_J_per_mol

    def u_liq(self, **kwargs):
        """Liquid phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_liq, tau)
        return u * self.param.uc_kJ_per_kg_to_J_per_mol

    def u_vap(self, **kwargs):
        """Vapor phase internal energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_vap, tau)
        return u * self.param.uc_kJ_per_kg_to_J_per_mol

    def g(self, **kwargs):
        """Mixed phase Gibb's free energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_liq, tau) * (1 - x) + blk.g_func(c, delta_vap, tau) * x
        return g * self.param.uc_kJ_per_kg_to_J_per_mol

    def g_liq(self, **kwargs):
        """Liquid phase Gibb's free energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_liq, tau)
        return g * self.param.uc_kJ_per_kg_to_J_per_mol

    def g_vap(self, **kwargs):
        """Vapor phase Gibb's free energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_vap, tau)
        return g * self.param.uc_kJ_per_kg_to_J_per_mol

    def f(self, **kwargs):
        """Mixed phase Helmholtz free energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_liq, tau) * (1 - x) + blk.f_func(c, delta_vap, tau) * x
        return f * self.param.uc_kJ_per_kg_to_J_per_mol

    def f_liq(self, **kwargs):
        """Liquid phase Helmholtz free energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_liq, tau)
        return f * self.param.uc_kJ_per_kg_to_J_per_mol

    def f_vap(self, **kwargs):
        """Vapor phase Helmholtz free energy"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_vap, tau)
        return f * self.param.uc_kJ_per_kg_to_J_per_mol

    def p(self, **kwargs):
        """Pressure"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["p_func"])
        # The following line looks a bit weird, but it is okay.  When in the
        # two-phase region the pressure for both phases is the same
        p = blk.p_func(c, delta_liq, tau) * (1 - x) + blk.p_func(c, delta_vap, tau) * x
        return p * self.param.uc_kPa_to_Pa

    def v(self, **kwargs):
        """Mixed phase molar volume"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = (
            ((1 - x) / delta_liq + x / delta_vap)
            / self.param.dens_mass_crit
            * self.param.mw
        )
        return v

    def v_liq(self, **kwargs):
        """Liquid phase molar volume"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = 1.0 / delta_liq / self.param.dens_mass_crit * self.param.mw
        return v

    def v_vap(self, **kwargs):
        """Vapor phase molar volume"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = 1.0 / delta_vap / self.param.dens_mass_crit * self.param.mw
        return v

    def x(self, **kwargs):
        """Vapor faction"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return x

    def T(self, **kwargs):
        """Temperature"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return self.param.temperature_crit / tau

    def tau(self, **kwargs):
        """Critical Temperature (K)/Temperature (K)"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return tau

    def delta_liq(self, **kwargs):
        """Return liquid phase reduced density (dens/critical dens) expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq

    def rho_liq(self, **kwargs):
        """Return liquid phase mass density expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq * self.param.dens_mass_crit

    def rho_mol_liq(self, **kwargs):
        """Return liquid phase molar density expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq * self.param.dens_mass_crit / self.param.mw

    def delta_vap(self, **kwargs):
        """Return vapor phase reduced density (dens/critical dens) expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap

    def rho_vap(self, **kwargs):
        """Return vapor phase mass density expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap * self.param.dens_mass_crit

    def rho_mol_vap(self, **kwargs):
        """Return vapor phase molar density expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap * self.param.dens_mass_crit / self.param.mw

    def cv_mol_liq(self, **kwargs):
        """Return liquid phase constant volume heat capacity expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        return blk.func_cv(delta_liq, tau) * self.param.uc_kJ_per_kgK_to_J_per_molK

    def cv_mol_vap(self, **kwargs):
        """Return vapor phase constant volume heat capacity expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_cv"])
        return blk.func_cv(delta_vap, tau) * self.param.uc_kJ_per_kgK_to_J_per_molK

    def w_liq(self, **kwargs):
        """Return liquid phase speed of sound expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_w"])
        return blk.func_w(delta_liq, tau)

    def w_vap(self, **kwargs):
        """Return vapor phase speed of sound expression"""
        pyo.unblk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["func_w"])
        return blk.func_w(delta_vap, tau)

    def p_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_crit / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("p_sat expression requires either T or tau arg")
        return self.blk.p_sat_func(self.param.pure_component, tau) * self.param.uc_kPa_to_Pa

    def T_sat(self, p):
        """Return saturation temperature as a function of p"""
        p *= self.param.uc_Pa_to_kPa
        return self.param.temperature_crit / self.blk.tau_sat_func(self.param.pure_component, p)

    def tau_sat(self, p):
        """Return saturation tau as a function of p"""
        p *= self.param.uc_Pa_to_kPa
        return self.blk.func_tau_sat(self.param.pure_component, p)


@declare_process_block_class("HelmholtzStateBlock")
class HelmholtzStateBlockData(StateBlockData):
    """
    This is a base clase for Helmholtz equations of state using IDAES standard
    Helmholtz EOS external functions written in C++.
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "pure_component",
        ConfigValue(
            default=None,
            domain=str,
            description="Pure chemical component",
            doc="Pure component to calculate properies for",
        ),
    )

    def build(self):
        # set the component_list as required for the generic IDAES properties
        self.componet_list = pyo.Set(initialize=[self.config.pure_component])
        # sinice this a only a pure component package, have a specific
        # pure_component attirbute
        self.pure_component = self.config.pure_component
        # To ensure consistency pull parameters from external functions
        add_helmholtz_external_functions(
            self,
            [
                # Constants
                "sgc_func",  # specific gas constant
                "mw_func",  # molecular weight
                # Critical properties
                "pc_func",  # Critical pressure
                "tc_func",  # Critical temperature
                "rhoc_func",  # Critical density
                # Tripple point properties
                "pt_func",  # Tripple point pressure
                "tt_func",  # Tripple point temperature
                "rhot_l_func",  # Tripple point liquid density
                "rhot_v_func",  # Tripple point vapor density
                # Tripple point properties
                "pt_func",  # Tripple point pressure
                "tt_func",  # Tripple point temperature
                "rhot_l_func",  # Tripple point liquid density
                "rhot_v_func",  # Tripple point vapor density
                # Bounds
                "pmin_func",  # pmin
                "tmin_func",  # tmin
                "pmax_func",  # pmax
                "tmax_func",  # tmax
            ],
        )
        # The parameters are constants and we don't want to call the external
        # functions more than once, so define Pyomo parameters
        self.sgc = pyo.Param(
            initialize=self.sgc_func(self.pure_component) * 1000,
            units=pyo.units.J / pyo.units.kg / pyo.units.K,
        )
        self.mw = pyo.Param(
            initialize=self.mw_func(self.pure_component) / 1000,
            units=pyo.units.kg / pyo.units.mol,
        )
        self.pressure_crit = pyo.Param(
            initialize=self.pc_func(self.pure_component) * 1000, units=pyo.units.Pa
        )
        self.temperature_crit = pyo.Param(
            initialize=self.tc_func(self.pure_component), units=pyo.units.K
        )
        self.pressure_trip = pyo.Param(
            initialize=self.pt_func(self.pure_component) * 1000, units=pyo.units.Pa
        )
        self.temperature_trip = pyo.Param(
            initialize=self.tt_func(self.pure_component), units=pyo.units.K
        )
        self.dens_mass_crit = pyo.Param(
            initialize=self.rhoc_func(self.pure_component),
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.uc_J_per_mol_to_kJ_per_kg = pyo.units.kJ / self.mw / 1000 / pyo.units.J
        self.uc_kJ_per_kg_to_J_per_mol = self.mw * pyo.units.J * 1000 / pyo.units.kJ
        # Entropy type units are same as energy type units since K is left as is
        self.uc_J_per_molK_to_kJ_per_kgK = self.uc_J_per_mol_to_kJ_per_kg
        self.uc_kJ_per_kgK_to_J_per_molK = self.uc_kJ_per_kg_to_J_per_mol
        # Pressure conversions
        self.uc_kPa_to_Pa = 1000 * pyo.units.Pa / pyo.units.kPa
        self.uc_Pa_to_kPa = pyo.units.kPa / pyo.units.Pa / 1000

    def initialize(self, *args, **kwargs):
        pass

    def ph_diagram(self):
        add_helmholtz_external_functions(
            self,
            [
                "p_sat_func",
                "tau_sat_func",
                "delta_sat_l_func",
                "delta_sat_v_func",
                "h_func",
                "hlpt_func",
                "hvpt_func",
                "delta_liq_func",
                "delta_vap_func",
            ]
        )

        tau_sat_vec = np.linspace(1, pyo.value(self.temperature_crit)/pyo.value(self.temperature_trip), 200)
        p_sat_vec = [None]*len(tau_sat_vec)
        delta_sat_v_vec = [None]*len(tau_sat_vec)
        delta_sat_l_vec = [None]*len(tau_sat_vec)
        h_sat_v_vec = [None]*len(tau_sat_vec)
        h_sat_l_vec = [None]*len(tau_sat_vec)

        for i, tau in enumerate(tau_sat_vec):
            p_sat_vec[i] = pyo.value(self.p_sat_func(self.pure_component, tau))
            delta_sat_l_vec[i] = pyo.value(self.delta_sat_l_func(self.pure_component, tau))
            delta_sat_v_vec[i] = pyo.value(self.delta_sat_v_func(self.pure_component, tau))
            h_sat_v_vec[i] = pyo.value(self.h_func(self.pure_component, delta_sat_v_vec[i], tau))
            h_sat_l_vec[i] = pyo.value(self.h_func(self.pure_component, delta_sat_l_vec[i], tau))

        plt.yscale("log")
        plt.plot(h_sat_l_vec, p_sat_vec, c="b", label="sat liquid")
        plt.plot(h_sat_v_vec, p_sat_vec, c="r", label="sat vapor")

        t_vec = [pyo.value(self.temperature_trip), 300, 400, 500, 600, pyo.value(self.temperature_crit)]
        #t_vec = [pyo.value(self.temperature_trip), 300, 400, 500, 600]
        p = {}
        h_l = {}
        h_v = {}

        # plot isotherms in sat region
        for t in t_vec:
            tau = pyo.value(self.temperature_crit)/t
            p[t] = pyo.value(self.p_sat_func(self.pure_component, tau))
            delta_l = pyo.value(self.delta_sat_l_func(self.pure_component, tau))
            delta_v = pyo.value(self.delta_sat_v_func(self.pure_component, tau))
            h_v[t] = pyo.value(self.h_func(self.pure_component, delta_v, tau))
            h_l[t] = pyo.value(self.h_func(self.pure_component, delta_l, tau))
            plt.plot([h_l[t], h_v[t]], [p[t], p[t]], c='g')
            plt.text(h_l[t]/2 + h_v[t]/2, p[t], f"T = {t} K", ha='center')

        for t in t_vec:
            tau = pyo.value(self.temperature_crit)/t
            p_vec = np.linspace(p[t], 1e6, 100)
            h_vec = [None]*len(p_vec)
            for i, pv in enumerate(p_vec):
                h_vec[i] = pyo.value(self.hlpt_func(self.pure_component, pv, tau))
            plt.plot(h_vec, p_vec, c="g")

        for t in t_vec:
            tau = pyo.value(self.temperature_crit)/t
            p_vec = np.linspace(0.1, p[t], 50)
            h_vec = [None]*len(p_vec)
            for i, pv in enumerate(p_vec):
                h_vec[i] = pyo.value(self.hvpt_func(self.pure_component, pv, tau))
            plt.plot(h_vec, p_vec, c="g")

        plt.xlabel("Enthalpy (kJ/kg)")
        plt.ylabel("Pressure (kPa)")






        plt.show()

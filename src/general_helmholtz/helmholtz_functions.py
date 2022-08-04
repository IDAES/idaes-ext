import enum

from matplotlib import pyplot as plt
import numpy as np

import pyomo.environ as pyo
from pyomo.core.base.units_container import InconsistentUnitsError
from pyomo.common.fileutils import find_library
from pyomo.common.config import ConfigValue, In

from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from idaes.core import (
    StateBlock,
    StateBlockData,
    PhysicalParameterBlock,
    LiquidPhase,
    VaporPhase,
    Phase,
    Component,
)

_flib = find_library("general_helmholtz_external.so")

def available():
    """Return whether the shared library is available. If it is not this cannot
    be used."""
    return _flib is not None


class StateVars(enum.Enum):
    """
    State variable set options
    """

    PH = 1  # Pressure-Enthalpy
    TPX = 2  # Temperature-Pressure-Quality


class PhaseType(enum.Enum):
    """
    Ways to present phases to the framework
    """

    MIX = 1  # Looks like a single phase called mixed with a vapor fraction
    LG = 2  # Looks like two phases vapor and liquid
    L = 3  # Assume only liquid is present
    G = 4  # Assume only vapor is pressent


class AmountBasis(enum.Enum):
    """
    Mass or mole basis
    """
    MOLE = 1
    MASS = 2


dimensionless = pyo.units.dimensionless


_external_function_map = {
    # Thermo properties as a function of delta and tau
    "p_func": {  # pressure
        "fname": "p",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "u_func": {  # internal energy
        "fname": "u",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "s_func": {  # entropy
        "fname": "s",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "h_func": {  # enthaply
        "fname": "h",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "g_func": {  # Gibbs free energy
        "fname": "g",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "f_func": {  # Helmholtz free energy
        "fname": "f",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "cv_func": {  # entropy
        "fname": "cv",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "cp_func": {  # entropy
        "fname": "cp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    # Dimensionless Helmholtz energy to calculate other thermo properties
    "phi0_func": {  # ideal part
        "fname": "phi0",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_d_func": {  # ideal part derivative wrt delta
        "fname": "phi0_d",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_dd_func": {  # ideal part second derivative wrt delta
        "fname": "phi0_dd",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_t_func": {  # ideal part derivative wrt tau
        "fname": "phi0_t",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_dt_func": {  # ideal part second derivative wrt delta and tau
        "fname": "phi0_dt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_tt_func": {  # ideal part second derivative wrt tau
        "fname": "phi0_tt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_func": {  # residual part
        "fname": "phir",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_d_func": {  # residual part derivative wrt delta
        "fname": "phir_d",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_dd_func": {  # residual part second derivative wrt delta
        "fname": "phir_dd",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_t_func": {  # residual part derivative wrt tau
        "fname": "phir_t",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_dt_func": {  # residual part second derivative wrt delta and tau
        "fname": "phir_dt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_tt_func": {  # residual part second derivative wrt tau
        "fname": "phir_tt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    # Phase specific functions of pressure and tau
    "hvpt_func": {  # vapor enthalpy
        "fname": "hvpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "hlpt_func": {  # liquid enthalpy
        "fname": "hlpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "svpt_func": {  # vapor entropy
        "fname": "svpt",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "slpt_func": {  # liquid entropy
        "fname": "slpt",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "uvpt_func": {  # vapor internal energy
        "fname": "uvpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "ulpt_func": {  # liquid internal energy
        "fname": "ulpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "delta_liq_func": {  # liquid density
        "fname": "delta_liq",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "delta_vap_func": {  # vapor density
        "fname": "delta_vap",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    # state variable change functions
    "tau_func": {  # tau as a function of h, p
        "fname": "tau",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
    },
    "vf_func": {  # vapor fraction as a function of h, p
        "fname": "vf",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
    },
    "taus_func": {  # tau as a function of s, p
        "fname": "taus",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    "vfs_func": {  # vapor fraction as a function of s, p
        "fname": "vfs",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    "tauu_func": {  # tau as a function of u, p
        "fname": "tauu",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    "vfu_func": {  # vapor fraction as a function of u, p
        "fname": "vfu",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg / pyo.units.K, pyo.units.kPa],
    },
    # saturation curve as a function of tau
    "p_sat_func": {
        "fname": "p_sat",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless, dimensionless],
    },
    "delta_sat_v_func": {
        "fname": "delta_sat_v",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless],
    },
    "delta_sat_l_func": {
        "fname": "delta_sat_l",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless],
    },
    # saturation curve as a function of p
    "tau_sat_func": {
        "fname": "tau_sat",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kPa],
    },
    # Parameters (these functions take no arguments)
    "mw_func": {
        "fname": "mw",
        "units": pyo.units.g / pyo.units.mol,
        "arg_units": [dimensionless],
    },
    "sgc_func": {
        "fname": "sgc",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless],
    },
    "pc_func": {
        "fname": "pc",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tc_func": {
        "fname": "tc",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "rhoc_func": {
        "fname": "rhoc",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "pt_func": {
        "fname": "pt",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tt_func": {
        "fname": "tt",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "rhot_l_func": {
        "fname": "rhot_l",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "rhot_v_func": {
        "fname": "rhot_v",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "pmin_func": {
        "fname": "pmin",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tmin_func": {
        "fname": "tmin",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "pmax_func": {
        "fname": "pmax",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tmax_func": {
        "fname": "tmax",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
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

    def __init__(self, blk, parameters, amount_basis=AmountBasis.MOLE):
        """Create a new thermodynamic property expression writer class.

        Args:
            blk: the block to attach the external functions to
            parameters: property parameter block

        Returns:
            HelmholtzThermoExpressions
        """
        self.param = parameters
        self.blk = blk
        self.amount_basis = amount_basis

    @staticmethod
    def _sv_str(**kwargs):
        a = [x for x in kwargs if kwargs[x] is not None]
        return ", ".join(a)

    def add_funcs(self, names=None):
        add_helmholtz_external_functions(self.blk, names=names)

    def basic_calculations(
        self, h=None, s=None, p=None, T=None, u=None, x=None, tau=None
    ):
        """This function is called as the basis for most thermo expression
        writer functions.  It takes the given state variables and returns
        expressions for liquid density, vapor density, vapor fraction and
        temperature, which can be used to write an expression for any thermo
        quantity.
        """
        mw = self.param.mw
        c = self.param.pure_component

        # 1.) convert units to those expected by external functions
        if self.amount_basis==AmountBasis.MOLE:
            if h is not None:
                h *= self.param.uc["J/mol to kJ/kg"]
            if u is not None:
                u *= self.param.uc["J/mol to kJ/kg"]
            if s is not None:
                s *= self.param.uc["J/mol/K to kJ/kg/K"]
        else:
            if h is not None:
                h *= self.param.uc["J/kg to kJ/kg"]
            if u is not None:
                u *= self.param.uc["J/kg to kJ/kg"]
            if s is not None:
                s *= self.param.uc["J/kg/K to kJ/kg/K"]
        if p is not None:
            p *= self.param.uc["Pa to kPa"]
        if T is not None:
            tau = self.param.temperature_crit / T
        if tau is not None:
            T = self.param.temperature_crit / tau
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
        if self.amount_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def s_liq(self, **kwargs):
        """Liquid phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def s_vap(self, **kwargs):
        """Vapor phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def h(self, **kwargs):
        """Mixed phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_liq, tau) * (1 - x) + blk.h_func(c, delta_vap, tau) * x
        if self.amount_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg/K to J/mol/K"]
        return h * self.param.uc["kJ/kg/K to J/kg/K"]

    def h_liq(self, **kwargs):
        """Liquid phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg/K to J/mol/K"]
        return h * self.param.uc["kJ/kg/K to J/kg/K"]

    def h_vap(self, **kwargs):
        """Vapor phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg/K to J/mol/K"]
        return h * self.param.uc["kJ/kg/K to J/kg/K"]

    def u(self, **kwargs):
        """Mixed phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_liq, tau) * (1 - x) + blk.u_func(c, delta_vap, tau) * x
        if self.amount_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg/K to J/mol/K"]
        return u * self.param.uc["kJ/kg/K to J/kg/K"]

    def u_liq(self, **kwargs):
        """Liquid phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg/K to J/mol/K"]
        return u * self.param.uc["kJ/kg/K to J/kg/K"]

    def u_vap(self, **kwargs):
        """Vapor phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg/K to J/mol/K"]
        return u * self.param.uc["kJ/kg/K to J/kg/K"]

    def g(self, **kwargs):
        """Mixed phase Gibb's free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_liq, tau) * (1 - x) + blk.g_func(c, delta_vap, tau) * x
        if self.amount_basis == AmountBasis.MOLE:
            return g * self.param.uc["kJ/kg/K to J/mol/K"]
        return g * self.param.uc["kJ/kg/K to J/kg/K"]

    def g_liq(self, **kwargs):
        """Liquid phase Gibb's free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return g * self.param.uc["kJ/kg/K to J/mol/K"]
        return g * self.param.uc["kJ/kg/K to J/kg/K"]

    def g_vap(self, **kwargs):
        """Vapor phase Gibb's free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return g * self.param.uc["kJ/kg/K to J/mol/K"]
        return g * self.param.uc["kJ/kg/K to J/kg/K"]

    def f(self, **kwargs):
        """Mixed phase Helmholtz free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_liq, tau) * (1 - x) + blk.f_func(c, delta_vap, tau) * x
        if self.amount_basis == AmountBasis.MOLE:
            return f * self.param.uc["kJ/kg/K to J/mol/K"]
        return f * self.param.uc["kJ/kg/K to J/kg/K"]

    def f_liq(self, **kwargs):
        """Liquid phase Helmholtz free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return f * self.param.uc["kJ/kg/K to J/mol/K"]
        return f * self.param.uc["kJ/kg/K to J/kg/K"]

    def f_vap(self, **kwargs):
        """Vapor phase Helmholtz free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return f * self.param.uc["kJ/kg/K to J/mol/K"]
        return f * self.param.uc["kJ/kg/K to J/kg/K"]

    def p(self, **kwargs):
        """Pressure"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["p_func"])
        # The following line looks a bit weird, but it is okay.  When in the
        # two-phase region the pressure for both phases is the same
        p = blk.p_func(c, delta_liq, tau) * (1 - x) + blk.p_func(c, delta_vap, tau) * x
        return p * self.param.uc_kPa_to_Pa

    def v_mol(self, **kwargs):
        """Mixed phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = (
            ((1 - x) / delta_liq + x / delta_vap)
            / self.param.dens_mass_crit
            * self.param.mw
        )
        return v

    def v_mol_liq(self, **kwargs):
        """Liquid phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = self.param.mw / delta_liq / self.param.dens_mass_crit
        return v

    def v_mol_vap(self, **kwargs):
        """Vapor phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = self.param.mw / delta_vap / self.param.dens_mass_crit
        return v

    def x(self, **kwargs):
        """Vapor faction"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return x

    def T(self, **kwargs):
        """Temperature"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return self.param.temperature_crit / tau

    def tau(self, **kwargs):
        """Critical Temperature (K)/Temperature (K)"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return tau

    def delta_liq(self, **kwargs):
        """Return liquid phase reduced density (dens/critical dens) expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq

    def rho_liq(self, **kwargs):
        """Return liquid phase mass density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq * self.param.dens_mass_crit

    def rho_mol_liq(self, **kwargs):
        """Return liquid phase molar density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq * self.param.dens_mass_crit / self.param.mw

    def delta_vap(self, **kwargs):
        """Return vapor phase reduced density (dens/critical dens) expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap

    def rho_vap(self, **kwargs):
        """Return vapor phase mass density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap * self.param.dens_mass_crit

    def rho_mol_vap(self, **kwargs):
        """Return vapor phase molar density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap * self.param.dens_mass_crit / self.param.mw

    def cv_liq(self, **kwargs):
        """Return liquid phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cv = blk.func_cv(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return cv * self.param.uc["kJ/kg/K to J/mol/K"]
        return cv * self.param.uc["kJ/kg/K to J/kg/K"]

    def cv_vap(self, **kwargs):
        """Return vapor phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cv = blk.cv_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return cv * self.param.uc["kJ/kg/K to J/mol/K"]
        return cv * self.param.uc["kJ/kg/K to J/kg/K"]

    def cp_liq(self, **kwargs):
        """Return liquid phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cp = blk.func_cp(c, delta_liq, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return cp * self.param.uc["kJ/kg/K to J/mol/K"]
        return cp * self.param.uc["kJ/kg/K to J/kg/K"]

    def cp_vap(self, **kwargs):
        """Return vapor phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cp = blk.cp_func(c, delta_vap, tau)
        if self.amount_basis == AmountBasis.MOLE:
            return cp * self.param.uc["kJ/kg/K to J/mol/K"]
        return cp * self.param.uc["kJ/kg/K to J/kg/K"]

    def w_liq(self, **kwargs):
        """Return liquid phase speed of sound expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["w_func"])
        return blk.func_w(c, delta_liq, tau)

    def w_vap(self, **kwargs):
        """Return vapor phase speed of sound expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["w_func"])
        return blk.func_w(delta_vap, tau)

    def p_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_crit / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("p_sat expression requires either T or tau arg")
        self.add_funcs(names=["p_sat_func"])
        return (
            self.blk.p_sat_func(self.param.pure_component, tau)
            * self.param.uc["kPa to Pa"]
        )

    def delta_liq_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_crit / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("delta_liq_sat expression requires either T or tau arg")
        return self.blk.delta_sat_l_func(self.param.pure_component, tau)

    def delta_vap_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_crit / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("delta_vap_sat expression requires either T or tau arg")
        return self.blk.delta_sat_v_func(self.param.pure_component, tau)

    def T_sat(self, p):
        """Return saturation temperature as a function of p"""
        p *= self.param.uc_Pa_to_kPa
        return self.param.temperature_crit / self.blk.tau_sat_func(
            self.param.pure_component, p
        )

    def tau_sat(self, p):
        """Return saturation tau as a function of p"""
        p *= self.param.uc_Pa_to_kPa
        return self.blk.func_tau_sat(self.param.pure_component, p)


@declare_process_block_class("HelmholtzParameterBlock")
class HelmholtzParameterBlockData(PhysicalParameterBlock):
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
    CONFIG.declare(
        "phase_presentation",
        ConfigValue(
            default=PhaseType.MIX,
            domain=In(PhaseType),
            description="Set the way phases are presented to models",
            doc="""Set the way phases are presented to models. The MIX option
appears to the framework to be a mixed phase containing liquid and/or vapor.
The mixed option can simplify calculations at the unit model level since it can
be treated as a single phase, but unit models such as flash vessels will not be
able to treat the phases independently. The LG option presents as two separate
phases to the framework. The L or G options can be used if it is known for sure
that only one phase is present.
**default** - PhaseType.MIX
**Valid values:** {
**PhaseType.MIX** - Present a mixed phase with liquid and/or vapor,
**PhaseType.LG** - Present a liquid and vapor phase,
**PhaseType.L** - Assume only liquid can be present,
**PhaseType.G** - Assume only vapor can be present}""",
        ),
    )

    CONFIG.declare(
        "state_vars",
        ConfigValue(
            default=StateVars.PH,
            domain=In(StateVars),
            description="State variable set",
            doc="""The set of state variables to use. Depending on the use, one
state variable set or another may be better computationally. Usually pressure
and enthalpy are the best choice because they are well behaved during a phase
change.
**default** - StateVars.PH
**Valid values:** {
**StateVars.PH** - Pressure-Enthalpy,
**StateVars.TPX** - Temperature-Pressure-Quality}""",
        ),
    )

    CONFIG.declare(
        "ammount_basis",
        ConfigValue(
            default=AmountBasis.MOLE,
            domain=In(AmountBasis),
            description="State variable set",
            doc="""The set of state variables to use. Depending on the use, one
state variable set or another may be better computationally. Usually pressure
and enthalpy are the best choice because they are well behaved during a phase
change.
**default** - StateVars.PH
**Valid values:** {
**StateVars.PH** - Pressure-Enthalpy,
**StateVars.TPX** - Temperature-Pressure-Quality}""",
        ),
    )

    def htpx(self, T=None, p=None, x=None, units=None, amount_basis=AmountBasis.MOLE):
        """
        Convenience function to calculate enthalpy from temperature and either
        pressure or vapor fraction. This function can be used for inlet streams and
        initialization where temperature is known instead of enthalpy.
        User must provide values for one of these sets of values: {T, P}, {T, x},
        or {P, x}.
        Args:
            T: Temperature
            P: Pressure, None if saturated
            x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
               superheated or subcooled
        Returns:
            Total molar enthalpy [J/mol].
        """
        if units is None:
            if amount_basis==AmountBasis.MOLE:
                units=pyo.units.J/pyo.units.mol
            else:
                units=pyo.units.J/pyo.units.kg
        te = HelmholtzThermoExpressions(self, self, amount_basis=amount_basis)
        tmin = pyo.value(self.temperature_min)
        tmax = pyo.value(self.temperature_max)
        pmin = pyo.value(self.pressure_min)
        pmax = pyo.value(self.pressure_max)

        if not sum((p is None, T is None, x is None)) == 1:
            raise RuntimeError(
                "htpx function must be provided exaclty two of the arguments T, p, x"
            )
        if T is not None:
            T = pyo.units.convert(T, to_units=pyo.units.K)
            if not tmin <= pyo.value(T) <= tmax:
                raise RuntimeError(
                    f"T = {pyo.value(T)}, ({tmin} K <= T <= {tmax} K)"
                )
        if x is not None:
            if not 0 <= pyo.value(x) <= 1:
                raise RuntimeError(f"x = {pyo.value(x)}, (0 K <= x <= 1)")
        if p is not None:
            p = pyo.units.convert(p, to_units=pyo.units.Pa)
            if not pmin <= pyo.value(p) <= pmax:
                raise RuntimeError(
                    f"p = {pyo.value(p)}, ({pmin} kPa <= p <= {pmax} kPa)"
                )
            if T is not None:
                # P, T may be underspecified, but assume you know it's clearly a
                # vapor or liquid so figure out which and set x.
                psat = te.p_sat(T)
                if pyo.value(p) < pyo.value(psat):
                    x = 1
                else:
                    x = 0
        return pyo.value(pyo.units.convert(te.h(T=T, p=p, x=x), units))


    def _set_default_scaling(self):
        """Set default scaling parameters to be used if not otherwise set"""

        # TODO<jce> leaving flow to not break things, but plan to remove it and
        #     make the user specify with no default.
        self.set_default_scaling("flow_mol", 1e-4)
        self.set_default_scaling("flow_mol_comp", 1e-4)
        self.set_default_scaling("flow_vol", 100)
        self.set_default_scaling("flow_mass", 1)

        # Set some scalings with reasonable a priori values
        self.set_default_scaling("temperature_crit", 1e-2)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("temperature", 1e-1)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("vapor_frac", 1e1)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)  # Not liq
        self.set_default_scaling("pressure_crit", 1e-6)
        self.set_default_scaling("dens_mass_crit", 1e-2)
        self.set_default_scaling("mw", 1e3)
        self.set_default_scaling("temperature_sat", 1e-2)
        self.set_default_scaling("temperature_red", 1)
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)
        self.set_default_scaling("dens_phase_red", 1)
        self.set_default_scaling("enth_mol_sat_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mol_sat_phase", 1e-4, index="Vap")
        self.set_default_scaling("dh_vap_mol", 1e-4)
        self.set_default_scaling("energy_internal_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("energy_internal_mol_phase", 1e-4)
        self.set_default_scaling("enth_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("entr_mol_phase", 1e-1, index="Liq")
        self.set_default_scaling("entr_mol_phase", 1e-1, index="Vap")
        self.set_default_scaling("cp_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("cp_mol_phase", 1e-2, index="Vap")
        self.set_default_scaling("cv_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("cv_mol_phase", 1e-2, index="Vap")
        self.set_default_scaling("dens_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("phase_frac", 10)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("energy_internal_mol", 1e-3)
        self.set_default_scaling("entr_mol", 1e-1)
        self.set_default_scaling("cp_mol", 1e-2)
        self.set_default_scaling("cv_mol", 1e-2)
        self.set_default_scaling("dens_mass", 1)
        self.set_default_scaling("dens_mol", 1e-3)
        self.set_default_scaling("heat_capacity_ratio", 1e1)
        self.set_default_scaling("enth_mass", 1)

    def _create_component_and_phase_objects(self):
        # Create chemical component objects
        for c in self.component_list:
            setattr(self, str(c), Component(_component_list_exists=True))
        # Create phase objects
        self.private_phase_list = pyo.Set(initialize=["Vap", "Liq"])
        # Make mixed phase object for mixed phase option
        if self.config.phase_presentation == PhaseType.MIX:
            self.Mix = Phase()
        # Make liquid phase object for L and LG options
        if (
            self.config.phase_presentation == PhaseType.LG
            or self.config.phase_presentation == PhaseType.L
        ):
            self.Liq = LiquidPhase()
        # Make vapor phase object for G and LG options
        if (
            self.config.phase_presentation == PhaseType.LG
            or self.config.phase_presentation == PhaseType.G
        ):
            self.Vap = VaporPhase()
        # State var set
        self.state_vars = self.config.state_vars

    def build(self):
        super().build()
        # set the component_list as required for the generic IDAES properties
        self.component_list = pyo.Set(initialize=[self.config.pure_component])
        # sinice this a only a pure component package, have a specific
        # pure_component attirbute
        self.pure_component = self.config.pure_component
        # Add default scaling factors for property expressions or variables
        self._set_default_scaling()
        # Add idaes component and phase objects
        self._create_component_and_phase_objects()
        # To ensure consistency pull parameters from external functions
        add_helmholtz_external_functions(  # Add parameter external functions
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
        # functions more than once, so define Pyomo parameters with the values,
        # use the external function here to avoid defining the parameters twice
        pu = pyo.units
        cmp = self.pure_component
        self.add_param(
            "mw",
            pu.convert(self.mw_func(cmp), pu.kg / pu.mol),
        )
        self.add_param(
            "sgc",
            pyo.units.convert(self.sgc_func(cmp), pu.J / pu.kg / pu.K),
        )
        self.add_param(
            "sgc_mol",
            pu.convert(self.sgc * self.mw, pu.J / pu.mol / pu.K),
        )
        self.add_param(
            "pressure_crit",
            pu.convert(self.pc_func(cmp), pu.Pa),
        )
        self.add_param(
            "pressure_trip",
            pu.convert(self.pt_func(cmp), pu.Pa),
        )
        self.add_param(
            "pressure_min",
            pu.convert(self.pmin_func(cmp), pu.Pa),
        )
        self.add_param(
            "pressure_max",
            pu.convert(self.pmax_func(cmp), pu.Pa),
        )
        self.add_param(
            "temperature_crit",
            pu.convert(self.tc_func(cmp), pu.K),
        )
        self.add_param(
            "temperature_trip",
            pu.convert(self.tt_func(cmp), pu.K),
        )
        self.add_param(
            "temperature_min",
            pu.convert(self.tmin_func(cmp), pu.K),
        )
        self.add_param(
            "temperature_max",
            pu.convert(self.tmax_func(cmp), pu.K),
        )
        self.add_param(
            "dense_mass_crit",
            pu.convert(self.rhoc_func(cmp), pu.kg / pu.m**3),
        )
        self.add_param(
            "dense_mol_crit",
            pu.convert(self.dense_mass_crit / self.mw, pu.mol / pu.m**3),
        )

        self.uc = {
            "J/mol to kJ/kg": (pyo.units.kJ / 1000 / pyo.units.J) / self.mw,
            "J/mol to J/kg": 1.0 / self.mw,
            "kJ/kg to J/mol": (pyo.units.J * 1000 / pyo.units.kJ) * self.mw,
            "J/mol/K to kJ/kg/K": (pyo.units.kJ / 1000 / pyo.units.J) / self.mw,
            "J/mol/K to J/kg/K": 1.0 / self.mw,
            "kJ/kg/K to J/mol/K": (pyo.units.J * 1000 / pyo.units.kJ) * self.mw,
            "J/kg to kJ/kg": (pyo.units.kJ / 1000 / pyo.units.J),
            "kJ/kg to J/kg": (pyo.units.J * 1000 / pyo.units.kJ),
            "J/kg/K to kJ/kg/K": (pyo.units.kJ / 1000 / pyo.units.J),
            "kJ/kg/K to J/kg/K": (pyo.units.J * 1000 / pyo.units.kJ),
            "kPa to Pa": (pyo.units.Pa * 1000 / pyo.units.kPa),
            "Pa to kPa": (pyo.units.kPa / 1000 / pyo.units.Pa),
        }

        if self.config.state_vars == StateVars.TPX:
            self._tpx_complimentarity_form_smoothing()

    def _tpx_complimentarity_form_smoothing(self):
        # Smoothing arameters for
        self.smoothing_pressure_over = Param(
            mutable=True,
            initialize=1e-4,
            doc="Smooth max parameter (pressure over)",
            units=pyo.units.kPa,
        )

        self.smoothing_pressure_under = Param(
            mutable=True,
            initialize=1e-4,
            doc="Smooth max parameter (pressure under)",
            units=pyo.units.kPa,
        )

    def add_param(self, name, expr):
        self.add_component(
            name,
            pyo.Param(
                initialize=pyo.value(expr),
                units=pyo.units.get_units(expr),
            ),
        )

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
            ],
        )

        pc = pyo.value(self.pressure_crit)
        pt = pyo.value(self.pressure_trip)
        tc = pyo.value(self.temperature_crit)
        tt = pyo.value(self.temperature_trip)
        pmax = pyo.value(self.pressure_max)

        tau_sat_vec = np.linspace(
            1,
            pyo.value(self.temperature_crit) / (pyo.value(self.temperature_trip)),
            200,
        )
        p_sat_vec = [None] * len(tau_sat_vec)
        delta_sat_v_vec = [None] * len(tau_sat_vec)
        delta_sat_l_vec = [None] * len(tau_sat_vec)
        h_sat_v_vec = [None] * len(tau_sat_vec)
        h_sat_l_vec = [None] * len(tau_sat_vec)

        for i, tau in enumerate(tau_sat_vec):
            p_sat_vec[i] = pyo.value(self.p_sat_func(self.pure_component, tau))
            delta_sat_l_vec[i] = pyo.value(
                self.delta_sat_l_func(self.pure_component, tau)
            )
            delta_sat_v_vec[i] = pyo.value(
                self.delta_sat_v_func(self.pure_component, tau)
            )
            h_sat_v_vec[i] = pyo.value(
                self.h_func(self.pure_component, delta_sat_v_vec[i], tau)
            )
            h_sat_l_vec[i] = pyo.value(
                self.h_func(self.pure_component, delta_sat_l_vec[i], tau)
            )

        plt.yscale("log")
        plt.plot(h_sat_l_vec, p_sat_vec, c="b", label="sat liquid")
        plt.plot(h_sat_v_vec, p_sat_vec, c="r", label="sat vapor")

        t_vec = np.linspace(
            pyo.value(self.temperature_trip), pyo.value(self.temperature_crit), 6
        )
        # t_vec = [340]
        p = {}
        h_l = {}
        h_v = {}

        # plot isotherms in sat region
        for t in t_vec:
            tau = pyo.value(self.temperature_crit) / t
            p[t] = pyo.value(self.p_sat_func(self.pure_component, tau))
            delta_l = pyo.value(self.delta_sat_l_func(self.pure_component, tau))
            delta_v = pyo.value(self.delta_sat_v_func(self.pure_component, tau))
            h_v[t] = pyo.value(self.h_func(self.pure_component, delta_v, tau))
            h_l[t] = pyo.value(self.h_func(self.pure_component, delta_l, tau))
            plt.plot([h_l[t], h_v[t]], [p[t], p[t]], c="g")
            plt.text(h_l[t] / 2 + h_v[t] / 2, p[t], f"T = {t} K", ha="center")

        for t in t_vec:
            tau = pyo.value(self.temperature_crit) / t
            p_vec = np.logspace(np.log10(p[t]), np.log10(pmax / 1000), 50)
            h_vec = [None] * len(p_vec)
            for i, pv in enumerate(p_vec):
                h_vec[i] = pyo.value(self.hlpt_func(self.pure_component, pv, tau))
            plt.plot(h_vec, p_vec, c="g")

        for t in t_vec:
            tau = pyo.value(self.temperature_crit) / t
            p_vec = np.logspace(np.log10(pt / 1000 / 10), np.log10(p[t]), 50)

            h_vec = [None] * len(p_vec)
            for i, pv in enumerate(p_vec):
                h_vec[i] = pyo.value(self.hvpt_func(self.pure_component, pv, tau))
            plt.plot(h_vec, p_vec, c="g")

        deltat_l = pyo.value(self.delta_liq_func(self.pure_component, pt, tc / tt))
        hc = pyo.value(self.h_func(self.pure_component, 1, 1))
        ht = pyo.value(self.h_func(self.pure_component, deltat_l, tc / tt))
        plt.scatter([hc], [pc / 1000])
        plt.scatter([ht], [pt / 1000])

        plt.title(f"P-H Diagram for {self.pure_component}")
        plt.xlabel("Enthalpy (kJ/kg)")
        plt.ylabel("Pressure (kPa)")

        plt.show()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "temperature_crit": {"method": None, "units": "K"},
                "pressure_crit": {"method": None, "units": "Pa"},
                "dens_mass_crit": {"method": None, "units": "kg/m^3"},
                "specific_gas_constant": {"method": None, "units": "J/kg.K"},
                "mw": {"method": None, "units": "kg/mol"},
                "temperature_sat": {"method": "None", "units": "K"},
                "flow_mol": {"method": None, "units": "mol/s"},
                "flow_mass": {"method": None, "units": "kg/s"},
                "flow_vol": {"method": None, "units": "m^3/s"},
                "temperature": {"method": None, "units": "K"},
                "pressure": {"method": None, "units": "Pa"},
                "vapor_frac": {"method": None, "units": None},
                "dens_mass_phase": {"method": None, "units": "kg/m^3"},
                "temperature_red": {"method": None, "units": None},
                "pressure_sat": {"method": None, "units": "kPa"},
                "energy_internal_mol_phase": {"method": None, "units": "J/mol"},
                "enth_mol_phase": {"method": None, "units": "J/mol"},
                "entr_mol_phase": {"method": None, "units": "J/mol.K"},
                "cp_mol_phase": {"method": None, "units": "J/mol.K"},
                "cv_mol_phase": {"method": None, "units": "J/mol.K"},
                "speed_sound_phase": {"method": None, "units": "m/s"},
                "dens_mol_phase": {"method": None, "units": "mol/m^3"},
                "therm_cond_phase": {"method": None, "units": "W/m.K"},
                "visc_d_phase": {"method": None, "units": "Pa.s"},
                "visc_k_phase": {"method": None, "units": "m^2/s"},
                "phase_frac": {"method": None, "units": None},
                "flow_mol_comp": {"method": None, "units": "mol/s"},
                "energy_internal_mol": {"method": None, "units": "J/mol"},
                "enth_mol": {"method": None, "units": "J/mol"},
                "entr_mol": {"method": None, "units": "J/mol.K"},
                "cp_mol": {"method": None, "units": "J/mol.K"},
                "cv_mol": {"method": None, "units": "J/mol.K"},
                "heat_capacity_ratio": {"method": None, "units": None},
                "dens_mass": {"method": None, "units": "kg/m^3"},
                "dens_mol": {"method": None, "units": "mol/m^3"},
                "dh_vap_mol": {"method": None, "units": "J/mol"},
            }
        )

        obj.add_default_units(
            {
                "time": pyo.units.s,
                "length": pyo.units.m,
                "mass": pyo.units.kg,
                "amount": pyo.units.mol,
                "temperature": pyo.units.K,
            }
        )

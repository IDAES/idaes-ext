#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""Generic Helmholtz EOS StateBlock Class
"""
__author__ = "John Eslick"


import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.core.base.units_container import InconsistentUnitsError
from pyomo.common.fileutils import find_library
from pyomo.common.config import ConfigValue, In

from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from idaes.core import (
    StateBlock,
    StateBlockData,
)
from helmholtz_functions import (
    add_helmholtz_external_functions,
    HelmholtzParameterBlock,
    HelmholtzThermoExpressions,
    AmountBasis,
    available,
)


class _StateBlock(StateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    @staticmethod
    def _set_fixed(v, f):
        if f:
            v.fix()
        else:
            v.unfix()

    @staticmethod
    def _set_not_fixed(v, state, key, hold):
        if state is not None:
            if not v.fixed:
                try:
                    v.value = state[key]
                except KeyError:
                    pass
            if hold:
                v.fix()

    def initialize(self, *args, **kwargs):
        flags = {}
        hold_state = kwargs.pop("hold_state", False)
        state_args = kwargs.pop("state_args", None)
        for i, v in self.items():
            pp = self[i].config.parameters.config.phase_presentation
            sv = self[i].state_vars
            ab = self[i].amount_basis
            if sv == StateVars.PH and ab == AmountBasis.MOLE:
                flags[i] = (v.flow_mol.fixed, v.enth_mol.fixed, v.pressure.fixed)
                _set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                _set_not_fixed(v.enth_mol, state_args, "enth_mol", hold_state)
                _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PH and ab == AmountBasis.MASS:
                flags[i] = (v.flow_mass.fixed, v.enth_mass.fixed, v.pressure.fixed)
                _set_not_fixed(v.flow_mass, state_args, "flow_mass", hold_state)
                _set_not_fixed(v.enth_mass, state_args, "enth_mass", hold_state)
                _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PS and ab == AmountBasis.MOLE:
                flags[i] = (v.flow_mol.fixed, v.entr_mol.fixed, v.pressure.fixed)
                _set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                _set_not_fixed(v.entr_mol, state_args, "enth_mol", hold_state)
                _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PS and ab == AmountBasis.MASS:
                flags[i] = (v.flow_mass.fixed, v.entr_mass.fixed, v.pressure.fixed)
                _set_not_fixed(v.flow_mass, state_args, "flow_mass", hold_state)
                _set_not_fixed(v.entr_mass, state_args, "enth_mass", hold_state)
                _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PU and ab == AmountBasis.MOLE:
                flags[i] = (v.flow_mol.fixed, v.energy_internal_mol.fixed, v.pressure.fixed)
                _set_not_fixed(v.energy_internal_mol, state_args, "energy_internal_mol", hold_state)
                _set_not_fixed(v.energy_internal_mol, state_args, "energy_internal_mol", hold_state)
                _set_not_fixed(v.energy_internal_mol, state_args, "energy_internal_mol", hold_state)
            elif sv == StateVars.PU and ab == AmountBasis.MASS:
                flags[i] = (v.flow_mass.fixed, v.energy_internal_mass.fixed, v.pressure.fixed)
                _set_not_fixed(v.energy_internal_mass, state_args, "energy_internal_mass", hold_state)
                _set_not_fixed(v.energy_internal_mass, state_args, "energy_internal_mass", hold_state)
                _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.TPX and ab == AmountBasis.MOLE:
                # Hold the T-P-x vars
                if pp in (PhaseType.MIX, PhaseType.LG):
                    flags[i] = (
                        v.flow_mol.fixed,
                        v.temperature.fixed,
                        v.pressure.fixed,
                        v.vapor_frac.fixed,
                    )
                    _set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                    _set_not_fixed(v.temperature, state_args, "temperaure", hold_state)
                    _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
                    _set_not_fixed(v.vapor_frac, state_args, "vapor_frac", hold_state)
                else:
                    flags[i] = (v.flow_mol.fixed, v.temperature.fixed, v.pressure.fixed)
                    _set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                    _set_not_fixed(v.temperature, state_args, "temperaure", hold_state)
                    _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.TPX and ab == AmountBasis.MASS:
                # Hold the T-P-x vars
                if pp in (PhaseType.MIX, PhaseType.LG):
                    flags[i] = (
                        v.flow_mass.fixed,
                        v.temperature.fixed,
                        v.pressure.fixed,
                        v.vapor_frac.fixed,
                    )
                    _set_not_fixed(v.flow_mass, state_args, "flow_mass", hold_state)
                    _set_not_fixed(v.temperature, state_args, "temperaure", hold_state)
                    _set_not_fixed(v.pressure, state_args, "pressure", hold_state)
                    _set_not_fixed(v.vapor_frac, state_args, "vapor_frac", hold_state)
                else:
                    flags[i] = (v.flow_mol.fixed, v.temperature.fixed, v.pressure.fixed)
                    _set_not_fixed(v.flow_mass, state_args, "flow_mass", hold_state)
                    _set_not_fixed(v.temperature, state_args, "temperaure", hold_state)
                    _set_not_fixed(v.pressure, state_args, "pressure", hold_state)

        # Call initialize on each data element
        for i in self:
            self[i].initialize(*args, **kwargs)
        return flags

    def release_state(self, flags, **kwargs):
        for i, f in flags.items():
            pp = self[i].config.parameters.config.phase_presentation
            sv = self[i].state_vars
            ab = self[i].amount_basis
            if sv == StateVars.PH and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].enth_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PH and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].enth_mass, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PS and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].entr_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PS and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].entr_mass, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PU and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].energy_internal_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PU and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].energy_internal_mass, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.TPX and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].temperature, f[1])
                self._set_fixed(self[i].pressure, f[2])
                if pp in (PhaseType.MIX, PhaseType.LG):
                    self._set_fixed(self[i].vapor_frac, f[3])
            elif sv == StateVars.TPX and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].temperature, f[1])
                self._set_fixed(self[i].pressure, f[2])
                if pp in (PhaseType.MIX, PhaseType.LG):
                    self._set_fixed(self[i].vapor_frac, f[3])


class HelmholtzStateBlockData(StateBlockData):
    """
    This is a base clase for Helmholtz equations of state using IDAES standard
    Helmholtz EOS external functions written in C++.
    """

    def initialize(self, *args, **kwargs):
        # This particular property package has no need for initialization
        pass

    def _state_vars(self):
        """Create the state variables"""
        params = self.config.parameters
        cmp = params.pure_component
        # Add flow variable
        if self.amount_basis == AmountBasis.MOLE:
            self.flow_mol = Var(
                initialize=1, doc="Total mole flow", units=pyunits.mol / pyunits.s
            )
            self.flow_mass = Expression(
                expr=self.mw * self.flow_mol, doc="Total mass flow"
            )
        else:
            self.flow_mass = Var(
                initialize=1, doc="Total mass flow", units=pyunits.kg / pyunits.s
            )
            self.flow_mol = Expression(
                expr=self.flow_mass / self.mw, doc="Total mole flow"
            )
        # All supported state variable sets include pressure
        self.pressure = Var(
            domain=PositiveReals,
            initialize=params.default_pressure_value,
            doc="Pressure",
            bounds=params.default_pressure_bounds,
            units=pyunits.Pa,
        )
        P = self.pressure * params.uc["Pa to kPa"]
        # If it's single phase use provide fixed expressions for vapor frac
        if phase_set == PhaseType.L:
            self.vapor_frac = Expression(
                expr=0.0, doc="Vapor mole fraction (mol vapor/mol total)"
            )
        elif phase_set == PhaseType.G:
            self.vapor_frac = Expression(
                expr=1.0,
                doc="Vapor mole fraction (mol vapor/mol total)",
            )
        # Add other state var
        if self.state_vars == StateVars.PH and self.amount_basis == AmountBasis.MOLE:
            self.enth_mol = Var(
                initialize=params.default_enthalpy_mol_value,
                doc="Total molar enthalpy",
                bounds=params.default_enthalpy_mol_bounds,
                units=pyunits.J / pyunits.mol,
            )
            h_mass = self.enth_mol * params.uc["J/mol to kJ/kg"]
            self.temperature = Expression(
                expr=self.temperature_crit / self.func_tau(cmp, h_mass, P), doc="Temperature"
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vf(cmp, h_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "enth_mol": self.enth_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.enth_mol, self.pressure))
        elif self.state_vars == StateVars.PH and self.amount_basis == AmountBasis.MASS:
            self.enth_mass = Var(
                initialize=params.default_enthalpy_mass_value,
                doc="Total enthalpy per mass",
                bounds=params.default_enthalpy_mass_bounds,
                units=pyunits.J / pyunits.kg,
            )
            h_mass = self.enth_mass * params.uc["J/kg to kJ/kg"]
            self.temperature = Expression(
                expr=self.temperature_crit / self.func_tau(cmp, h_mass, P), doc="Temperature"
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vf(cmp, h_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mass": self.flow_mol,
                "enth_mass": self.enth_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mass,))
            self.intensive_set = ComponentSet((self.enth_mass, self.pressure))
        elif self.state_vars == StateVars.PS and self.amount_basis == AmountBasis.MOLE:
            self.entr_mol = Var(
                initialize=params.default_entropy_mol_value,
                doc="Total molar entropy",
                bounds=params.default_entropy_mol_bounds,
                units=pyunits.J / pyunits.mol / pyunits.K,
            )
            s_mass = self.entr_mol * params.uc["J/mol/K to kJ/kg/K"]
            self.temperature = Expression(
                expr=self.temperature_crit / self.func_taus(cmp, s_mass, P), doc="Temperature"
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vfs(cmp, s_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "entr_mol": self.entr_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.entr_mol, self.pressure))
        elif self.state_vars == StateVars.PS and self.amount_basis == AmountBasis.MASS:
            self.entr_mass = Var(
                initialize=params.default_entropy_mass_value,
                doc="Total entropy per mass",
                bounds=params.default_entropy_mass_bounds,
                units=pyunits.J / pyunits.kg / pyunits.K,
            )
            s_mass = self.entr_mass * params.uc["J/kg/K to kJ/kg/K"]
            self.temperature = Expression(
                expr=self.temperature_crit / self.func_taus(cmp, s_mass, P), doc="Temperature"
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vfs(cmp, s_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mass,
                "entr_mol": self.entr_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mass,))
            self.intensive_set = ComponentSet((self.entr_mass, self.pressure))
        elif self.state_vars == StateVars.PU and self.amount_basis == AmountBasis.MOLE:
            self.energy_internal_mol = Var(
                initialize=params.default_internal_energy_mol_value,
                doc="Total molar internal energy",
                bounds=params.default_internal_energy_mol_bounds,
                units=pyunits.J / pyunits.mol,
            )
            u_mass = self.energy_internal_mol *  params.uc["J/mol to kJ/kg"]
            self.temperature = Expression(
                expr=self.temperature_crit / self.func_tauu(cmp, u_mass, P), doc="Temperature"
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vfu(cmp, s_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "energy_internal_mol": self.energy_internal_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.energy_internal_mol, self.pressure))
        elif self.state_vars == StateVars.PU and self.amount_basis == AmountBasis.MASS:
            self.energy_internal_mass = Var(
                initialize=params.default_internal_energy_mass_value,
                doc="Total internal energy per mass",
                bounds=params.default_internal_energy_mass_bounds,
                units=pyunits.J / pyunits.kg,
            )
            u_mass = self.energy_internal_mass * params.uc["J/kg to kJ/kg"]
            self.temperature = Expression(
                expr=self.temperature_crit / self.func_tauu(cmp, u_mass, P), doc="Temperature"
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = Expression(
                    expr=self.func_vfu(cmp, s_mass, P),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mass": self.flow_mass,
                "energy_internal_mass": self.energy_internal_mass,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mass,))
            self.intensive_set = ComponentSet((self.energy_internal_mass, self.pressure))
        if self.state_vars == StateVars.TPX:
            self.temperature = Var(
                domain=PositiveReals,
                initialize=params.default_temperature_value,
                doc="Temperature",
                bounds=params.default_temperature_bounds,
                units=pyunits.K,
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG
                self.vapor_frac = Var(
                    initialize=0.0,
                    units=pyunits.dimensionless,
                    doc="Vapor fraction"
                    # No bounds here, since it is often (usually) on it's bound
                    # and that's not the best for IPOPT
                )
                self.intensive_set = ComponentSet(
                    (self.temperature, self.pressure, self.vapor_frac)
                )
                self._state_vars_dict = {
                    "temperature": self.temperature,
                    "pressure": self.pressure,
                    "vapor_frac": self.vapor_frac,
                }
            else:
                self.intensive_set = ComponentSet((self.temperature, self.pressure))
                self._state_vars_dict = {
                    "temperature": self.temperature,
                    "pressure": self.pressure,
                }
            if self.amount_basis == AmountBasis.MOLE:
                self.extensive_set = ComponentSet((self.flow_mol,))
                _state_vars_dict["flow_mol"] = self.flow_mol
            else:
                self.extensive_set = ComponentSet((self.flow_mass,))
                _state_vars_dict["flow_mass"] = self.flow_mass


    def _tpx_phase_eq(self):
        # Saturation pressure
        params = self.config.parameters
        cmp = params.pure_component
        eps_pu = params.smoothing_pressure_under
        eps_po = params.smoothing_pressure_over
        priv_plist = params.private_phase_list
        plist = params.phase_list
        rhoc = params.dens_mass_crit
        # Convert presssures to kPa for external functions and nicer scaling in
        # the compimentarity-type constraints.
        P = self.pressure * params.uc["Pa to kPa"]
        Psat = self.pressure_sat * params.uc["Pa to kPa"]
        vf = self.vapor_frac
        tau = self.tau
        # Terms for determining if you are above, below, or at the Psat
        self.P_under_sat = Expression(
            expr=smooth_max(0, Psat - P, eps_pu),
            doc="pressure above Psat, 0 if liqid exists [kPa]",
        )
        self.P_over_sat = Expression(
            expr=smooth_max(0, P - Psat, eps_po),
            doc="pressure below Psat, 0 if vapor exists [kPa]",
        )
        # Calculate liquid and vapor density.  If the phase doesn't exist,
        # density will be calculated at the saturation or critical pressure
        def rule_dens_mass(b, p):
            if p == "Liq":
                return rhoc * self.func_delta_liq(cmp, P + self.P_under_sat, tau)
            else:
                return rhoc * self.func_delta_vap(cmp, P - self.P_over_sat, tau)

        self.dens_mass_phase = Expression(priv_plist, rule=rule_dens_mass)
        # Reduced Density (no _mass_ identifier because mass or mol is same)
        def rule_dens_red(b, p):
            return self.dens_mass_phase[p] / rhoc
        self.dens_phase_red = Expression(
            priv_plist, rule=rule_dens_red, doc="reduced density"
        )
        if not self.config.defined_state:
            self.eq_complementarity = Constraint(
                expr=0 == (vf * self.P_over_sat - (1 - vf) * self.P_under_sat)
            )
        # eq_sat can activated to force the pressure to be the saturation
        # pressure, if you use this constraint deactivate eq_complementarity
        self.eq_sat = Constraint(expr=P / 1000.0 == Psat / 1000.0)
        self.eq_sat.deactivate()

    def build(self, *args):
        """
        Callable method for Block construction
        """
        # Call base class build
        super().build(*args)
        # Short path to the parameter block
        params = self.config.parameters
        self.params = params
        # Check if the library is available, and add external functions.
        self.available = params.available
        if not self.available:
            _log.error(
                "Shared lib 'general_helmholtz_external' not found. Is it installed?")
        add_helmholtz_external_functions(self)
        cmp = params.pure_component
        # Which state vars to use
        self.state_vars = params.state_vars
        # private phase list
        phlist = params.private_phase_list
        # Public phase list
        pub_phlist = params.phase_list
        component_list = params.component_list
        phase_set = params.config.phase_presentation
        self.phase_equilibrium_list = params.phase_equilibrium_list
        # Expressions that link to some parameters in the param block, which
        # are commonly needed, this lets you get the parameters with scale
        # factors directly from the state block
        self.temperature_crit = Expression(
            expr=params.temperature_crit, doc="critical temperaure"
        )
        self.pressure_crit = Expression(
            expr=params.pressure_crit, doc="critical pressure"
        )
        self.dens_mass_crit = Expression(
            expr=params.dens_mass_crit, doc="critical mass density"
        )
        self.dens_mol_crit = Expression(
            expr=params.dens_mass_crit/params.mw, doc="critical mole density"
        )
        self.mw = Expression(expr=params.mw, doc="molecular weight")
        # create the appropriate state variables and expressions for anything
        # that could be a state variable (H, S, U, P, T, X) on a mass and mole
        # basis if using TPx as state variables it also adds the complimentarity
        # constraints. Beyond this everything else is common.
        self._state_vars()

        # Some parameters/variables show up in several expressions, so to
        # enhance readability and compactness, give them short aliases
        Tc = params.temperature_crit
        rhoc = params.dens_mass_crit
        mw = self.mw
        # P is pressure in kPa for external function calls
        P = self.pressure * params.uc["Pa to kPa"]
        T = self.temperature
        vf = self.vapor_frac
        # Saturation temperature expression
        self.temperature_sat = Expression(
            expr=Tc / self.func_tau_sat(cmp, P), doc="Stauration temperature"
        )
        # Saturation tau (tau = Tc/T)
        self.tau_sat = Expression(expr=self.func_tau_sat(cmp, P))
        # Reduced temperature
        self.temperature_red = Expression(
            expr=T / Tc, doc="reduced temperature T/Tc"
        )
        self.tau = Expression(expr=Tc / T, doc="Tc/T")
        tau = self.tau
        # Saturation pressure
        self.pressure_sat = Expression(
            expr=self.func_p_sat(cmp, tau) * params.uc["kPa to Pa"],
            doc="Saturation pressure",
        )

        if self.state_vars != StateVars.TPX or len(phlist) == 1:
            # If we aren't using the TPX with complimentarity, just directly
            # calculate the density and reduced density, this includes single
            # phase TPX
            def rule_dens_mass(b, p):
                if p == "Liq":
                    return rhoc * self.func_delta_liq(cmp, P, tau)
                else:
                    return rhoc * self.func_delta_vap(cmp, P, tau)
            self.dens_mass_phase = Expression(
                phlist, rule=rule_dens_mass, doc="Mass density by phase"
            )
            # Reduced Density (no _mass_ identifier as mass or mol is same)
            def rule_dens_red(b, p):
                return self.dens_mass_phase[p] / rhoc
            self.dens_phase_red = Expression(
                phlist, rule=rule_dens_red, doc="Reduced density"
            )
        else:
            self._tpx_phase_eq()

        # the reduced density (delta) and Tc/T (tau) we can get rest of props
        delta = self.dens_phase_red

        # Saturated Enthalpy molar
        def rule_enth_mol_sat_phase(b, p):
            if p == "Liq":
                return (
                    self.func_hlpt(cmp, P, self.tau_sat) * params.uc["kJ/kg to J/mol"]
                )
            else:
                return (
                    self.func_hvpt(cmp, P, self.tau_sat) * params.uc["kJ/kg to J/mol"]
                )
        self.enth_mol_sat_phase = Expression(
            phlist,
            rule=rule_enth_mol_sat_phase,
            doc="Saturated enthalpy of the phases at pressure",
        )
        # Saturated Enthalpy mass
        def rule_enth_mass_sat_phase(b, p):
            if p == "Liq":
                return (
                    self.func_hlpt(cmp, P, self.tau_sat) * params.uc["kJ/kg to J/kg"]
                )
            else:
                return (
                    self.func_hvpt(cmp, P, self.tau_sat) * params.uc["kJ/kg to J/kg"]
                )
        self.enth_mass_sat_phase = Expression(
            phlist,
            rule=rule_enth_mass_sat_phase,
            doc="Saturated enthalpy of the phases at pressure",
        )
        # Saturated Entropy molar
        def rule_entr_mol_sat_phase(b, p):
            if p == "Liq":
                return (
                    self.func_slpt(cmp, P, self.tau_sat) * params.uc["kJ/kg/K to J/mol/K"]
                )
            else:
                return (
                    self.func_svpt(cmp, P, self.tau_sat) * params.uc["kJ/kg/K to J/mol/K"]
                )
        self.entr_mol_sat_phase = Expression(
            phlist,
            rule=rule_entr_mol_sat_phase,
            doc="Saturated entropy of the phases at pressure",
        )
        # Saturated Entropy mass
        def rule_entr_mass_sat_phase(b, p):
            if p == "Liq":
                return (
                    self.func_slpt(cmp, P, self.tau_sat) * params.uc["kJ/kg/K to J/kg/K"]
                )
            else:
                return (
                    self.func_svpt(cmp, P, self.tau_sat) * params.uc["kJ/kg/K to J/kg/K"]
                )
        self.entr_mass_sat_phase = Expression(
            phlist,
            rule=rule_entr_mol_sat_phase,
            doc="Saturated entropy of the phases at pressure",
        )

        if len(phlist) == 2:
            # delta h vap at P
            self.dh_vap_mol = Expression(
                expr=(self.enth_mol_sat_phase["Vap"] - self.enth_mol_sat_phase["Liq"]),
                doc="Enthaply of vaporization at pressure and saturation temperature",
            )
            # delta s vap at P
            self.ds_vap_mol = Expression(
                expr=(self.entr_mol_sat_phase["Vap"] - self.entr_mol_sat_phase["Liq"]),
                doc="Entropy of vaporization at pressure and saturation temperature",
            )
            # delta h vap at P
            self.dh_vap_mass = Expression(
                expr=(self.enth_mass_sat_phase["Vap"] - self.enth_mass_sat_phase["Liq"]),
                doc="Enthaply of vaporization at pressure and saturation temperature",
            )
            # delta s vap at P
            self.ds_vap_mol = Expression(
                expr=(self.entr_mass_sat_phase["Vap"] - self.entr_mass_sat_phase["Liq"]),
                doc="Entropy of vaporization at pressure and saturation temperature",
            )


        # Phase fraction
        def rule_phase_frac(b, p):
            if p == "Vap":
                return vf
            elif p == "Liq":
                return 1.0 - vf
        self.phase_frac = Expression(
            phlist, rule=rule_phase_frac, doc="Phase fraction"
        )
        # Phase Internal Energy
        def rule_energy_internal_mol_phase(b, p):
            return self.func_u(cmp, delta[p], tau) * params.uc["kJ/kg to J/mol"]
        self.energy_internal_mol_phase = Expression(
            phlist,
            rule=rule_energy_internal_mol_phase,
            doc="Phase internal energy",
        )
        def rule_energy_internal_mass_phase(b, p):
            return self.func_u(cmp, delta[p], tau) * params.uc["kJ/kg to J/kg"]
        self.energy_internal_mass_phase = Expression(
            phlist,
            rule=rule_energy_internal_mass_phase,
            doc="Phase internal energy",
        )
        # Phase Enthalpy
        def rule_enth_mol_phase(b, p):
            return self.func_h(cmp, delta[p], tau) * params.uc["kJ/kg to J/mol"]
        self.enth_mol_phase = Expression(
            phlist,
            rule=rule_enth_mol_phase,
            doc="Phase enthalpy",
        )
        def rule_enth_mass_phase(b, p):
            return self.func_h(cmp, delta[p], tau) * params.uc["kJ/kg to J/kg"]
        self.enth_mass_phase = Expression(
            phlist,
            rule=rule_enth_mass_phase,
            doc="Phase enthalpy",
        )
        # Phase Entropy
        def rule_entr_mol_phase(b, p):
            return self.func_s(cmp, delta[p], tau) * params.uc["kJ/kg/K to J/mol/K"]
        self.entr_mol_phase = Expression(
            phlist,
            rule=rule_entr_mol_phase,
            doc="Phase entropy",
        )
        def rule_entr_mass_phase(b, p):
            return self.func_s(cmp, delta[p], tau) * params.uc["kJ/kg/K to J/kg/K"]
        self.entr_mol_phase = Expression(
            phlist,
            rule=rule_entr_mass_phase,
            doc="Phase entropy",
        )
        # Phase constant pressure heat capacity, cp
        def rule_cp_mol_phase(b, p):
            return self.func_cp(cmp, delta[p], tau) * params.uc["kJ/kg/K to J/mol/K"]
        self.cp_mol_phase = Expression(
            phlist,
            rule=rule_cp_mol_phase,
            doc="Phase isobaric heat capacity",
        )
        def rule_cp_mass_phase(b, p):
            return self.func_cp(cmp, delta[p], tau) * params.uc["kJ/kg/K to J/kg/K"]
        self.cp_mass_phase = Expression(
            phlist,
            rule=rule_cp_mass_phase,
            doc="Phase isobaric heat capacity",
        )
        # Phase constant volume heat capacity, cv
        def rule_cv_mol_phase(b, p):
            return self.func_cv(cmp, delta[p], tau) * params.uc["kJ/kg/K to J/mol/K"]
        self.cv_mol_phase = Expression(
            phlist,
            rule=rule_cv_mol_phase,
            doc="Phase isochoric heat capacity",
        )
        def rule_cv_mass_phase(b, p):
            return self.func_cv(cmp, delta[p], tau) * params.uc["kJ/kg/K to J/mass/K"]
        self.cv_mass_phase = Expression(
            phlist,
            rule=rule_cv_mass_phase,
            doc="Phase isochoric heat capacity",
        )
        # Phase speed of sound
        """
        def rule_speed_sound_phase(b, p):
            return self.func_w(cmp, delta[p], tau)

        self.speed_sound_phase = Expression(
            phlist,
            rule=rule_speed_sound_phase,
            doc="Phase speed of sound or saturated if phase doesn't exist",
        )
        """
        # Phase Mole density
        def rule_dens_mol_phase(b, p):
            return self.dens_mass_phase[p] / mw
        self.dens_mol_phase = Expression(
            phlist,
            rule=rule_dens_mol_phase,
            doc="Phase mole density",
        )
        # Component flow (for units that need it)
        def component_flow_mol(b, i):
            return self.flow_mol
        self.flow_mol_comp = Expression(
            component_list,
            rule=component_flow,
            doc="Total flow (both phases) of component",
        )
        def component_flow_mass(b, i):
            return self.flow_mass
        self.flow_mass_comp = Expression(
            component_list,
            rule=component_flow,
            doc="Total flow (both phases) of component",
        )

        #
        # Total (mixed phase) properties
        #

        # Enthalpy
        if not is_property_constructed(enth_mol):
            self.enth_mol = Expression(
                expr=sum(self.phase_frac[p] * self.enth_mol_phase[p] for p in phlist)
            )
        if not is_property_constructed(enth_mass):
            self.enth_mass = Expression(
                expr=sum(self.phase_frac[p] * self.enth_mass_phase[p] for p in phlist)
            )
        # Entropy
        if not is_property_constructed(entr_mol):
            self.entr_mol = Expression(
                expr=sum(self.phase_frac[p] * self.entr_mol_phase[p] for p in phlist)
            )
        if not is_property_constructed(entr_mass):
            self.entr_mass = Expression(
                expr=sum(self.phase_frac[p] * self.entr_mass_phase[p] for p in phlist)
            )
        # Internal Energy
        if not is_property_constructed(energy_internal_mol):
            self.energy_internal_mol = Expression(
                expr=sum(
                    self.phase_frac[p] * self.energy_internal_mol_phase[p] for p in phlist
                )
            )
        if not is_property_constructed(energy_internal_mass):
            self.energy_internal_mass = Expression(
                expr=sum(
                    self.phase_frac[p] * self.energy_internal_mass_phase[p] for p in phlist
                )
            )
        # cp
        self.cp_mol = Expression(
            expr=sum(self.phase_frac[p] * self.cp_mol_phase[p] for p in phlist)
        )
        self.cp_mass = Expression(
            expr=sum(self.phase_frac[p] * self.cp_mass_phase[p] for p in phlist)
        )
        # cv
        self.cv_mol = Expression(
            expr=sum(self.phase_frac[p] * self.cv_mol_phase[p] for p in phlist)
        )
        self.cv_mass = Expression(
            expr=sum(self.phase_frac[p] * self.cv_mass_phase[p] for p in phlist)
        )

        # mass density
        self.dens_mass = Expression(
            expr=1.0
            / sum(self.phase_frac[p] * 1.0 / self.dens_mass_phase[p] for p in phlist)
        )
        # mole density
        self.dens_mol = Expression(
            expr=1.0
            / sum(self.phase_frac[p] * 1.0 / self.dens_mol_phase[p] for p in phlist)
        )
        # heat capacity ratio
        self.heat_capacity_ratio = Expression(expr=self.cp_mol / self.cv_mol)
        # Flows
        self.flow_vol = Expression(
            expr=self.flow_mol / self.dens_mol,
            doc="Total liquid + vapor volumetric flow (m3/s)",
        )



        self.enth_mass = Expression(expr=self.enth_mol / mw, doc="Mass enthalpy (J/kg)")

        # Set the state vars dictionary
        if self.state_vars == StateVars.PH:
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "enth_mol": self.enth_mol,
                "pressure": self.pressure,
            }
        elif self.state_vars == StateVars.TPX and phase_set in (
            PhaseType.MIX,
            PhaseType.LG,
        ):
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure,
                "vapor_frac": self.vapor_frac,
            }
        elif self.state_vars == StateVars.TPX and phase_set in (
            PhaseType.G,
            PhaseType.L,
        ):
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure,
            }

        def rule_mole_frac_phase_comp(b, p, j):
            if p == "Mix":
                return 1.0
            else:
                return self.phase_frac[p]

        self.mole_frac_phase_comp = Expression(
            pub_phlist, component_list, rule=rule_mole_frac_phase_comp
        )

        # Define some expressions for the balance terms returned by functions
        # This is just to allow assigning scale factors to the expressions
        # returned
        #
        # Marterial flow term exprsssions
        def rule_material_flow_terms(b, p):
            if p == "Mix":
                return self.flow_mol
            else:
                return self.flow_mol * self.phase_frac[p]

        self.material_flow_terms = Expression(pub_phlist, rule=rule_material_flow_terms)

        # Enthaply flow term expressions
        def rule_enthalpy_flow_terms(b, p):
            if p == "Mix":
                return self.enth_mol * self.flow_mol
            else:
                return self.enth_mol_phase[p] * self.phase_frac[p] * self.flow_mol

        self.enthalpy_flow_terms = Expression(pub_phlist, rule=rule_enthalpy_flow_terms)

        # Energy density term expressions
        def rule_energy_density_terms(b, p):
            if p == "Mix":
                return self.dens_mol * self.energy_internal_mol
            else:
                return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]

        self.energy_density_terms = Expression(
            pub_phlist, rule=rule_energy_density_terms
        )

    def get_material_flow_terms(self, p, j):
        return self.material_flow_terms[p]

    def get_enthalpy_flow_terms(self, p):
        return self.enthalpy_flow_terms[p]

    def get_material_density_terms(self, p, j):
        if p == "Mix":
            return self.dens_mol
        else:
            return self.dens_mol_phase[p]

    def get_energy_density_terms(self, p):
        return self.energy_density_terms[p]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return self._state_vars_dict

    def define_display_vars(self):
        return {
            "Molar Flow (mol/s)": self.flow_mol,
            "Mass Flow (kg/s)": self.flow_mass,
            "T (K)": self.temperature,
            "P (Pa)": self.pressure,
            "Vapor Fraction": self.vapor_frac,
            "Molar Enthalpy (J/mol)": self.enth_mol_phase,
        }

    def extensive_state_vars(self):
        return self.extensive_set

    def intensive_state_vars(self):
        return self.intensive_set

    def model_check(self):
        pass

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        sf_flow = iscale.get_scaling_factor(self.flow_mol, default=1)
        sf_enth = iscale.get_scaling_factor(self.enth_mol, default=1)
        sf_inte = iscale.get_scaling_factor(self.energy_internal_mol, default=1)
        sf_dens = iscale.get_scaling_factor(self.dens_mol, default=1)
        sf_pres = iscale.get_scaling_factor(self.pressure, default=1)
        for v in self.material_flow_terms.values():
            iscale.set_scaling_factor(v, sf_flow)
        for v in self.enthalpy_flow_terms.values():
            iscale.set_scaling_factor(v, sf_enth * sf_flow)
        for k, v in self.energy_density_terms.items():
            if k == "Mix":
                iscale.set_scaling_factor(v, sf_inte * sf_dens)
            else:
                sf_inte_p = iscale.get_scaling_factor(
                    self.energy_internal_mol_phase[k], default=1
                )
                sf_dens_p = iscale.get_scaling_factor(self.dens_mol_phase[k], default=1)
                iscale.set_scaling_factor(v, sf_inte_p * sf_dens_p)
        try:
            iscale.set_scaling_factor(self.eq_sat, sf_pres / 1000.0)
        except AttributeError:
            pass  # may not have eq_sat, and that's ok
        try:
            iscale.set_scaling_factor(self.eq_complementarity, sf_pres / 10)
        except AttributeError:
            pass  # may not have eq_complementarity which is fine

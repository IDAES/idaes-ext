

__author__ = "John Eslick"

import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import (
    ControlVolume0DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
import idaes.logger as idaeslog
from idaes.core.util import (
    from_json,
    to_json,
    StoreSpec
)

from helmholtz_functions import (
    HelmholtzParameterBlock,
    HelmholtzParameterBlockData,
    HelmholtzThermoExpressions,
)
from helmholtz_state import (
    HelmholtzStateBlock,
)


_log = idaeslog.getModelLogger("RefrigerationUnits")

def is_helmholtz_parameter_block(val):
    if isinstance(val, HelmholtzParameterBlockData) or val == useDefault:
        return val
    else:
        raise ValueError("Property package must be a Helmholtz EOS")


def _add_property_pacakge_options(config):
    """Each unit model deffined here has the same property pacakage args so
    define a general function for add them."""

    config.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_helmholtz_parameter_block,
            description="Property package, must be Helmholtz EOS type",
            doc="Property package, must be Helmholtz EOS type"
                "**default** - useDefault."
                "**Valid values:** {"
                "**useDefault** - use default package from parent model or flowsheet"
                "**PropertyParameterObject** - a PropertyParameterBlock object.}",
        )
    )
    config.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="A ConfigBlock with arguments to be passed to a property block(s)"
                "and used when constructing these,"
                "**default** - None."
                "**Valid values:** {"
                "see property package for documentation.}",
        )
    )

def _add_property_pacakge_options_hx(config):
    """Each unit model deffined here has the same property pacakage args so
    define a general function for add them."""

    config.declare(
        "property_package_hot",
        ConfigValue(
            default=useDefault,
            domain=is_helmholtz_parameter_block,
            description="Property package, must be Helmholtz EOS type",
            doc="Property package, must be Helmholtz EOS type"
                "**default** - useDefault."
                "**Valid values:** {"
                "**useDefault** - use default package from parent model or flowsheet"
                "**PropertyParameterObject** - a PropertyParameterBlock object.}",
        )
    )
    config.declare(
        "property_package_args_hot",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="A ConfigBlock with arguments to be passed to a property block(s)"
                "and used when constructing these,"
                "**default** - None."
                "**Valid values:** {"
                "see property package for documentation.}",
        )
    )

    config.declare(
        "property_package_cold",
        ConfigValue(
            default=useDefault,
            domain=is_helmholtz_parameter_block,
            description="Property package, must be Helmholtz EOS type",
            doc="Property package, must be Helmholtz EOS type"
                "**default** - useDefault."
                "**Valid values:** {"
                "**useDefault** - use default package from parent model or flowsheet"
                "**PropertyParameterObject** - a PropertyParameterBlock object.}",
        )
    )
    config.declare(
        "property_package_args_cold",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="A ConfigBlock with arguments to be passed to a property block(s)"
                "and used when constructing these,"
                "**default** - None."
                "**Valid values:** {"
                "see property package for documentation.}",
        )
    )


@declare_process_block_class("Compressor")
class CompressorData(UnitModelBlockData):
    # Add unit configuration options
    CONFIG = UnitModelBlockData.CONFIG()
    _add_property_pacakge_options(CONFIG)

    def build(self):
        """
        Add equations to the unit model. This is called by a default block
        construnction rule when the unit model is created.
        """
        super().build() # Basic unit model build/read config
        config = self.config # shorter config pointer

        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=config.property_package,
            property_package_args=config.property_package_args
        )
        # Using the mixed state version of the properties so we want the
        # unit models to see this as one mixed phase, and leave equilibrium
        # to the property package.
        self.control_volume.add_state_blocks(has_phase_equilibrium=False)
        self.add_inlet_port()
        self.add_outlet_port()

        eff = self.efficiency = pyo.Var(initialize=0.9)
        self.efficiency.fix()
        self.efficiency.latex_symbol = "\eta"

        # The locations of properies can be a little cumbersome for writing
        # clear readale constraints, so make some short pointers to help,
        # for the compressor we'll use three states: 0, inlet, 1 outlet, and s,
        # isentropic
        prop_0 = self.control_volume.properties_in[0]
        prop_1 = self.control_volume.properties_out[0]
        s_0 = prop_0.entr_mass
        s_0.latex_symbol = "s_0"
        p_1 = prop_1.pressure
        p_1.latex_symbol = "p_1"
        h_0 = prop_0.enth_mass
        h_0.latex_symbol = "h_0"
        h_1 = prop_1.enth_mass
        h_1.latex_symbol = "h_1"
        F_0 = prop_0.flow_mass
        F_0.latex_symbol = "F_0"
        F_1 = prop_1.flow_mass
        F_1.latex_symbol = "F_1"

        # Thermo expression writer, this will create thermo expression that look
        # like thermodynaic function calls when writing the constraints.
        te = HelmholtzThermoExpressions(
            blk=self,
            parameters=config.property_package
        )

        # Expression for isentropic enthalpy
        h_s = self.h_s = pyo.Expression(expr=te.h(s=s_0, p=p_1))
        h_s.latex_nice_expr = f"h(c={self.config.property_package.pure_component}, s=s_0, p=p_1)"
        h_s.latex_symbol = "h_s"
        # Calculate work
        self.specific_work =  pyo.Expression(expr=(h_s - h_0)/eff, doc="Work per mass")
        self.specific_work.latex_symbol = "w"
        self.work = pyo.Expression(expr=self.specific_work*F_0)
        self.work.latex_symbol = "\dot{w}"
        # Energy balance
        self.eq_enth_out = pyo.Constraint(expr=h_1 == self.specific_work + h_0)
        # Mass balance
        self.eq_flow_out = pyo.Constraint(expr=F_1 == F_0)

    def initialize_build(
        self,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6},
    ):
        """
        Simple initialization, assume that a good guess for the inlet,
        efficiency, and outlet pressure are provided and solve.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Set solver options
        solver = pyo.SolverFactory(solver)
        solver.options = optarg
        prop_0 = self.control_volume.properties_in[0]
        prop_1 = self.control_volume.properties_out[0]

        # Set flows
        if prop_0.flow_mass.fixed and prop_1.flow_mass.fixed:
            raise RuntimeError("Compressor inlet and outlet flow are fixed")
        elif prop_1.flow_mass.fixed:
            prop_0.flow_mass.value = pyo.value(prop_1.flow_mass)
        else:
            prop_1.flow_mass.value = pyo.value(prop_0.flow_mass)

        # Set enthalpy
        if prop_0.enth_mass.fixed and prop_1.enth_mass.fixed:
            raise RuntimeError("Compressor inlet and outlet flow are fixed")
        elif prop_1.enth_mass.fixed:
            prop_0.enth_mass.value = pyo.value(prop_1.enth_mass)
        else:
            prop_1.enth_mass.value = pyo.value(prop_0.enth_mass)

        # Save model state, so fixed/free, active/inative, and fixed values
        # can be restored, ensuring the propblem specs don't change due to
        # initialization
        #
        # sp is what to save/load from json, the only_fixed option only loads
        # originally fixed values, so unfixed values can change in initialization,
        # this also saves which variables are fixed and which constraints are
        # active.
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # Fix inlet and efficiency and solve
        self.inlet.fix()
        self.outlet.unfix()
        self.efficiency.fix()
        self.outlet.pressure.fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)

        # restore specs
        from_json(self, sd=istate, wts=sp)
        init_log.info_high("Compressor Initialization Complete")


@declare_process_block_class("Valve")
class ValveData(UnitModelBlockData):
    # Add unit configuration options
    CONFIG = UnitModelBlockData.CONFIG()
    _add_property_pacakge_options(CONFIG)

    def build(self):
        """
        Add equations to the unit model. This is called by a default block
        construnction rule when the unit model is created.
        """
        super().build() # Basic unit model build/read config
        config = self.config # shorter config pointer

        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=config.property_package,
            property_package_args=config.property_package_args
        )
        # Using the mixed state version of the properties so we want the
        # unit models to see this as one mixed phase, and leave equilibrium
        # to the property package.
        self.control_volume.add_state_blocks(has_phase_equilibrium=False)
        self.add_inlet_port()
        self.add_outlet_port()

        # The locations of properies can be a little cumbersome for writing
        # clear readale constraints, so make some short pointers to help,
        # for the compressor we'll use three states: 0, inlet, 1 outlet, and s,
        # isentropic
        prop_0 = self.control_volume.properties_in[0]
        prop_1 = self.control_volume.properties_out[0]
        p_1 = prop_1.pressure
        h_0 = prop_0.enth_mass
        h_1 = prop_1.enth_mass
        F_0 = prop_0.flow_mass
        F_1 = prop_1.flow_mass

        # Thermo expression writer, this will create thermo expression that look
        # like thermodynaic function calls when writing the constraints.
        te = HelmholtzThermoExpressions(
            blk=self,
            parameters=config.property_package
        )
        # Energy
        self.eq_enth_out = pyo.Constraint(expr=h_1 == h_0)
        # Mass balance
        self.eq_flow_out = pyo.Constraint(expr=F_0 == F_1)

    def initialize_build(
        self,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6},
    ):
        """
        Simple initialization, assume that a good guess for the inlet and
        outlet pressure are provided and solve.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Set solver options
        solver = pyo.SolverFactory(solver)
        solver.options = optarg
        prop_0 = self.control_volume.properties_in[0]
        prop_1 = self.control_volume.properties_out[0]

        # Set flows
        if prop_0.flow_mass.fixed and prop_1.flow_mass.fixed:
            raise RuntimeError("Compressor inlet and outlet flow are fixed")
        elif prop_1.flow_mass.fixed:
            prop_0.flow_mass.value = pyo.value(prop_1.flow_mass)
        else:
            prop_1.flow_mass.value = pyo.value(prop_0.flow_mass)

        # Set enthalpy
        if prop_0.enth_mass.fixed and prop_1.enth_mass.fixed:
            raise RuntimeError("Compressor inlet and outlet flow are fixed")
        elif prop_1.enth_mass.fixed:
            prop_0.enth_mass.value = pyo.value(prop_1.enth_mass)
        else:
            prop_1.enth_mass.value = pyo.value(prop_0.enth_mass)

        # Save model state, so fixed/free, active/inative, and fixed values
        # can be restored, ensuring the propblem specs don't change due to
        # initialization
        #
        # sp is what to save/load from json, the only_fixed option only loads
        # originally fixed values, so unfixed values can change in initialization,
        # this also saves which variables are fixed and which constraints are
        # active.
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # Fix inlet and solve
        self.inlet.fix()
        self.outlet.unfix()
        self.outlet.pressure.fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)

        # restore specs
        from_json(self, sd=istate, wts=sp)
        init_log.info_high("Valve Initialization Complete")


@declare_process_block_class("HeatExchanger")
class HeatExchangerData(UnitModelBlockData):
    # Add unit configuration options
    CONFIG = UnitModelBlockData.CONFIG()
    _add_property_pacakge_options_hx(CONFIG)

    def build(self):
        """
        Add equations to the unit model. This is called by a default block
        construnction rule when the unit model is created.
        """
        super().build() # Basic unit model build/read config
        config = self.config # shorter config pointer

        self.hot_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=config.property_package_hot,
            property_package_args=config.property_package_args_hot
        )

        self.cold_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=config.property_package_cold,
            property_package_args=config.property_package_args_cold
        )

        self.hot_side.add_state_blocks(has_phase_equilibrium=False)
        self.cold_side.add_state_blocks(has_phase_equilibrium=False)

        self.add_inlet_port(name="hot_inlet", block=self.hot_side)
        self.add_inlet_port(name="cold_inlet", block=self.cold_side)
        self.add_outlet_port(name="hot_outlet", block=self.hot_side)
        self.add_outlet_port(name="cold_outlet", block=self.cold_side)

        prop_0h = self.hot_side.properties_in[0]
        prop_1h = self.hot_side.properties_out[0]
        h_0h = prop_0h.enth_mass
        h_0h.latex_symbol = "h_{\\textrm{hot},0}"
        h_1h = prop_1h.enth_mass
        h_1h.latex_symbol = "h_{\\textrm{hot},1}"
        F_0h = prop_0h.flow_mass
        F_0h.latex_symbol = "F_{\\textrm{hot},0}"
        F_1h = prop_1h.flow_mass
        F_1h.latex_symbol = "F_{\\textrm{hot},1}"
        T_0h = prop_0h.temperature
        T_0h.latex_symbol = "T_{\\textrm{hot},0}"
        T_1h = prop_1h.temperature
        T_1h.latex_symbol = "T_{\\textrm{hot},1}"
        P_0h = prop_0h.pressure
        P_0h.latex_symbol = "P_{\\textrm{hot},0}"
        P_1h = prop_1h.pressure
        P_1h.latex_symbol = "P_{\\textrm{hot},1}"

        prop_0c = self.cold_side.properties_in[0]
        prop_1c = self.cold_side.properties_out[0]
        h_0c = prop_0c.enth_mass
        h_0c.latex_symbol = "h_{\\textrm{cold},0}"
        h_1c = prop_1c.enth_mass
        h_1c.latex_symbol = "h_{\\textrm{cold},1}"
        F_0c = prop_0c.flow_mass
        F_0c.latex_symbol = "F_{\\textrm{cold},0}"
        F_1c = prop_1c.flow_mass
        F_1c.latex_symbol = "F_{\\textrm{cold},1}"
        T_0c = prop_0c.temperature
        T_0c.latex_symbol = "T_{\\textrm{cold},0}"
        T_1c = prop_1c.temperature
        T_1c.latex_symbol = "T_{\\textrm{cold},1}"
        P_0c = prop_0c.pressure
        P_0c.latex_symbol = "P_{\\textrm{cold},0}"
        P_1c = prop_1c.pressure
        P_1c.latex_symbol = "P_{\\textrm{cold},1}"

        self.U = pyo.Var(initialize=100, units=pyo.units.W/pyo.units.m**2/pyo.units.K)
        self.U.latex_symbol = "U"
        self.A = pyo.Var(initialize=600, units=pyo.units.m**2)
        self.A.latex_symbol = "A"

        self.U.fix()
        self.A.fix()

        # I know it's a condenser and this isn't quite right, but well assume
        # countercurrent for now
        dTA = T_1h - T_0c
        dTB = T_0h - T_1c
        self.lmtd = pyo.Expression(expr=(dTB - dTA)/(pyo.log(dTB/dTA)))
        self.lmtd.latex_symbol = "\Delta T"

        # Energy Balances
        self.Q = pyo.Expression(expr=self.lmtd*self.U*self.A)
        self.Q.latex_symbol = "Q"
        self.hot_energy_eq = pyo.Constraint(expr=h_1h == h_0h - self.Q/F_0h)
        self.cold_energy_eq = pyo.Constraint(expr=h_1c == h_0c + self.Q/F_0c)

        # Mass Balances
        self.hot_mass_eq = pyo.Constraint(expr=F_1h == F_0h)
        self.cold_mass_eq = pyo.Constraint(expr=F_1c == F_0c)

        # Pressure Drop
        self.hot_pressure_eq = pyo.Constraint(expr=P_1h == P_0h)
        self.cold_pressure_eq = pyo.Constraint(expr=P_1c == P_0c)

    def initialize_build(
        self,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6},
    ):
        """
        Simple initialization, assume that a good guess for the inlet and
        outlet pressure are provided and solve.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Set solver options
        solver = pyo.SolverFactory(solver)
        solver.options = optarg
        prop_0h = self.hot_side.properties_in[0]
        prop_1h = self.hot_side.properties_out[0]
        prop_0c = self.cold_side.properties_in[0]
        prop_1c = self.cold_side.properties_out[0]

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        prop_1h.flow_mass.value = pyo.value(prop_0h.flow_mass)
        prop_1c.flow_mass.value = pyo.value(prop_0c.flow_mass)
        prop_1h.enth_mass.value = pyo.value(prop_0h.enth_mass - 1)
        prop_1c.enth_mass.value = pyo.value(prop_0c.enth_mass + 1)
        prop_1h.pressure.value = pyo.value(prop_0h.pressure)
        prop_1c.pressure.value = pyo.value(prop_0c.pressure)

        self.cold_inlet.fix()
        self.hot_inlet.fix()
        self.cold_outlet.unfix()
        self.hot_outlet.unfix()
        self.U.fix()
        self.A.fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(self, tee=slc.tee)

        from_json(self, sd=istate, wts=sp)
        init_log.info_high(f"HeatExchanger '{self.name}' Initialization Complete")

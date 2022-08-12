__author__ = "John Eslick"

import os

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir

from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
import idaes.core.util.tables as tables
from idaes.core.util.tags import svg_tag
import idaes.core.util as iutil

from ref_units import Compressor, HeatExchanger, Valve
from helmholtz_functions import HelmholtzParameterBlock, AmountBasis
from idaes.core.util.initialization import propagate_state

_log = idaeslog.getModelLogger("RefrigerationFlowsheet")


def main():
    """ Create the refrigeration cycle flowsheet """
    solver = pyo.SolverFactory("ipopt")
    solver.options = {"tol": 1e-6, "halt_on_ampl_error": "no", "max_iter": 50}

    # Create the Pyomo model with a flowsheet and property parameter block
    m = pyo.ConcreteModel(name="Example Refrigeration Cycle")
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_ref = HelmholtzParameterBlock(pure_component="r1234ze", amount_basis=AmountBasis.MASS)
    m.fs.prop_h2o = HelmholtzParameterBlock(pure_component="h2o", amount_basis=AmountBasis.MASS)
    m.fs.compressor = Compressor(property_package=m.fs.prop_ref)

    m.fs.compressor.inlet.flow_mass.fix(1*pyo.units.kg/pyo.units.s)
    m.fs.compressor.inlet.enth_mass.fix(m.fs.prop_ref.htpx(p=300*pyo.units.kPa, x=1.0))
    m.fs.compressor.inlet.pressure.fix(300*pyo.units.kPa)
    m.fs.compressor.efficiency.fix(0.80)
    m.fs.compressor.outlet.pressure.fix(1400*pyo.units.kPa)

    m.fs.compressor.initialize()

    m.fs.condenser = HeatExchanger(
        property_package_hot=m.fs.prop_ref,
        property_package_cold=m.fs.prop_h2o,
    )

    m.fs.condenser.cold_inlet.flow_mass.fix(5*pyo.units.kg/pyo.units.s)
    m.fs.condenser.cold_inlet.pressure.fix(200*pyo.units.kPa)
    m.fs.condenser.cold_inlet.enth_mass.fix(m.fs.prop_h2o.htpx(p=200*pyo.units.kPa, T=320.0*pyo.units.K))
    m.fs.condenser.A.fix(100*pyo.units.m**2)
    m.fs.sB = Arc(
        source=m.fs.compressor.outlet,
        destination=m.fs.condenser.hot_inlet
    )
    propagate_state(arc=m.fs.sB)
    m.fs.condenser.initialize()

    m.fs.valve = Valve(property_package=m.fs.prop_ref)
    m.fs.valve.outlet.pressure.fix(300*pyo.units.kPa)
    m.fs.sD = Arc(
        source=m.fs.condenser.hot_outlet,
        destination=m.fs.valve.inlet
    )
    propagate_state(arc=m.fs.sD)
    m.fs.valve.initialize()

    m.fs.evaporator = HeatExchanger(
        property_package_hot=m.fs.prop_h2o,
        property_package_cold=m.fs.prop_ref,
    )

    m.fs.evaporator.hot_inlet.flow_mass.fix(1*pyo.units.kg/pyo.units.s)
    m.fs.evaporator.hot_inlet.pressure.fix(200*pyo.units.kPa)
    m.fs.evaporator.hot_inlet.enth_mass.fix(m.fs.prop_h2o.htpx(p=200*pyo.units.kPa, T=320.0*pyo.units.K))
    m.fs.evaporator.A.fix(50*pyo.units.m**2)

    m.fs.sE = Arc(
        source=m.fs.valve.outlet,
        destination=m.fs.evaporator.cold_inlet
    )
    propagate_state(arc=m.fs.sE)
    m.fs.evaporator.initialize()


    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    m.fs.compressor.inlet.flow_mass.unfix()
    m.fs.cond_eq = pyo.Constraint(expr=m.fs.condenser.hot_side.properties_out[0].enth_mass/1000 == m.fs.condenser.hot_side.properties_out[0].enth_mass_sat_phase["Liq"]/1000)
    solver.solve(m, tee=True)

    m.fs.sA = Arc(
        source=m.fs.evaporator.cold_outlet,
        destination=m.fs.compressor.inlet
    )
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    m.fs.compressor.inlet.flow_mass.unfix()
    m.fs.compressor.inlet.enth_mass.unfix()
    m.fs.compressor.inlet.pressure.unfix()

    m.fs.evaporator.hot_inlet.flow_mass.unfix()
    m.fs.evap_eq = pyo.Constraint(expr=m.fs.evaporator.cold_side.properties_out[0].enth_mass/1000 == m.fs.evaporator.cold_side.properties_out[0].enth_mass_sat_phase["Vap"]/1000)


    # Flow constraint from arc is redundant, since the mass balances
    # provide the same
    m.fs.sA.expanded_block.flow_mass_equality.deactivate()

    solver.solve(m, tee=True)
    _add_tags(m.fs)
    return m, solver

def _add_tags(self):
    tag_group = iutil.ModelTagGroup()
    self.tags_streams = tag_group
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            self,
            additional={
                "w1": self.evaporator.hot_inlet,
                "w2": self.evaporator.hot_outlet,
                "w3": self.condenser.cold_inlet,
                "w4": self.condenser.cold_outlet,
            },
            descend_into=False,
        )
    )
    for i, s in stream_states.items():
        tag_group[f"{i}_F"] = iutil.ModelTag(
            doc=f"{i}: mass flow",
            expr=s.flow_mass,
            format_string="{:.3f}",
            display_units=pyo.units.kg / pyo.units.s,
        )
        tag_group[f"{i}_Fvol"] = iutil.ModelTag(
            doc=f"{i}: volumetric flow",
            expr=s.flow_vol,
            format_string="{:.3f}",
            display_units=pyo.units.m**3 / pyo.units.s,
        )
        tag_group[f"{i}_p"] = iutil.ModelTag(
            doc=f"{i}: pressure",
            expr=s.pressure,
            format_string="{:.3f}",
            display_units=pyo.units.kPa,
        )
        tag_group[f"{i}_T"] = iutil.ModelTag(
            doc=f"{i}: temperature",
            expr=s.temperature,
            format_string="{:.2f}",
            display_units=pyo.units.K,
        )
        tag_group[f"{i}_x"] = iutil.ModelTag(
            doc=f"{i}: vapor fraction",
            expr=s.phase_frac["Vap"],
            format_string="{:.3f}",
            display_units=None,
        )
        tag_group[f"{i}_h"] = iutil.ModelTag(
            doc=f"{i}: specific enthalpy",
            expr=s.enth_mass,
            format_string="{:.3f}",
            display_units=pyo.units.kJ / pyo.units.kg,
        )
        tag_group[f"{i}_s"] = iutil.ModelTag(
            doc=f"{i}: specific entropy",
            expr=s.entr_mass,
            format_string="{:.3f}",
            display_units=pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    tag_group = iutil.ModelTagGroup()
    self.tags_unit_output = tag_group
    tag_group["cmp_eff"] = iutil.ModelTag(
        doc="Compressor isentropic efficiency",
        expr=self.compressor.efficiency,
        format_string="{:.3f}",
        display_units=None
    )
    tag_group["cmp_work"] = iutil.ModelTag(
        doc="Compressor work",
        expr=self.compressor.work,
        format_string="{:.1f}",
        display_units=pyo.units.kW
    )
    tag_group["evap_A"] = iutil.ModelTag(
        doc="Evaporator area",
        expr=self.evaporator.A,
        format_string="{:.1f}",
        display_units=pyo.units.m**2
    )
    tag_group["evap_Q"] = iutil.ModelTag(
        doc="Evaporator heat transfer",
        expr=self.evaporator.Q,
        format_string="{:.1f}",
        display_units=pyo.units.kW
    )
    tag_group["evap_U"] = iutil.ModelTag(
        doc="Evaporator heat transfer coefficient",
        expr=self.evaporator.U,
        format_string="{:.1f}",
        display_units=pyo.units.W/pyo.units.m**2/pyo.units.K
    )
    tag_group["evap_lmtd"] = iutil.ModelTag(
        doc="Evaporator LMTD",
        expr=self.evaporator.lmtd,
        format_string="{:.1f}",
        display_units=pyo.units.K
    )
    tag_group["cond_A"] = iutil.ModelTag(
        doc="Condenser area",
        expr=self.condenser.A,
        format_string="{:.1f}",
        display_units=pyo.units.m**2
    )
    tag_group["cond_Q"] = iutil.ModelTag(
        doc="Condenser heat transfer",
        expr=self.condenser.Q,
        format_string="{:.1f}",
        display_units=pyo.units.kW
    )
    tag_group["cond_U"] = iutil.ModelTag(
        doc="Condenser heat transfer coefficient",
        expr=self.condenser.U,
        format_string="{:.1f}",
        display_units=pyo.units.W/pyo.units.m**2/pyo.units.K
    )
    tag_group["cond_lmtd"] = iutil.ModelTag(
        doc="Condenser LMTD",
        expr=self.condenser.lmtd,
        format_string="{:.1f}",
        display_units=pyo.units.K
    )

def write_pfd(self, fname=None):
    """Add model results to the flowsheet template.  If fname is specified,
    this saves the resulting svg to a file.  If fname is not specified, it
    returns the svg string.
    Args:
        fname: Name of file to save svg.  If None, return the svg string
    Returns: (None or Str)
    """
    infilename = os.path.join(this_file_dir(), "ref_diagram.svg")
    with open(infilename, "r") as f:
        s = svg_tag(svg=f, tag_group=self.tags_streams)
        s = svg_tag(svg=s, tag_group=self.tags_unit_output, outfile=fname)
    if fname is None:
        return s

if __name__ == "__main__":
    m, solver = main()

    m.fs.compressor.control_volume.properties_in[0].temperature.display()
    m.fs.compressor.control_volume.properties_out[0].temperature.display()
    m.fs.compressor.control_volume.properties_out[0].temperature_sat.display()
    m.fs.compressor.control_volume.properties_out[0].vapor_frac.display()
    m.fs.compressor.specific_work.display()

    pA = m.fs.compressor.control_volume.properties_in[0]
    pB = m.fs.compressor.control_volume.properties_out[0]
    pD = m.fs.condenser.hot_side.properties_out[0]
    pE = m.fs.valve.control_volume.properties_out[0]
    pA2 = m.fs.evaporator.cold_side.properties_out[0]
    pCW = m.fs.evaporator.hot_side.properties_out[0]



    ph_plot = {
        "A": (pyo.value(pA.enth_mass*m.fs.prop_ref.uc["J/kg to kJ/kg"]), pyo.value(pA.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "B": (pyo.value(pB.enth_mass*m.fs.prop_ref.uc["J/kg to kJ/kg"]), pyo.value(pB.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "C": (pyo.value(pB.enth_mass_sat_phase["Vap"]*m.fs.prop_ref.uc["J/kg to kJ/kg"]), pyo.value(pB.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "D": (pyo.value(pD.enth_mass*m.fs.prop_ref.uc["J/kg to kJ/kg"]), pyo.value(pD.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "E": (pyo.value(pD.enth_mass*m.fs.prop_ref.uc["J/kg to kJ/kg"]), pyo.value(pE.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
    }

    m.fs.prop_ref.ph_diagram(points=ph_plot)

    pA.flow_mass.display()
    pE.temperature.display()

    pCW.temperature.display()

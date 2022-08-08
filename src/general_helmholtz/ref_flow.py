__author__ = "John Eslick"

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog

from ref_units import Compressor, HeatExchanger, Valve
from helmholtz_functions import HelmholtzParameterBlock
from idaes.core.util.initialization import propagate_state

_log = idaeslog.getModelLogger("RefrigerationFlowsheet")


def main():
    """ Create the refrigeration cycle flowsheet """
    solver = pyo.SolverFactory("ipopt")
    solver.options = {"tol": 1e-6, "halt_on_ampl_error": "no", "max_iter": 50}

    # Create the Pyomo model with a flowsheet and property parameter block
    m = pyo.ConcreteModel(name="Example Refrigeration Cycle")
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_ref = HelmholtzParameterBlock(pure_component="r1234ze")
    m.fs.prop_h2o = HelmholtzParameterBlock(pure_component="h2o")
    m.fs.compressor = Compressor(property_package=m.fs.prop_ref)

    m.fs.compressor.inlet.flow_mol.fix(10*pyo.units.mol/pyo.units.s)
    m.fs.compressor.inlet.enth_mol.fix(m.fs.prop_ref.htpx(p=300*pyo.units.kPa, x=1.0))
    m.fs.compressor.inlet.pressure.fix(300*pyo.units.kPa)
    m.fs.compressor.efficiency.fix(0.80)
    m.fs.compressor.outlet.pressure.fix(1400*pyo.units.kPa)

    m.fs.compressor.initialize()

    m.fs.condenser = HeatExchanger(
        property_package_hot=m.fs.prop_ref,
        property_package_cold=m.fs.prop_h2o,
    )

    m.fs.condenser.cold_inlet.flow_mol.fix(100*pyo.units.mol/pyo.units.s)
    m.fs.condenser.cold_inlet.pressure.fix(200*pyo.units.kPa)
    m.fs.condenser.cold_inlet.enth_mol.fix(m.fs.prop_h2o.htpx(p=200*pyo.units.kPa, T=320.0*pyo.units.K))

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

    m.fs.evaporator.hot_inlet.flow_mol.fix(65*pyo.units.mol/pyo.units.s)
    m.fs.evaporator.hot_inlet.pressure.fix(200*pyo.units.kPa)
    m.fs.evaporator.hot_inlet.enth_mol.fix(m.fs.prop_h2o.htpx(p=200*pyo.units.kPa, T=320.0*pyo.units.K))
    m.fs.evaporator.A.fix(45)

    m.fs.sE = Arc(
        source=m.fs.valve.outlet,
        destination=m.fs.evaporator.cold_inlet
    )
    propagate_state(arc=m.fs.sE)
    m.fs.evaporator.initialize()


    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    m.fs.compressor.inlet.flow_mol.unfix()
    m.fs.cond_eq = pyo.Constraint(expr=m.fs.condenser.hot_side.properties_out[0].enth_mol == m.fs.condenser.hot_side.properties_out[0].enth_mol_sat_phase["Liq"])
    solver.solve(m, tee=True)

    m.fs.sA = Arc(
        source=m.fs.evaporator.cold_outlet,
        destination=m.fs.compressor.inlet
    )
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    m.fs.compressor.inlet.flow_mol.unfix()
    m.fs.compressor.inlet.enth_mol.unfix()
    m.fs.compressor.inlet.pressure.unfix()

    # Flow constraint from arc is redundant, since the mass balances
    # provide the same
    m.fs.sA.expanded_block.flow_mol_equality.deactivate()

    solver.solve(m, tee=True)

    return m, solver

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
        "A": (pyo.value(pA.enth_mol*m.fs.prop_ref.uc["J/mol to kJ/kg"]), pyo.value(pA.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "B": (pyo.value(pB.enth_mol*m.fs.prop_ref.uc["J/mol to kJ/kg"]), pyo.value(pB.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "C": (pyo.value(pB.enth_mol_sat_phase["Vap"]*m.fs.prop_ref.uc["J/mol to kJ/kg"]), pyo.value(pB.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "D": (pyo.value(pD.enth_mol*m.fs.prop_ref.uc["J/mol to kJ/kg"]), pyo.value(pD.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
        "E": (pyo.value(pD.enth_mol*m.fs.prop_ref.uc["J/mol to kJ/kg"]), pyo.value(pE.pressure*m.fs.prop_ref.uc["Pa to kPa"])),
    }

    m.fs.prop_ref.ph_diagram(points=ph_plot)

    pA.flow_mol.display()
    pE.temperature.display()

    pCW.temperature.display()

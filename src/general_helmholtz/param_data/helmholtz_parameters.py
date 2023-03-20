import json

import pyomo.environ as pyo
import logging

_log = logging.getLogger("idaes.helmholtz_parameters")

class WriteParameters(object):

    variables = {
        "delta": 0,
        "tau": 1,
        "p": 2,
        "T": 3,
    }

    expressions = { # these are expressions in the thermo model. Other models 
        # are separate 
        # phi is dimensionless Helmholtz free energy
        "phii": 0,  # ideal part of phi(delta, tau)
        "phii_d": 1,  # ideal part of phi partial wrt delta
        "phii_dd": 2,  # ideal part of phi partial wrt delta and delta
        "phii_t": 3,  # ideal part of phi partial wrt tau
        "phii_tt": 4,  # ideal part of phi partial tau delta and tau
        "phii_dt": 5,  # ideal part of phi partial wrt delta and tau
        "phir": 6,  # residual part of phi(delta, tau)
        "phir_d": 7,  # residual part of phi partial wrt delta
        "phir_dd": 8,  # residual part of phi partial wrt delta and delta
        "phir_t": 9,  # residual part of phi partial wrt tau
        "phir_tt": 10,  # residual part of phi partial wrt tau and tau
        "phir_dt": 11,  # residual part of phi partial wrt delta and tau
        # Functions for initial guess on saturated liquid and vapor reduced 
        #   density. These aid the phase equilibrium calculations, but their
        #   results are not directly used.
        "delta_v_sat_approx": 12, # approximate delta_sat_vapor(tau)
        "delta_l_sat_approx": 13, # approximate delta_sat_liquid(tau)
    }

    optional_expressions = {
        "thermal_conductivity": "tcx",
        "surface_tension": "st",
        "viscosity": "visc",
    }

    def __init__(
        self,
        comp,
        R,  # specific ideal gas constant [kJ/kg/K]
        MW,  # molecular weight [g/mol]
        T_star,  # for calculating tau = T_star/T [K]
        rho_star,  # for calculating delta = rho_star/rho [kg/m3]
        Tc,  # critical temperature [K]
        rhoc,  # critical density [kg/m3]
        Pc,  # critical pressure [kPa]
        Tt,  # triple point temperature [K]
        Pt,  # triple point pressure [K]
        rhot_l,  # liquid triple point density [kg/m3]
        rhot_v,  # vapor triple point density [kg/m3]
        P_min,  # minimum pressure [kPa]
        P_max,  # maximum pressure [kPa]
        rho_max,  # maximum density [kg/m3]
        T_min,  # minimum temperature [kPa]
        T_max,  # maximum temperature [kPa]
    ):
        self.comp = comp
        self.R = R
        self.MW = MW
        self.T_star = T_star
        self.rho_star = rho_star
        self.Tc = Tc
        self.rhoc = rhoc
        self.Pc = Pc
        self.Tt = Tt
        self.Pt = Pt
        self.rhot_l = rhot_l
        self.rhot_v = rhot_v
        self.P_min = P_min
        self.P_max = P_max
        self.rho_max = rho_max
        self.T_min = T_min
        self.T_max = T_max
    
        # Main thermo model
        self.model = self.make_model("delta", "tau")
        # Transport models (Thermal Conductivity, Viscosity, and Surface Tension)
        # keep these separated from other models, because they may need to call
        # other models as part of the calculation and don't want them interfering
        # with each other. 
        self.model_tcx = self.make_model("delta", "tau")
        self.model_visc = self.make_model("delta", "tau")
        self.model_st = self.make_model("tau") 
        self.has_expression = []

    def make_model(self, *args):
        m = pyo.ConcreteModel()
        for a in args:
            setattr(m, a, pyo.Var())
        m.R = pyo.Param(initialize=self.R)
        m.MW = pyo.Param(initialize=self.MW)
        m.T_star = pyo.Param(initialize=self.T_star)
        m.rho_star = pyo.Param(initialize=self.rho_star)
        m.Tc = pyo.Param(initialize=self.Tc)
        m.rhoc = pyo.Param(initialize=self.rhoc)
        m.Pc = pyo.Param(initialize=self.Pc)
        m.Tt = pyo.Param(initialize=self.Tt)
        m.Pt = pyo.Param(initialize=self.Pt)
        m.rhot_l = pyo.Param(initialize=self.rhot_l)
        m.rhot_v = pyo.Param(initialize=self.rhot_v)
        m.P_min = pyo.Param(initialize=self.P_min)
        m.P_max = pyo.Param(initialize=self.P_max)
        m.rho_max = pyo.Param(initialize=self.rho_max)
        m.T_min = pyo.Param(initialize=self.T_min)
        m.T_max = pyo.Param(initialize=self.T_max)
        return m

    def add(self, expressions):
        for name, expr in expressions.items():
            if name in self.expressions: # check if in thermo model
                m = self.model
            elif name == "thermal_conductivity":
                m = self.model_tcx
            elif name == "surface_tension":
                m = self.model_st
            elif name == "viscosity":
                m = self.model_visc
            else:
                raise RuntimeError(f"Unknown expression {name}")
            if callable(expr):
                setattr(m, name, pyo.Objective(rule=expr))
            else:
                setattr(m, name, pyo.Objective(expr=expr))
            self.has_expression.append(name)
            
    def write_model(self, model, model_name, expressions=None):
        nl_file, smap_id = model.write(f"{self.comp}_expressions_{model_name}.nl")
        smap = model.solutions.symbol_map[smap_id]
        var_map = [1000]*4
        for s, c in smap.bySymbol.items():
            if s.startswith("v"):
                j = int(s[1:])
                var_map[j] = self.variables[c().name]
        if expressions is not None:
            expr_map = [0] * len(expressions)
            for s, c in smap.bySymbol.items():
                if s.startswith("o"):
                    i = expressions[c().name]
                    j = int(s[1:])
                    expr_map[i] = j
            return nl_file, expr_map, var_map
        return nl_file, var_map

    def write(self):
        for name in self.expressions:
            if name not in self.has_expression:
                raise RuntimeError(f"Required expression {name} not provided.")
        nl_file, expr_map, var_map = self.write_model(self.model, "eos", self.expressions)
        param_dict = {
            "nl_file": nl_file,
            "expr_map": expr_map,
            "var_map": var_map,
            "param": {
                "R": self.R,
                "MW": self.MW,
                "T_star": self.T_star,
                "rho_star": self.rho_star,
                "Tc": self.Tc,
                "rhoc": self.rhoc,
                "Pc": self.Pc,
                "Tt": self.Tt,
                "Pt": self.Pt,
                "rhot_l": self.rhot_l,
                "rhot_v": self.rhot_v,
                "P_min": self.P_min,
                "P_max": self.P_max,
                "rho_max": self.rho_max,
                "T_min": self.T_min,
                "T_max": self.T_max,
            },
        }

        # Add optional models
        for name, short_name in self.optional_expressions.items():
            if name in self.has_expression:
                model = getattr(self, f"model_{short_name}")
                nl_file, var_map = self.write_model(model, short_name)
                param_dict[f"nl_file_{short_name}"] = nl_file
                param_dict[f"var_map_{short_name}"] = var_map
                param_dict[f"have_{short_name}"] = True
            else:
                _log.warning(f"Missing optional expression {name}")
                param_dict[f"have_{short_name}"] = False

        with open(f"{self.comp}_parameters.json", "w") as f:
            json.dump(param_dict, f, indent=4)

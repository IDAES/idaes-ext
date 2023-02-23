import pyomo.environ as pyo
import json


class WriteParameters(object):
    is_omitted_index = 1000

    variables = {
        "delta": 0,
        "tau": 1,
        "p": 2,
    }

    expressions = {
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
        "delta_v_sat_approx": 12, # approximate delta_sat_vapor(tau)
        "delta_l_sat_approx": 13, # approximate delta_sat_liquid(tau)
        "viscosity": 14, # in P*s as function of delta and tau
        "thermal_conductivity": 15, # in W/m/K as function of delta and tau
        "surface_tension": 16, # in N/m as function of delta and tau
    }

    parameters = [
        "R",  # specific ideal gas constant [kJ/kg/K]
        "MW",  # molecular weight [g/mol]
        "T_star",  # for calculating tau = T_star/T [K]
        "rho_star",  # for calculating delta = rho_star/rho [kg/m3]
        "Tc",  # critical temperature [K]
        "rhoc",  # critical density [kg/m3]
        "Pc",  # critical pressure [kPa]
        "Tt",  # triple point temperature [K]
        "Pt",  # triple point pressure [K]
        "rhot_l",  # liquid triple point density [kg/m3]
        "rhot_v",  # vapor triple point density [kg/m3]
        "P_min",  # minimum pressure [kPa]
        "P_max",  # maximum pressure [kPa]
        "rho_max",  # maximum density [kg/m3]
        "T_min",  # minimum temperature [kPa]
        "T_max",  # maximum temperature [kPa]
    ]

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
        self.m = pyo.ConcreteModel()
        self.m.delta = pyo.Var()  # rho / rho_star
        self.m.tau = pyo.Var()  # T_star / T
        self.m.p = pyo.Var()  # pressure
        self.m.R = pyo.Param(initialize=R)
        self.m.MW = pyo.Param(initialize=MW)
        self.m.T_star = pyo.Param(initialize=T_star)
        self.m.rho_star = pyo.Param(initialize=rho_star)
        self.m.Tc = pyo.Param(initialize=Tc)
        self.m.rhoc = pyo.Param(initialize=rhoc)
        self.m.Pc = pyo.Param(initialize=Pc)
        self.m.Tt = pyo.Param(initialize=Tt)
        self.m.Pt = pyo.Param(initialize=Pt)
        self.m.rhot_l = pyo.Param(initialize=rhot_l)
        self.m.rhot_v = pyo.Param(initialize=rhot_v)
        self.m.P_min = pyo.Param(initialize=P_min)
        self.m.P_max = pyo.Param(initialize=P_max)
        self.m.rho_max = pyo.Param(initialize=rho_max)
        self.m.T_min = pyo.Param(initialize=T_min)
        self.m.T_max = pyo.Param(initialize=T_max)

    @property
    def model(self):
        return self.m

    @model.setter
    def set_model(self, v):
        raise RuntimeError("model is not settable")

    def add(self, expressions):
        for name in self.expressions:
            try:
                expr = expressions[name]
            except:
                continue
            if callable(expr):
                setattr(self.m, name, pyo.Objective(rule=expr))
            else:
                setattr(self.m, name, pyo.Objective(expr=expr))

    def write(self):
        nl_file, smap_id = self.m.write(f"{self.comp}_expressions.nl")
        smap = self.model.solutions.symbol_map[smap_id]
        expr_map = [self.is_omitted_index] * len(self.expressions)
        var_map = [self.is_omitted_index] * len(self.variables)
        for s, c in smap.bySymbol.items():
            if s.startswith("o"):
                i = self.expressions[c().name]
                j = int(s[1:])
                expr_map[i] = j
            elif s.startswith("v"):
                i = self.variables[c().name]
                j = int(s[1:])
                var_map[i] = j
        param_dict = {
            "nl_file": nl_file,
            "expr_map": expr_map,
            "var_map": var_map,
            "param": {
                "R": pyo.value(self.m.R),
                "MW": pyo.value(self.m.MW),
                "T_star": pyo.value(self.m.T_star),
                "rho_star": pyo.value(self.m.rho_star),
                "Tc": pyo.value(self.m.Tc),
                "rhoc": pyo.value(self.m.rhoc),
                "Pc": pyo.value(self.m.Pc),
                "Tt": pyo.value(self.m.Tt),
                "Pt": pyo.value(self.m.Pt),
                "rhot_l": pyo.value(self.m.rhot_l),
                "rhot_v": pyo.value(self.m.rhot_v),
                "P_min": pyo.value(self.m.P_min),
                "P_max": pyo.value(self.m.P_max),
                "rho_max": pyo.value(self.m.rho_max),
                "T_min": pyo.value(self.m.T_min),
                "T_max": pyo.value(self.m.T_max),
            },
        }
        with open(f"{self.comp}_parameters.json", "w") as f:
            json.dump(param_dict, f, indent=4)

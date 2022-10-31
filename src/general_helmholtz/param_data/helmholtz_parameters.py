import pyomo.environ as pyo
import json


class WriteParameters(object):
    is_omitted_index = 1000

    required_variables = {
        "delta": 0,
        "tau": 1,
        "p": 2,
    }

    required_expressions = {
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
        "delta_v_sat_approx": 12,
        "delta_l_sat_approx": 13,
    }

    required_parameters = [
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

    @property
    def model(self):
        return self.m

    @model.setter
    def set_model(self, v):
        raise RuntimeError("model is not settable")

    def add(self, expressions):
        for name in self.required_expressions:
            try:
                expr = expressions[name]
            except:
                # raise RuntimeError(f"Missing expected expression '{name}'")
                continue
            setattr(self.m, name, pyo.Objective(expr=expr))

    def write(self):
        nl_file, smap_id = self.m.write(f"{self.comp}_expressions.nl", format="nl_v2")
        smap = self.model.solutions.symbol_map[smap_id]
        expr_map = [self.is_omitted_index] * len(self.required_expressions)
        var_map = [self.is_omitted_index] * len(self.required_variables)
        for s, c in smap.bySymbol.items():
            if s.startswith("o"):
                i = self.required_expressions[c().name]
                j = int(s[1:])
                expr_map[i] = j
            elif s.startswith("v"):
                i = self.required_variables[c().name]
                j = int(s[1:])
                var_map[i] = j
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
        with open(f"{self.comp}_parameters.json", "w") as f:
            json.dump(param_dict, f, indent=4)

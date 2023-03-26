import json
import numpy
import pyomo.environ as pyo
import logging

_log = logging.getLogger("idaes.helmholtz_parameters")


def _parse_int_key(pairs):
    d = {}
    for i, x in enumerate(pairs):
        try:
            if x[0][0] == "(" and x[0][-1] == ")":
                d[tuple(map(int, map(str.strip, "(4, 5)"[1:-1].split(","))))] = x[1]
            else:
                d[int(x[0])] = x[1]
        except ValueError:
            d[x[0]] = x[1]
    return d


class WriteParameters(object):

    variables = {
        "delta": 0,
        "tau": 1,
        "p": 2,
        "T": 3,
    }

    expressions = {  # these are expressions in the thermo model. Other models
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
        "delta_v_sat_approx": 12,  # approximate delta_sat_vapor(tau)
        "delta_l_sat_approx": 13,  # approximate delta_sat_liquid(tau)
    }

    optional_expressions = {
        "thermal_conductivity": "tcx",
        "surface_tension": "st",
        "viscosity": "visc",
    }

    def __init__(
        self,
        parameters,
    ):
        # if parameters is a string read from json
        if isinstance(parameters, str):
            with open(parameters, "r") as fp:
                parameters = json.load(fp, object_pairs_hook=_parse_int_key)
        self.parameters = parameters
        self.comp = parameters["comp"]
        self.R = parameters["basic"]["R"]
        self.MW = parameters["basic"]["MW"]
        self.T_star = parameters["basic"]["T_star"]
        self.rho_star = parameters["basic"]["rho_star"]
        self.Tc = parameters["basic"]["Tc"]
        self.rhoc = parameters["basic"]["rhoc"]
        self.Pc = parameters["basic"]["Pc"]
        self.Tt = parameters["basic"]["Tt"]
        self.Pt = parameters["basic"]["Pt"]
        self.rhot_l = parameters["basic"]["rhot_l"]
        self.rhot_v = parameters["basic"]["rhot_v"]
        self.P_min = parameters["basic"]["P_min"]
        self.P_max = parameters["basic"]["P_max"]
        self.rho_max = parameters["basic"]["rho_max"]
        self.T_min = parameters["basic"]["T_min"]
        self.T_max = parameters["basic"]["T_max"]

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

        phi_ideal_type = parameters["eos"].get("phi_ideal_type", 0)
        if phi_ideal_type > 0:
            self.add(
                phi_ideal_types[phi_ideal_type](model=self.model, parameters=parameters)
            )
        phi_residual_type = parameters["eos"].get("phi_residual_type", 0)
        if phi_residual_type > 0:
            self.add(
                phi_residual_types[phi_residual_type](
                    model=self.model, parameters=parameters
                )
            )

        aux_parameters = parameters.get("aux", None)
        if aux_parameters is not None:
            delta_l_sat_parameters = parameters["aux"].get("delta_l_sat_approx", None)
            delta_v_sat_parameters = parameters["aux"].get("delta_v_sat_approx", None)
            if delta_l_sat_parameters is not None:
                etype = delta_l_sat_parameters.get("type", 0)
                if etype:
                    self.add(
                        {
                            "delta_l_sat_approx": delta_sat_types[etype](
                                model=self.model,
                                name="delta_l_sat_approx",
                                parameters=parameters,
                            ),
                        }
                    )
            if delta_v_sat_parameters is not None:
                etype = delta_v_sat_parameters.get("type", 0)
                if etype:
                    self.add(
                        {
                            "delta_v_sat_approx": delta_sat_types[etype](
                                model=self.model,
                                name="delta_v_sat_approx",
                                parameters=parameters,
                            ),
                        }
                    )

        try:
            etype = parameters["transport"]["surface_tension"]["type"]
        except KeyError:  # No surface tension to add
            etype = 0
        if etype > 0:
            self.add(
                {
                    "surface_tension": surface_tension_types[etype](
                        model=self.model_st, parameters=parameters
                    ),
                }
            )

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
            if name in self.expressions:  # check if in thermo model
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

    def approx_sat_curves(self, trange, off=0):
        """This provides an easy way to verify that the approximate saturated density
        curves are correct, since they are not directly returned by any EOS property
        functions.
        """
        print("\n=====================================================================")
        print(" Check approx sat delta curves")
        print("=====================================================================")
        print(
            f"{'T [K]':7s}  {'T [C]':8s}  {'~rho_l [kg/m3]':14s}  {'~rho_v [kg/m3]':14s}  {'~P [kPa]':9s}"
        )
        print("---------------------------------------------------------------------")
        for T in trange:
            self.model.tau.fix(self.model.T_star / (T + off))
            delta_l = pyo.value(self.model.delta_l_sat_approx)
            delta_v = pyo.value(self.model.delta_v_sat_approx)
            rho_l = delta_l * self.model.rho_star
            rho_v = delta_v * self.model.rho_star
            self.model.delta.fix(delta_v)
            P_v = pyo.value(
                rho_v * self.R * (T + off) * (1 + self.model.delta * self.model.phir_d)
            )
            print(
                f"{T:7.3f}, {T - 273.15: 8.3f}, {rho_l:14.4f}, {rho_v:14.4f}, {P_v:9.4f}"
            )
        print("=====================================================================\n")

    def write_model(self, model, model_name, expressions=None):
        nl_file, smap_id = model.write(f"{self.comp}_expressions_{model_name}.nl")
        smap = model.solutions.symbol_map[smap_id]
        var_map = [1000] * 4
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
        print(
            "\n======================================================================="
        )
        print(f" Writing expressions for {self.comp}")
        print("=======================================================================")
        for name in self.expressions:
            if name not in self.has_expression:
                raise RuntimeError(f"Required expression {name} not provided.")
        nl_file, expr_map, var_map = self.write_model(
            self.model, "eos", self.expressions
        )
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

        trng = []
        linspc = numpy.linspace(self.Tt, self.T_star, 15)
        for i, t in enumerate(linspc):
            if i != 0 and i != len(linspc) - 1:
                t = round(t, 3)
            trng.append(t)
        self.approx_sat_curves(trng)


def surface_tension_type01(model, parameters):
    s = parameters["transport"]["surface_tension"]["s"]
    n = parameters["transport"]["surface_tension"]["n"]
    tc = parameters["transport"]["surface_tension"]["Tc"]
    return sum(s[i] * (1 - model.T_star / model.tau / tc) ** n[i] for i in s)


def sat_delta_type01(model, name, parameters):
    c = parameters["aux"][name]["c"]
    n = parameters["aux"][name]["n"]
    t = parameters["aux"][name]["t"]
    return c + sum(n[i] * (1 - 1 / model.tau) ** t[i] for i in n)


def sat_delta_type02(model, name, parameters):
    c = parameters["aux"][name]["c"]
    n = parameters["aux"][name]["n"]
    t = parameters["aux"][name]["t"]
    return c * pyo.exp(sum(n[i] * (1 - 1 / model.tau) ** t[i] for i in n))


def phi_ideal_expressions_type01(model, parameters):
    last_term = parameters["eos"]["last_term_ideal"]
    n0 = parameters["eos"]["n0"]
    g0 = parameters["eos"]["g0"]
    offset = parameters["eos"].get("reference_state_offset", {1:0.0, 2:0.0})
    print(offset)
    rng = range(4, last_term + 1)
    return {
        "phii": pyo.log(model.delta)
        + (n0[1] + offset[1])
        + (n0[2] + offset[2]) * model.tau
        + n0[3] * pyo.log(model.tau)
        + sum(n0[i] * pyo.log(1 - pyo.exp(-g0[i] * model.tau)) for i in rng),
        "phii_d": 1.0 / model.delta,
        "phii_dd": -1.0 / model.delta**2,
        "phii_t": (n0[2] + offset[2])
        + n0[3] / model.tau
        + sum(n0[i] * g0[i] / (pyo.exp(g0[i] * model.tau) - 1) for i in rng),
        "phii_tt": -n0[3] / model.tau**2
        - sum(
            n0[i]
            * g0[i] ** 2
            * pyo.exp(-g0[i] * model.tau)
            / (1 - pyo.exp(-g0[i] * model.tau)) ** 2
            for i in rng
        ),
        "phii_dt": 0,
    }


def phi_ideal_expressions_type02(model, parameters):
    last_term = parameters["eos"]["last_term_ideal"]
    n0 = parameters["eos"]["n0"]
    g0 = parameters["eos"]["g0"]
    offset = parameters["eos"].get("reference_state_offset", {1:0.0, 2:0.0})
    rng1 = range(4, last_term[0] + 1)
    rng2 = range(last_term[0] + 1, last_term[1] + 1)
    return {
        "phii": pyo.log(model.delta)
        + n0[1] + offset[1]
        + (n0[2] + offset[2]) * model.tau
        + n0[3] * pyo.log(model.tau)
        + sum(n0[i] * model.tau ** g0[i] for i in rng1)
        + sum(n0[i] * pyo.log(1 - pyo.exp(-g0[i] * model.tau)) for i in rng2),
        "phii_d": 1.0 / model.delta,
        "phii_dd": -1.0 / model.delta**2,
        "phii_t": (n0[2] + offset[2])
        + n0[3] / model.tau
        + sum(n0[i] * g0[i] * model.tau ** (g0[i] - 1) for i in rng1)
        + sum(n0[i] * g0[i] / (pyo.exp(g0[i] * model.tau) - 1) for i in rng2),
        "phii_tt": -n0[3] / model.tau**2
        + sum(n0[i] * g0[i] * (g0[i] - 1) * model.tau ** (g0[i] - 2) for i in rng1)
        - sum(
            n0[i]
            * g0[i] ** 2
            * pyo.exp(g0[i] * model.tau)
            / (pyo.exp(g0[i] * model.tau) - 1) ** 2
            for i in rng2
        ),
        "phii_dt": 0,
    }


def phi_ideal_expressions_type03(model, parameters):
    last_term = parameters["eos"]["last_term_ideal"]
    n0 = parameters["eos"]["n0"]
    g0 = parameters["eos"]["g0"]
    offset = parameters["eos"].get("reference_state_offset", {1:0.0, 2:0.0})
    rng1 = range(4, last_term + 1)
    return {
        "phii": pyo.log(model.delta)
        + n0[1] + offset[1]
        + (n0[2] + offset[2]) * model.tau
        + n0[3] * pyo.log(model.tau)
        + sum(n0[i] * model.tau ** g0[i] for i in rng1),
        "phii_d": 1.0 / model.delta,
        "phii_dd": -1.0 / model.delta**2,
        "phii_t": (n0[2] + offset[2])
        + n0[3] / model.tau
        + sum(n0[i] * g0[i] * model.tau ** (g0[i] - 1) for i in rng1),
        "phii_tt": -n0[3] / model.tau**2
        + sum(n0[i] * g0[i] * (g0[i] - 1) * model.tau ** (g0[i] - 2) for i in rng1),
        "phii_dt": 0,
    }


def phi_residual_expressions_type01(model, parameters):
    last_term = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    c = parameters["eos"]["c"]
    first_term = 1
    rng = []
    for last_term in last_term:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 1)
            * model.tau ** t[i]
            * (d[i] - c[i] * model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 2)
            * model.tau ** t[i]
            * (
                (d[i] - c[i] * model.delta ** c[i])
                * (d[i] - 1 - c[i] * model.delta ** c[i])
                - c[i] ** 2 * model.delta ** c[i]
            )
            for i in rng[1]
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 1)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * (t[i] - 1)
            * model.delta ** d[i]
            * model.tau ** (t[i] - 2)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_dt": sum(
            n[i] * t[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** (d[i] - 1)
            * model.tau ** (t[i] - 1)
            * (d[i] - c[i] * model.delta ** c[i])
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
    }


def phi_residual_expressions_type02(model, parameters):
    last_term = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    c = parameters["eos"]["c"]
    a = parameters["eos"]["a"]
    b = parameters["eos"]["b"]
    e = parameters["eos"]["e"]
    g = parameters["eos"]["g"]
    first_term = 1
    rng = []
    for last_term in last_term:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            for i in rng[2]
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 1)
            * model.tau ** t[i]
            * (d[i] - c[i] * model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (d[i] / model.delta - 2 * a[i] * (model.delta - e[i]))
            for i in rng[2]
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 2)
            * model.tau ** t[i]
            * (
                (d[i] - c[i] * model.delta ** c[i])
                * (d[i] - 1 - c[i] * model.delta ** c[i])
                - c[i] ** 2 * model.delta ** c[i]
            )
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (
                -2.0 * a[i] * model.delta ** d[i]
                + 4 * a[i] ** 2 * model.delta ** d[i] * (model.delta - e[i]) ** 2
                - 4 * d[i] * a[i] * model.delta ** (d[i] - 1) * (model.delta - e[i])
                + d[i] * (d[i] - 1) * model.delta ** (d[i] - 2)
            )
            for i in rng[2]
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 1)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (t[i] / model.tau - 2 * b[i] * (model.tau - g[i]))
            for i in rng[2]
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * (t[i] - 1)
            * model.delta ** d[i]
            * model.tau ** (t[i] - 2)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (
                (t[i] / model.tau - 2 * b[i] * (model.tau - g[i])) ** 2
                - t[i] / model.tau**2
                - 2 * b[i]
            )
            for i in rng[2]
        ),
        "phir_dt": sum(
            n[i] * t[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** (d[i] - 1)
            * model.tau ** (t[i] - 1)
            * (d[i] - c[i] * model.delta ** c[i])
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (d[i] / model.delta - 2.0 * a[i] * (model.delta - e[i]))
            * (t[i] / model.tau - 2.0 * b[i] * (model.tau - g[i]))
            for i in rng[2]
        ),
    }


def phi_residual_expressions_type03(model, parameters):
    last_term = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    c = parameters["eos"]["c"]
    b = parameters["eos"]["b"]
    first_term = 1
    rng = []
    for last_term in last_term:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** c[i])
            * pyo.exp(-model.tau ** b[i])
            for i in rng[2]
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 1)
            * model.tau ** t[i]
            * (d[i] - c[i] * model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 1)
            * model.tau ** t[i]
            * (d[i] - c[i] * model.delta ** c[i])
            * pyo.exp(-model.tau ** b[i])
            for i in rng[2]
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 2)
            * model.tau ** t[i]
            * (
                (d[i] - c[i] * model.delta ** c[i])
                * (d[i] - 1 - c[i] * model.delta ** c[i])
                - c[i] ** 2 * model.delta ** c[i]
            )
            for i in rng[1]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 2)
            * model.tau ** t[i]
            * (
                (d[i] - c[i] * model.delta ** c[i])
                * (d[i] - 1 - c[i] * model.delta ** c[i])
                - c[i] ** 2 * model.delta ** c[i]
            )
            * pyo.exp(-model.tau ** b[i])
            for i in rng[2]
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 1)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 1)
            * pyo.exp(-model.delta ** c[i])
            * pyo.exp(-model.tau ** b[i])
            * (t[i] - b[i] * model.tau ** b[i])
            for i in rng[2]
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * (t[i] - 1)
            * model.delta ** d[i]
            * model.tau ** (t[i] - 2)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 2)
            * pyo.exp(-model.delta ** c[i])
            * pyo.exp(-model.tau ** b[i])
            * (
                (t[i] - b[i] * model.tau ** b[i])
                * (t[i] - 1 - b[i] * model.tau ** b[i])
                - b[i] ** 2 * model.tau ** b[i]
            )
            for i in rng[2]
        ),
        "phir_dt": sum(
            n[i] * t[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** (d[i] - 1)
            * model.tau ** (t[i] - 1)
            * (d[i] - c[i] * model.delta ** c[i])
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** (d[i] - 1)
            * model.tau ** (t[i] - 1)
            * (d[i] - c[i] * model.delta ** c[i])
            * (t[i] - b[i] * model.tau ** b[i])
            * pyo.exp(-model.delta ** c[i])
            * pyo.exp(-model.tau ** b[i])
            for i in rng[2]
        ),
    }


def phi_residual_expressions_type04(model, parameters):
    last_term = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    first_term = 1
    rng = []
    for last_term in last_term:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            pyo.exp(-model.delta**k)
            * sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[k])
            for k in range(1, len(rng))
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i]
                * (d[i] - k * model.delta**k)
                * model.delta ** (d[i] - 1)
                * model.tau ** t[i]
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i]
                * (
                    -k * model.delta**k * (d[i] - k * model.delta**k)
                    + (d[i] - 1) * (d[i] - k * model.delta**k)
                    - k**2 * model.delta**k
                )
                * model.delta ** (d[i] - 2)
                * model.tau ** t[i]
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1)
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_dt": sum(
            n[i] * d[i] * t[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i]
                * t[i]
                * (d[i] - k * model.delta**k)
                * model.delta ** (d[i] - 1)
                * model.tau ** (t[i] - 1)
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
    }


phi_residual_types = {
    0: None,  # custom
    1: phi_residual_expressions_type01,
    2: phi_residual_expressions_type02,
    3: phi_residual_expressions_type03,
    4: phi_residual_expressions_type04,
}

phi_ideal_types = {
    0: None,  # custom
    1: phi_ideal_expressions_type01,
    2: phi_ideal_expressions_type02,
    3: phi_ideal_expressions_type03,
}

delta_sat_types = {
    0: None,  # custom
    1: sat_delta_type01,
    2: sat_delta_type02,
}

surface_tension_types = {
    0: None,
    1: surface_tension_type01,
}

##################################################################################
#                                                                                #
# R227EA EOS Expressions and Parameters:                                         #
#                                                                                #
# Eric W. Lemmon and Roland Span. Thermodynamic Properties of R-227ea,           #
#     R-365mfc, R-115, and R13I1. J. Chem. Eng. Data, 2016, submitted.           #
#     doi:10.1021/acs.jced.5b00684.                                              #
#                                                                                #
##################################################################################

import pyomo.environ as pyo
import json
from helmholtz_parameters import WriteParameters


def main():
    m = pyo.ConcreteModel()
    m.delta = pyo.Var()
    m.tau = pyo.Var()

    c = {
        6: 1,
        7: 1,
        8: 1,
        9: 1,
        10: 2,
        11: 2,
    }

    d = {
        1: 1,
        2: 1,
        3: 2,
        4: 2,
        5: 4,
        6: 1,
        7: 3,
        8: 6,
        9: 6,
        10: 2,
        11: 3,
        12: 1,
        13: 2,
        14: 1,
        15: 1,
        16: 4,
        17: 2,
        18: 1,
    }

    t = {
        1: 0.34,
        2: 0.77,
        3: 0.36,
        4: 0.90,
        5: 1.00,
        6: 2.82,
        7: 2.10,
        8: 0.90,
        9: 1.13,
        10: 3.80,
        11: 2.75,
        12: 1.50,
        13: 2.50,
        14: 2.50,
        15: 5.40,
        16: 4.00,
        17: 1.00,
        18: 3.50,
    }

    n = {
        1: 2.024341,
        2: -2.605930,
        3: 0.4957216,
        4: -0.8240820,
        5: 0.06543703,
        6: -1.024610,
        7: 0.6247065,
        8: 0.2997521,
        9: -0.3539170,
        10: -1.232043,
        11: -0.8824483,
        12: 0.1349661,
        13: -0.2662928,
        14: 0.1764733,
        15: 0.01536163,
        16: -0.004667185,
        17: -11.70854,
        18: 0.9114512,
    }

    a = {
        12: 0.83,
        13: 2.19,
        14: 2.44,
        15: 3.65,
        16: 8.88,
        17: 8.23,
        18: 2.01,
    }

    b = {
        12: 1.72,
        13: 5.20,
        14: 2.31,
        15: 1.02,
        16: 5.63,
        17: 50.9,
        18: 1.56,
    }

    g = {
        12: 0.414,
        13: 1.051,
        14: 1.226,
        15: 1.700,
        16: 0.904,
        17: 1.420,
        18: 0.926,
    }

    e = {
        12: 1.13,
        13: 0.71,
        14: 1.20,
        15: 1.70,
        16: 0.546,
        17: 0.896,
        18: 0.747,
    }

    n0 = {
        1: -15.8291124137,  # aka a0[1]
        2: 11.0879509962,  # aka a0[2]
        3: 4.0 - 1,  # aka c0 - 1
        4: 11.43,  # aka v[1]
        5: 12.83,  # aka v[2]
    }

    g0 = {  # aka u[k] / Tc
        4: 403.0 / 374.9,
        5: 1428.0 / 374.9,
    }

    we = WriteParameters(
        comp="r227ea",
        R=8.3144621 / 170.02886,
        MW=170.02886,
        T_star=374.9,
        rho_star=3.495 * 170.02886,
        Tc=374.9,
        rhoc=3.495 * 170.02886,
        Pc=2925,
        Tt=146.35,
        Pt=7.33e-3,
        rhot_l=1878.3,
        rhot_v=0.0010245,
        P_min=1e-9,
        P_max=5e5,
        rho_max=2500.0,
        T_min=145,
        T_max=500,
    )
    m = we.model
    we.add(
        {
            "phii": pyo.log(m.delta)
            + n0[1]
            + n0[2] * m.tau
            + n0[3] * pyo.log(m.tau)
            + sum(
                n0[i] * pyo.log(1 - pyo.exp(-g0[i] * m.tau)) for i in range(4, 5 + 1)
            ),
            "phii_d": 1.0 / m.delta,
            "phii_dd": -1.0 / m.delta**2,
            "phii_t": n0[2]
            + n0[3] / m.tau
            + sum(
                n0[i] * g0[i] * ((1 - pyo.exp(-g0[i] * m.tau)) ** (-1) - 1)
                for i in range(4, 5 + 1)
            ),
            "phii_tt": -n0[3] / m.tau**2
            - sum(
                n0[i]
                * g0[i] ** 2
                * pyo.exp(-g0[i] * m.tau)
                * (1 - pyo.exp(-g0[i] * m.tau)) ** (-2)
                for i in range(4, 5 + 1)
            ),
            "phii_dt": 0,
            "phir": sum(n[i] * m.delta ** d[i] * m.tau ** t[i] for i in range(1, 5 + 1))
            + sum(
                n[i] * m.delta ** d[i] * m.tau ** t[i] * pyo.exp(-m.delta ** c[i])
                for i in range(6, 11 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                for i in range(12, 18 + 1)
            ),
            "phir_d": sum(
                n[i] * d[i] * m.delta ** (d[i] - 1) * m.tau ** t[i]
                for i in range(1, 5 + 1)
            )
            + sum(
                n[i]
                * pyo.exp(-m.delta ** c[i])
                * m.delta ** (d[i] - 1)
                * m.tau ** t[i]
                * (d[i] - c[i] * m.delta ** c[i])
                for i in range(6, 11 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (d[i] / m.delta - 2 * a[i] * (m.delta - e[i]))
                for i in range(12, 18 + 1)
            ),
            "phir_dd": sum(
                n[i] * d[i] * (d[i] - 1) * m.delta ** (d[i] - 2) * m.tau ** t[i]
                for i in range(1, 5 + 1)
            )
            + sum(
                n[i]
                * pyo.exp(-m.delta ** c[i])
                * m.delta ** (d[i] - 2)
                * m.tau ** t[i]
                * (
                    (d[i] - c[i] * m.delta ** c[i])
                    * (d[i] - 1 - c[i] * m.delta ** c[i])
                    - c[i] ** 2 * m.delta ** c[i]
                )
                for i in range(6, 11 + 1)
            )
            + sum(
                n[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (
                    -2.0 * a[i] * m.delta ** d[i]
                    + 4 * a[i] ** 2 * m.delta ** d[i] * (m.delta - e[i]) ** 2
                    - 4 * d[i] * a[i] * m.delta ** (d[i] - 1) * (m.delta - e[i])
                    + d[i] * (d[i] - 1) * m.delta ** (d[i] - 2)
                )
                for i in range(12, 18 + 1)
            ),
            "phir_t": sum(
                n[i] * t[i] * m.delta ** d[i] * m.tau ** (t[i] - 1)
                for i in range(1, 5 + 1)
            )
            + sum(
                n[i]
                * t[i]
                * m.delta ** d[i]
                * m.tau ** (t[i] - 1)
                * pyo.exp(-m.delta ** c[i])
                for i in range(6, 11 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (t[i] / m.tau - 2 * b[i] * (m.tau - g[i]))
                for i in range(12, 18 + 1)
            ),
            "phir_tt": sum(
                n[i] * t[i] * (t[i] - 1) * m.delta ** d[i] * m.tau ** (t[i] - 2)
                for i in range(1, 5 + 1)
            )
            + sum(
                n[i]
                * t[i]
                * (t[i] - 1)
                * m.delta ** d[i]
                * m.tau ** (t[i] - 2)
                * pyo.exp(-m.delta ** c[i])
                for i in range(6, 11 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (
                    (t[i] / m.tau - 2 * b[i] * (m.tau - g[i])) ** 2
                    - t[i] / m.tau**2
                    - 2 * b[i]
                )
                for i in range(12, 18 + 1)
            ),
            "phir_dt": sum(
                n[i] * t[i] * d[i] * m.delta ** (d[i] - 1) * m.tau ** (t[i] - 1)
                for i in range(1, 5 + 1)
            )
            + sum(
                n[i]
                * t[i]
                * m.delta ** (d[i] - 1)
                * m.tau ** (t[i] - 1)
                * (d[i] - c[i] * m.delta ** c[i])
                * pyo.exp(-m.delta ** c[i])
                for i in range(6, 11 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (d[i] / m.delta - 2.0 * a[i] * (m.delta - e[i]))
                * (t[i] / m.tau - 2.0 * b[i] * (m.tau - g[i]))
                for i in range(12, 18 + 1)
            ),
            "delta_v_sat_approx": pyo.exp(
                -109.367 * (1 - 1 / m.tau) ** 0.64
                + 322.88 * (1 - 1 / m.tau) ** 0.77
                - 485.87 * (1 - 1 / m.tau) ** 0.96
                + 417.10 * (1 - 1 / m.tau) ** 1.2
                - 174.52 * (1 - 1 / m.tau) ** 1.45
                - 52.695 * (1 - 1 / m.tau) ** 5.35
                - 114.41 * (1 - 1 / m.tau) ** 12.0
            ),
            "delta_l_sat_approx": (
                1.00
                - 0.29926 * (1 - 1 / m.tau) ** 0.15
                + 2.8025 * (1 - 1 / m.tau) ** 0.3
                - 1.9602 * (1 - 1 / m.tau) ** 0.44
                + 2.0784 * (1 - 1 / m.tau) ** 0.6
                + 0.21701 * (1 - 1 / m.tau) ** 2.75
            ),
        }
    )
    we.write()


if __name__ == "__main__":
    main()

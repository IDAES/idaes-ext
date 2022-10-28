import pyomo.environ as pyo
import json
from helmholtz_parameters import WriteParameters


def main():
    m = pyo.ConcreteModel()
    m.delta = pyo.Var()
    m.tau = pyo.Var()

    d = {
        1: 2,
        2: 1,
        3: 3,
        4: 6,
        5: 6,
        6: 1,
        7: 1,
        8: 2,
        9: 5,
        10: 2,
        11: 2,
        12: 4,
        13: 1,
        14: 4,
        15: 1,
        16: 2,
        17: 4,
        18: 1,
        19: 5,
        20: 3,
        21: 10,
    }

    t = {
        1: -0.5,
        2: 0,
        3: 0,
        4: 0,
        5: 1.5,
        6: 1.5,
        7: 2,
        8: 2,
        9: 1,
        10: 3,
        11: 5,
        12: 1,
        13: 5,
        14: 5,
        15: 6,
        16: 10,
        17: 10,
        18: 10,
        19: 18,
        20: 22,
        21: 50,
    }

    n = {
        1: 0.5586817e-1,
        2: 0.4982230e0,
        3: 0.2458698e-1,
        4: 0.8570145e-3,
        5: 0.4788584e-3,
        6: -0.1800808e1,
        7: 0.2671641e0,
        8: -0.4781652e-1,
        9: 0.1423987e-1,
        10: 0.3324062e0,
        11: -0.7485907e-2,
        12: 0.1017263e-3,
        13: -0.5184567e0,
        14: -0.8692288e-1,
        15: 0.2057144e0,
        16: -0.5000457e-2,
        17: 0.4603262e-3,
        18: -0.3497836e-2,
        19: 0.6995038e-2,
        20: -0.1452184e-1,
        21: -0.1285458e-3,
    }

    n0 = {
        1: -1.019535,
        2: 9.047135,
        3: -1.629789,
        4: -9.723916,
        5: -3.927170,
    }

    g0 = {
        4: -1.0 / 2.0,
        5: -3.0 / 4.0,
    }

    we = WriteParameters(
        comp="r134a",
        R=0.081488856,
        MW=102.032,
        T_star=374.18,
        rho_star=508.0,
        Tc=374.21,
        rhoc=511.95,
        Pc=4059.11,
        Tt=169.85,
        Pt=0.391,
        rhot_l=1591.1,
        rhot_v=0.028172,
        P_min=1e-9,
        P_max=7e4,
        rho_max=2500.0,
        T_min=168,
        T_max=460,
    )
    nk = [8, 11, 17, 20, 21]

    print(list(range(1, nk[0] + 1)))
    for k in range(1, 4 + 1):
        print(k)
        print(list(range(nk[k - 1] + 1, nk[k] + 1)))

    m = we.model
    we.add(
        {
            "phii": pyo.log(m.delta)
            + n0[1]
            + n0[2] * m.tau
            + n0[3] * pyo.log(m.tau)
            + sum(n0[i] * m.tau ** g0[i] for i in range(4, 5 + 1)),
            "phii_d": 1.0 / m.delta,
            "phii_dd": -1.0 / m.delta**2,
            "phii_t": n0[2]
            + n0[3] / m.tau
            + sum(n0[i] * g0[i] * m.tau ** (g0[i] - 1) for i in range(4, 5 + 1)),
            "phii_tt": -n0[3] / m.tau**2
            + sum(
                n0[i] * g0[i] * (g0[i] - 1) * m.tau ** (g0[i] - 2)
                for i in range(4, 5 + 1)
            ),
            "phii_dt": 0,
            "phir": sum(n[i] * m.delta ** d[i] * m.tau ** t[i] for i in range(1, 8 + 1))
            + pyo.exp(-m.delta)
            * sum(n[i] * m.delta ** d[i] * m.tau ** t[i] for i in range(9, 11 + 1))
            + pyo.exp(-m.delta**2)
            * sum(n[i] * m.delta ** d[i] * m.tau ** t[i] for i in range(12, 17 + 1))
            + pyo.exp(-m.delta**3)
            * sum(n[i] * m.delta ** d[i] * m.tau ** t[i] for i in range(18, 20 + 1))
            + pyo.exp(-m.delta**4) * n[21] * m.delta ** d[21] * m.tau ** t[21],
            "phir_d": sum(
                n[i] * d[i] * m.delta ** (d[i] - 1) * m.tau ** t[i]
                for i in range(1, nk[0] + 1)
            )
            + sum(
                pyo.exp(-m.delta**k)
                * sum(
                    n[i]
                    * (d[i] - k * m.delta**k)
                    * m.delta ** (d[i] - 1)
                    * m.tau ** t[i]
                    for i in range(nk[k - 1] + 1, nk[k] + 1)
                )
                for k in range(1, 4 + 1)
            ),
            "phir_dd": sum(
                n[i] * d[i] * (d[i] - 1) * m.delta ** (d[i] - 2) * m.tau ** t[i]
                for i in range(1, nk[0] + 1)
            )
            + sum(
                pyo.exp(-m.delta**k)
                * sum(
                    n[i]
                    * (
                        -k * m.delta**k * (d[i] - k * m.delta**k)
                        + (d[i] - 1) * (d[i] - k * m.delta**k)
                        - k**2 * m.delta**k
                    )
                    * m.delta ** (d[i] - 2)
                    * m.tau ** t[i]
                    for i in range(nk[k - 1] + 1, nk[k] + 1)
                )
                for k in range(1, 4 + 1)
            ),
            "phir_t": sum(
                n[i] * t[i] * m.delta ** d[i] * m.tau ** (t[i] - 1)
                for i in range(1, nk[0] + 1)
            )
            + sum(
                pyo.exp(-m.delta**k)
                * sum(
                    n[i] * t[i] * m.delta ** d[i] * m.tau ** (t[i] - 1)
                    for i in range(nk[k - 1] + 1, nk[k] + 1)
                )
                for k in range(1, 4 + 1)
            ),
            "phir_tt": sum(
                n[i] * t[i] * (t[i] - 1) * m.delta ** d[i] * m.tau ** (t[i] - 2)
                for i in range(1, nk[0] + 1)
            )
            + sum(
                pyo.exp(-m.delta**k)
                * sum(
                    n[i] * t[i] * (t[i] - 1) * m.delta ** d[i] * m.tau ** (t[i] - 2)
                    for i in range(nk[k - 1] + 1, nk[k] + 1)
                )
                for k in range(1, 4 + 1)
            ),
            "phir_dt": sum(
                n[i] * d[i] * t[i] * m.delta ** (d[i] - 1) * m.tau ** (t[i] - 1)
                for i in range(1, nk[0] + 1)
            )
            + sum(
                pyo.exp(-m.delta**k)
                * sum(
                    n[i]
                    * t[i]
                    * (d[i] - k * m.delta**k)
                    * m.delta ** (d[i] - 1)
                    * m.tau ** (t[i] - 1)
                    for i in range(nk[k - 1] + 1, nk[k] + 1)
                )
                for k in range(1, 4 + 1)
            ),
            "delta_v_sat_approx": 516.86
            / we.rho_star
            * pyo.exp(
                -2.837294 * (1 - 1 / m.tau) ** (1.0 / 3.0)
                - 7.875988 * (1 - 1 / m.tau) ** (2.0 / 3.0)
                + 4.478586 * (1 - 1 / m.tau) ** (1.0 / 2.0)
                - 14.140125 * (1 - 1 / m.tau) ** (9.0 / 4.0)
                - 52.361297 * (1 - 1 / m.tau) ** (11.0 / 2.0)
            ),
            "delta_l_sat_approx": (
                518.20
                + 884.13 * (1 - 1 / m.tau) ** (1.0 / 3.0)
                + 485.84 * (1 - 1 / m.tau) ** (2.0 / 3.0)
                + 193.29 * (1 - 1 / m.tau) ** (10.0 / 3.0)
            )
            / we.rho_star,
        }
    )
    we.write()


if __name__ == "__main__":
    main()

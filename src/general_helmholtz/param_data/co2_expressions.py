##################################################################################
#                                                                                #
# CO2 EOS Expressions and Parameters:                                            #
#                                                                                #
# Span, R., and W. Wanger (1996). "A New Equation of State for Carbon Dioxide    #
#     Covering the Fluid Region from the Triple-Point Temperature to 1100 K as   #
#     Pressures up to 800 MPa." Journal of Physical and Chemical Reference Data, #
#     25, 1509.                                                                  #
# Vesovic, V., W.A. Wakeham, G.A. Olchowy, J.V. Sengers, J.T.R. Watson, J.       #
#     Millat, (1990). "The transport properties of carbon dioxide." J. Phys.     #
#     Chem. Ref. Data, 19, 763-808.                                              #
# Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon       #
#     Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.                             #
#                                                                                #
##################################################################################

import pyomo.environ as pyo
from helmholtz_parameters import WriteParameters


def thermal_conductivity_rule(m):
    b = {
        0: 0.4226159,
        1: 0.6280115,
        2: -0.5387661,
        3: 0.6735941,
        4: 0,
        5: 0,
        6: -0.4362677,
        7: 0.2255388,
    }
    c = {
        1: 2.387869e-2,
        2: 4.350794,
        3: -10.33404,
        4: 7.981590,
        5: -1.940558,
    }
    d = {
        1: 2.447164e-5,
        2: 8.705605e-8,
        3: -6.547950e-11,
        4: 6.594919e-14,
    }
    T = m.T_star / m.tau
    Ts = T / 251.196
    G = sum(b[i] / Ts**i for i in b)
    cint_over_k = 1.0 + pyo.exp(-183.5 / T) * sum(
        c[i] * (T / 100) ** (2 - i) for i in c
    )
    return (
        (
            0.475598 * pyo.sqrt(T) * (1 + 2.0 / 5.0 * cint_over_k) / G
            + d[1] * m.delta
            + d[2] * m.delta**2
            + d[3] * m.delta**3
            + d[4] * m.delta**4
        )
        / 1e3
    )

def viscosity_rule(m):
    a = {
        0: 0.235156,
        1: -0.491266,
        2: 5.211155e-2,
        3: 5.347906e-2,
        4: -1.537102e-2,
    }
    d = {
        1: 0.4071119e-8 * pyo.value(m.rho_star),
        2: 0.7198037e-10 * pyo.value(m.rho_star) ** 2,
        3: 0.2411697e-22 * pyo.value(m.rho_star) ** 6,
        4: 0.2971072e-28 * pyo.value(m.rho_star) ** 8,
        5: -0.1627888e-28 * pyo.value(m.rho_star) ** 8,
    }
    T = m.T_star / m.tau
    Ts = T / 251.196
    return (
        (
            1.00697
            * pyo.sqrt(T)
            / pyo.exp(sum(a[i] * pyo.log(Ts) ** i for i in a))
            / 1e6
            + d[1] * m.delta
            + d[2] * m.delta**2
            + d[3] * m.delta**6 / Ts**3
            + d[4] * m.delta**8
            + d[5] * m.delta**8 / Ts
        )
    )

def main():
    m = pyo.ConcreteModel()
    m.delta = pyo.Var()
    m.tau = pyo.Var()

    c = {
        8: 1,
        9: 1,
        10: 1,
        11: 1,
        12: 1,
        13: 1,
        14: 1,
        15: 1,
        16: 1,
        17: 2,
        18: 2,
        19: 2,
        20: 2,
        21: 2,
        22: 2,
        23: 2,
        24: 3,
        25: 3,
        26: 3,
        27: 4,
        28: 4,
        29: 4,
        30: 4,
        31: 4,
        32: 4,
        33: 5,
        34: 6,
    }

    d = {
        1: 1,
        2: 1,
        3: 1,
        4: 1,
        5: 2,
        6: 2,
        7: 3,
        8: 1,
        9: 2,
        10: 4,
        11: 5,
        12: 5,
        13: 5,
        14: 6,
        15: 6,
        16: 6,
        17: 1,
        18: 1,
        19: 4,
        20: 4,
        21: 4,
        22: 7,
        23: 8,
        24: 2,
        25: 3,
        26: 3,
        27: 5,
        28: 5,
        29: 6,
        30: 7,
        31: 8,
        32: 10,
        33: 4,
        34: 8,
        35: 2,
        36: 2,
        37: 2,
        38: 3,
        39: 3,
    }

    t = {
        1: 0.00,
        2: 0.75,
        3: 1.00,
        4: 2.00,
        5: 0.75,
        6: 2.00,
        7: 0.75,
        8: 1.50,
        9: 1.50,
        10: 2.50,
        11: 0.00,
        12: 1.50,
        13: 2.00,
        14: 0.00,
        15: 1.00,
        16: 2.00,
        17: 3.00,
        18: 6.00,
        19: 3.00,
        20: 6.00,
        21: 8.00,
        22: 6.00,
        23: 0.00,
        24: 7.00,
        25: 12.00,
        26: 16.00,
        27: 22.00,
        28: 24.00,
        29: 16.00,
        30: 24.00,
        31: 8.00,
        32: 2.00,
        33: 28.00,
        34: 14.00,
        35: 1.00,
        36: 0.00,
        37: 1.00,
        38: 3.00,
        39: 3.00,
    }

    n = {
        1: 3.88568232031610e-01,
        2: 2.93854759427400e00,
        3: -5.58671885349340e00,
        4: -7.67531995924770e-01,
        5: 3.17290055804160e-01,
        6: 5.48033158977670e-01,
        7: 1.22794112203350e-01,
        8: 2.16589615432200e00,
        9: 1.58417351097240e00,
        10: -2.31327054055030e-01,
        11: 5.81169164314360e-02,
        12: -5.53691372053820e-01,
        13: 4.89466159094220e-01,
        14: -2.42757398435010e-02,
        15: 6.24947905016780e-02,
        16: -1.21758602252460e-01,
        17: -3.70556852700860e-01,
        18: -1.67758797004260e-02,
        19: -1.19607366379870e-01,
        20: -4.56193625087780e-02,
        21: 3.56127892703460e-02,
        22: -7.44277271320520e-03,
        23: -1.73957049024320e-03,
        24: -2.18101212895270e-02,
        25: 2.43321665592360e-02,
        26: -3.74401334234630e-02,
        27: 1.43387157568780e-01,
        28: -1.34919690832860e-01,
        29: -2.31512250534800e-02,
        30: 1.23631254929010e-02,
        31: 2.10583219729400e-03,
        32: -3.39585190263680e-04,
        33: 5.59936517715920e-03,
        34: -3.03351180556460e-04,
        35: -2.13654886883200e02,
        36: 2.66415691492720e04,
        37: -2.40272122045570e04,
        38: -2.83416034239990e02,
        39: 2.12472844001790e02,
        40: -6.66422765407510e-01,
        41: 7.26086323498970e-01,
        42: 5.50686686128420e-02,
    }

    a = {
        35: 25,
        36: 25,
        37: 25,
        38: 15,
        39: 20,
    }

    b = {
        35: 325,
        36: 300,
        37: 300,
        38: 275,
        39: 275,
    }

    g = {
        35: 1.16,
        36: 1.19,
        37: 1.19,
        38: 1.25,
        39: 1.22,
    }

    e = {
        35: 1,
        36: 1,
        37: 1,
        38: 1,
        39: 1,
    }

    n0 = {
        1: 8.37304456,
        2: -3.70454304,
        3: 2.5,
        4: 1.99427042,
        5: 0.62105248,
        6: 0.41195293,
        7: 1.04028922,
        8: 0.08327678,
    }

    g0 = {
        4: 3.15163,
        5: 6.11190,
        6: 6.77708,
        7: 11.32384,
        8: 27.08792,
    }

    we = WriteParameters(
        comp="co2",
        R=0.1889241,
        MW=44.0098,
        T_star=304.1282,
        rho_star=467.6,
        Tc=304.1282,
        rhoc=467.6,
        Pc=7377.3,
        Tt=216.592,
        Pt=517.95,
        rhot_l=1178.46,
        rhot_v=13.761,
        P_min=1e-9,
        P_max=5e5,
        rho_max=1600.0,
        T_min=216,
        T_max=1400,
    )
    m = we.model
    we.add(
        {
            "phii": pyo.log(m.delta)
            + n0[1]
            + n0[2] * m.tau
            + n0[3] * pyo.log(m.tau)
            + sum(
                n0[i] * pyo.log(1 - pyo.exp(-g0[i] * m.tau)) for i in range(4, 8 + 1)
            ),
            "phii_d": 1.0 / m.delta,
            "phii_dd": -1.0 / m.delta**2,
            "phii_t": n0[2]
            + n0[3] / m.tau
            + sum(
                n0[i] * g0[i] * ((1 - pyo.exp(-g0[i] * m.tau)) ** (-1) - 1)
                for i in range(4, 8 + 1)
            ),
            "phii_tt": -n0[3] / m.tau**2
            - sum(
                n0[i]
                * g0[i] ** 2
                * pyo.exp(-g0[i] * m.tau)
                * (1 - pyo.exp(-g0[i] * m.tau)) ** (-2)
                for i in range(4, 8 + 1)
            ),
            "phii_dt": 0,
            "phir": sum(n[i] * m.delta ** d[i] * m.tau ** t[i] for i in range(1, 7 + 1))
            + sum(
                n[i] * m.delta ** d[i] * m.tau ** t[i] * pyo.exp(-m.delta ** c[i])
                for i in range(8, 34 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                for i in range(35, 39 + 1)
            ),
            "phir_d": sum(
                n[i] * d[i] * m.delta ** (d[i] - 1) * m.tau ** t[i]
                for i in range(1, 7 + 1)
            )
            + sum(
                n[i]
                * pyo.exp(-m.delta ** c[i])
                * m.delta ** (d[i] - 1)
                * m.tau ** t[i]
                * (d[i] - c[i] * m.delta ** c[i])
                for i in range(8, 34 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (d[i] / m.delta - 2 * a[i] * (m.delta - e[i]))
                for i in range(35, 39 + 1)
            ),
            "phir_dd": sum(
                n[i] * d[i] * (d[i] - 1) * m.delta ** (d[i] - 2) * m.tau ** t[i]
                for i in range(1, 7 + 1)
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
                for i in range(8, 34 + 1)
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
                for i in range(35, 39 + 1)
            ),
            "phir_t": sum(
                n[i] * t[i] * m.delta ** d[i] * m.tau ** (t[i] - 1)
                for i in range(1, 7 + 1)
            )
            + sum(
                n[i]
                * t[i]
                * m.delta ** d[i]
                * m.tau ** (t[i] - 1)
                * pyo.exp(-m.delta ** c[i])
                for i in range(8, 34 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (t[i] / m.tau - 2 * b[i] * (m.tau - g[i]))
                for i in range(35, 39 + 1)
            ),
            "phir_tt": sum(
                n[i] * t[i] * (t[i] - 1) * m.delta ** d[i] * m.tau ** (t[i] - 2)
                for i in range(1, 7 + 1)
            )
            + sum(
                n[i]
                * t[i]
                * (t[i] - 1)
                * m.delta ** d[i]
                * m.tau ** (t[i] - 2)
                * pyo.exp(-m.delta ** c[i])
                for i in range(8, 34 + 1)
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
                for i in range(35, 39 + 1)
            ),
            "phir_dt": sum(
                n[i] * t[i] * d[i] * m.delta ** (d[i] - 1) * m.tau ** (t[i] - 1)
                for i in range(1, 7 + 1)
            )
            + sum(
                n[i]
                * t[i]
                * m.delta ** (d[i] - 1)
                * m.tau ** (t[i] - 1)
                * (d[i] - c[i] * m.delta ** c[i])
                * pyo.exp(-m.delta ** c[i])
                for i in range(8, 34 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (d[i] / m.delta - 2.0 * a[i] * (m.delta - e[i]))
                * (t[i] / m.tau - 2.0 * b[i] * (m.tau - g[i]))
                for i in range(35, 39 + 1)
            ),
            "delta_v_sat_approx": pyo.exp(
                -1.7074879 * (1 - 1.0 / m.tau) ** 0.340
                - 0.82274670 * (1 - 1.0 / m.tau) ** 0.5
                - 4.6008549 * (1 - 1.0 / m.tau) ** 1.0
                - 10.111178 * (1 - 1.0 / m.tau) ** 7.0 / 3.0
                - 29.742252 * (1 - 1.0 / m.tau) ** 14.0 / 3.0
            ),
            "delta_l_sat_approx": pyo.exp(
                1.9245108 * (1 - 1.0 / m.tau) ** 0.34
                - 0.62385555 * (1 - 1.0 / m.tau) ** 0.5
                - 0.32731127 * (1 - 1.0 / m.tau) ** 10.0 / 6.0
                + 0.39245142 * (1 - 1.0 / m.tau) ** 11.0 / 6.0
            ),
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule, 
        }
    )
    #print(we.m.thermal_conductivity.expr)
    we.write()


if __name__ == "__main__":
    main()

########################################################################################
#                                                                                      #
# R1234ZE EOS Expressions and Parameters:                                              #
#                                                                                      #
# Monika Thol and Eric W. Lemmon. "Equation of State for the Thermodynamic             #
#   Properties of trans-1,3,3,3-Tetrafluoropropene [R-1234ze(E)]." Int. J.             #
#   Thermophys, 37(3):1–16, 2016. doi:10.1007/s10765-016-2040-6.                       #
#                                                                                      #
# Richard A. Perkins and Marcia L. Huber. Measurement and Correlation of the Thermal   #
#   Conductivity of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf) and                        #
#   trans-1,3,3,3-Tetrafluoropropene (R1234ze(E)). J. Chem. Eng. Data, 56:4868–4874,   #
#   2011. doi:10.1021/je200811n.                                                       #
#                                                                                      #
# Huber ML, Assael MJ. Correlations for the Viscosity of 2,3,3,3-Tetrafluoroprop-1-ene # 
#   (R1234yf) and trans-1,3,3,3-Tetrafluoropropene (R1234ze(E)). Int J Refrig.         #
#   2016; 71:39-45. doi:10.1016/j.ijrefrig.2016.08.007                                 #
#                                                                                      #
########################################################################################

import pyomo.environ as pyo
from helmholtz_parameters import WriteParameters

def thermal_conductivity_rule(m):
    a = [
        -0.0103589,
        0.0308929,
        0.000230348,
    ]
    b = {
        (1, 1): -0.0428296,
        (1, 2): 0.0434288,
        (2, 1): 0.0927099,
        (2, 2): -0.0605844,
        (3, 1): -0.0702107,
        (3, 2): 0.0440187,
        (4, 1): 0.0249708,
        (4, 2): -0.0155082,
        (5, 1): -0.00301838,
        (5, 2): 0.00210190,
    }
    T = m.T_star / m.tau
    rho = m.delta * m.rho_star 
    Tred = 1.0 / m.tau
    l0 = sum(ai * Tred**i for i, ai in enumerate(a))
    lr = sum(m.delta**i * (b[i, 1] + b[i, 2] * Tred) for i in range(1, 6))
    lc = 0
    return 1000*(l0 + lr + lc)

def viscosity_rule(m):
    a = [
        -963382,
        9614.09,
        -13.233,
        0.0360562,
        122059,
        -224.741,
    ]
    b = [
        -19.572881,
        219.73999,
        -1015.3226,
        2471.0125,
        -3375.1717,
        2491.6597,
        -787.26086,
        14.085455,
        -0.34664158,
    ]
    c = [
        8.61691913,
        0,
        20.83024738,
        0,
        0.54243690,
        -10.49684841,
        -1.38137689,
        1,
        0,
    ]
    T = m.T_star / m.tau 
    rho = m.delta * m.rho_star / m.MW * 1000
    eok = 340
    sigma = 5.017e-1
    M = 114.0415928
    NA = 6.0221408e23
    Ts = T / eok
    eta0 = (a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3) / (
        a[4] + a[5] * T + T**2
    )
    Bs = (
        sum(b[i] * Ts ** (-0.25 * i) for i in range(0, 7))
        + b[7] * Ts**-2.5
        + b[8] * Ts**-5.5
    )
    B = NA * Bs * sigma**3 / 1e9**3
    eta1 = eta0 * B
    Tr = 1 / m.tau
    etar = (
        m.delta ** (2.0 / 3.0)
        * Tr**0.5
        * (
            c[0]
            + c[1] * m.delta
            + c[2] * m.delta**2
            + (
                c[3] * m.delta
                + c[4] * m.delta**6
                + c[5] * m.delta * Tr**2
                + c[6] * m.delta**5 * Tr
            )
            / (c[7] * Tr + c[8] * m.delta * Tr)
        )
    )
    return (eta0 + eta1 * rho + etar)


def main():
    m = pyo.ConcreteModel()
    m.delta = pyo.Var()
    m.tau = pyo.Var()

    c = {
        6: 2,
        7: 2,
        8: 1,
        9: 2,
        10: 1,
    }

    d = {
        1: 4,
        2: 1,
        3: 1,
        4: 2,
        5: 3,
        6: 1,
        7: 3,
        8: 2,
        9: 2,
        10: 7,
        11: 1,
        12: 1,
        13: 3,
        14: 3,
        15: 2,
        16: 1,
    }

    t = {
        1: 1.0,
        2: 0.223,
        3: 0.755,
        4: 1.24,
        5: 0.44,
        6: 2.0,
        7: 2.2,
        8: 1.2,
        9: 1.5,
        10: 0.9,
        11: 1.33,
        12: 1.75,
        13: 2.11,
        14: 1.0,
        15: 1.5,
        16: 1.0,
    }

    n = {
        1: 0.03982797,
        2: 1.812227,
        3: -2.537512,
        4: -0.5333254,
        5: 0.1677031,
        6: -1.323801,
        7: -0.6694654,
        8: 0.8072718,
        9: -0.7740229,
        10: -0.01843846,
        11: 1.407916,
        12: -0.4237082,
        13: -0.2270068,
        14: -0.805213,
        15: 0.00994318,
        16: -0.008798793,
    }

    a = {
        11: 1.0,
        12: 1.61,
        13: 1.24,
        14: 9.34,
        15: 5.78,
        16: 3.08,
    }

    b = {
        11: 1.21,
        12: 1.37,
        13: 0.98,
        14: 171.0,
        15: 47.4,
        16: 15.4,
    }

    g = {
        11: 0.943,
        12: 0.642,
        13: 0.59,
        14: 1.2,
        15: 1.33,
        16: 0.64,
    }

    e = {
        11: 0.728,
        12: 0.87,
        13: 0.855,
        14: 0.79,
        15: 1.3,
        16: 0.71,
    }

    n0 = {
        1: -12.558347537,  # aka a0[1]
        2: 8.7912297624,  # aka a0[2]
        3: 4.0 - 1,  # aka c0 - 1
        4: 9.3575,  # aka v[1]
        5: 10.717,  # aka v[2]
    }

    g0 = {  # aka u[k] / Tc
        4: 513 / 382.513,
        5: 1972 / 382.513,
    }

    we = WriteParameters(
        comp="r1234ze",
        R=0.07290727,
        MW=114.0416,
        T_star=382.513,
        rho_star=489.238,
        Tc=382.513,
        rhoc=489.238,
        Pc=3634.86,
        Tt=169.0,
        Pt=0.228564,
        rhot_l=1510.76,
        rhot_v=0.0185566,
        P_min=1e-9,
        P_max=5e5,
        rho_max=2500.0,
        T_min=150,
        T_max=1300,
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
                for i in range(6, 10 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                for i in range(11, 16 + 1)
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
                for i in range(6, 10 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (d[i] / m.delta - 2 * a[i] * (m.delta - e[i]))
                for i in range(11, 16 + 1)
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
                for i in range(6, 10 + 1)
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
                for i in range(11, 16 + 1)
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
                for i in range(6, 10 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (t[i] / m.tau - 2 * b[i] * (m.tau - g[i]))
                for i in range(11, 16 + 1)
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
                for i in range(6, 10 + 1)
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
                for i in range(11, 16 + 1)
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
                for i in range(6, 10 + 1)
            )
            + sum(
                n[i]
                * m.delta ** d[i]
                * m.tau ** t[i]
                * pyo.exp(-a[i] * (m.delta - e[i]) ** 2 - b[i] * (m.tau - g[i]) ** 2)
                * (d[i] / m.delta - 2.0 * a[i] * (m.delta - e[i]))
                * (t[i] / m.tau - 2.0 * b[i] * (m.tau - g[i]))
                for i in range(11, 16 + 1)
            ),
            "delta_v_sat_approx": pyo.exp(
                -1.0308 * (1 - 1 / m.tau) ** 0.24
                + -5.0422 * (1 - 1 / m.tau) ** 0.72
                + -11.5 * (1 - 1 / m.tau) ** 2.1
                + -37.499 * (1 - 1 / m.tau) ** 4.8
                + -77.945 * (1 - 1 / m.tau) ** 9.5
            ),
            "delta_l_sat_approx": (
                1.00
                + 1.1913 * (1 - 1 / m.tau) ** 0.27
                + 2.2456 * (1 - 1 / m.tau) ** 0.7
                - 1.7747 * (1 - 1 / m.tau) ** 1.25
                + 1.3096 * (1 - 1 / m.tau) ** 1.9
            ),
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule, 
        }
    )
    we.write()


if __name__ == "__main__":
    main()

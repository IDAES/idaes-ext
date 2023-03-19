######################################################################################
#                                                                                    #
# R134A EOS Expressions and Parameters:                                              #
#                                                                                    #
# Tillner-Roth, R.; Baehr, H.D., An International Standard Formulation for the       #
#    Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane (HFC-134a) for            #
#    Temperatures from 170 K to 455 K and Pressures up to 70 MPa, J. Phys. Chem.     #
#    Ref. Data, 1994, 23, 5, 657-729, https://doi.org/10.1063/1.555958               #
#                                                                                    #
# Perkins, R.A.; Laesecke, A.; Howley, J.; Ramires, M.L.V.; Gurova, A.N.; Cusco, L., # 
#    Experimental thermal conductivity values for the IUPAC round-robin sample of    #
#    1,1,1,2-tetrafluoroethane (R134a), NIST Interagency/Internal Report (NISTIR)    #
#    - 6605, 2000, https://doi.org/10.6028/NIST.IR.6605.                             #
#                                                                                    #
# Huber, M.L.; Laesecke, A.; Perkins, R.A., Model for the Viscosity and Thermal      #
#    Conductivity of Refrigerants, Including a New Correlation for the Viscosity     #
#    of R134a, Ind. Eng. Chem. Res., 2003, 42, 13, 3163-3178,                        #
#    https://doi.org/10.1021/ie0300880.                                              #
#                                                                                    #
# Thermal conductivity parameter errata correction from CoolProp parameter file.     #
# lambda^d.g. a1: 8.00982 -> 8.00982e-5                                              #
#                                                                                    #
######################################################################################

import math
import pyomo.environ as pyo
from idaes.core.util.math import smooth_max
from helmholtz_parameters import WriteParameters

def thermal_conductivity_rule(m):
    a = [
        -1.05248e-2,
        8.00982e-5,
    ]
    b = {
        1: 1.836526,
        2: 5.126143,
        3: -1.436883,
        4: 6.261441e-1,
    }

    T = m.T_star / m.tau
    MW = 0.102032 # kg/mol
    rho_star = m.rho_star / MW
    rho = m.delta * rho_star # mol/m^3
    tau = m.tau 
    delta = m.delta
    rho_crit = 5049.886  # mol/m^3
    lamb_red = 2.055e-3 # reducing tc W/m/K
    Tref = 561.411 # reference T K
    k = 1.380649e-23
    Pc = 4.05928 # MPa    
    big_gam = 0.0496 
    R0 = 1.03
    gamma = 1.239
    qd = 1892020000.0
    xi0 = 1.94e-10 # m
    nu = 0.63
    
    m.cp = pyo.ExternalFunction(library="", function="cp")
    m.cv = pyo.ExternalFunction(library="", function="cv")
    m.mu = pyo.ExternalFunction(library="", function="mu")
    m.itc = pyo.ExternalFunction(library="", function="itc")

    drho_dp = m.itc("r134a", delta, m.T_star/T)*rho_star*m.delta
    drho_dp_ref = m.itc("r134a", delta, m.T_star/Tref)*rho_star*m.delta

    deltchi = smooth_max(Pc * rho / rho_crit**2 * (drho_dp - drho_dp_ref * Tref/T), 0, eps=1e-8)
    xi = xi0*(deltchi/big_gam)**(nu/gamma)

    cp = m.cp("r134a", delta, tau)
    cv = m.cv("r134a", delta, tau)
    mu = m.mu("r134a", delta, tau) / 1e6
    y = qd * xi
    kappa_inv = cv/cp
    pi = math.pi
    Omega = 2.0 / pi * ((1.0 - kappa_inv) * pyo.atan(y) + kappa_inv * y)
    Omega0 = 2.0 / pi * (1.0 - pyo.exp(-1.0 / (1.0 / y + 1.0 / 3.0 * (y * rho_crit / rho)**2)))
    lambda_c = 1000*MW*cp*rho*R0*k*T/6.0/math.pi/mu/xi*(Omega - Omega0)
    lambda_dg = a[0] + a[1] * T
    lambda_r = lamb_red*sum(bi * (rho/rho_crit)**i for i, bi in b.items())
    return 1000*(lambda_dg + lambda_r + lambda_c)

def viscosity_rule(m):
    a = [
        0.355404,
        -0.464337,
        0.257353e-1,
    ]
    b = [
        -19.572881,
        219.73999,
        -1015.3226,
        2471.01251,
        -3375.1717,
        2491.6597,
        -787.26086,
        14.085455,
        -0.34664158,
    ]
    te = [
        0.0,
        -0.25,
        -0.50,
        -0.75,
        -1.00,
        -1.25,
        -1.50,
        -2.50,
        -5.50,
    ]
    c = {
        1: -20.6900719,
        2: 0.356029549,
        3: 2.11101816,
        4: 13.9601415,
        5: -4.5643502,
        6: -3.51593275,
        7: 214.76332,
        8: -0.890173375e-1,
        9: 0.100035295,
        10: 3.163695636,
    }
    T = m.T_star / m.tau
    rho = m.delta * m.rho_star / m.MW * 1000
    M = 102.031
    sigma = 0.46893
    eok = 299.363
    NA = 6.0221408e23
    Ts = T / eok
    vs = pyo.exp(sum(ai * pyo.log(Ts) ** i for i, ai in enumerate(a)))
    Bs = sum(bi * Ts**ti for bi, ti in zip(b, te))
    B = NA * sigma**3 * Bs / 1e9**3
    etas = 0.021357 * pyo.sqrt(M * T) / (sigma**2 * vs)
    tau = T / m.Tc
    delta0 = c[10] / (1 + c[8] * tau + c[9] * tau**2)
    delta = rho / 5017.053
    eta = (
        c[1] * delta
        + (c[2] / tau**6 + c[3] / tau**2 + c[4] / tau**0.5 + c[5] * tau**2)
        * delta**2
        + c[6] * delta**3
        + c[7] / (delta0 - delta)
        - c[7] / delta0
    )
    return (etas * (1 + B * rho) + eta)

def surface_tension_rule(m):
    x = smooth_max(1 - 1/m.tau*m.T_star/m.Tc, 0, 1e-8)
    sig = [0.05801]
    ns = [1.241]
    return sum(s*x**n for s, n in zip(sig, ns))

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
        2: 0.498223,
        3: 0.2458698e-1,
        4: 0.8570145e-3,
        5: 0.4788584e-3,
        6: -0.1800808e1,
        7: 0.2671641,
        8: -0.4781652e-1,
        9: 0.1423987e-1,
        10: 0.3324062,
        11: -0.7485907e-2,
        12: 0.1017263e-3,
        13: -0.5184567,
        14: -0.8692288e-1,
        15: 0.2057144,
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
        R=8.14888564372e-2,
        MW=102.032,
        T_star=374.18,
        rho_star=508.0,
        Tc=374.21,
        rhoc=511.89995169599996,
        Pc=4059.28,
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
            / m.rho_star
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
            / m.rho_star,
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule,
            "surface_tension": surface_tension_rule,
        }
    )
    we.write()


if __name__ == "__main__":
    main()

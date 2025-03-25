import os
import math
import pytest
import parameterized
import itertools
import scipy.sparse as sparse
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
import pyomo.common.unittest as unittest
from pyomo.common.collections import ComponentMap
from pyomo.repn.util import FileDeterminism
from pyomo.util.subsystems import create_subsystem_block
from pyomo.contrib.sensitivity_toolbox.sens import sensitivity_calculation, _add_sensitivity_suffixes
from pyomo.contrib.sensitivity_toolbox.k_aug import InTempDir
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from contextlib import nullcontext


"""In addition to the above Python dependencies, the following environment must be
set up to run these tests:
- IDAES binaries have been downloaded to $HOME/.idaes/bin
- This directory should be on: PATH and (DY)LD_LIBRARY_PATH
- idaes-local-*.tar.gz is unpacked (somewhere)
- idaes-local/lib/pkgconfig is on PKG_CONFIG_PATH
- CyIpopt is installed with `pip install cyipopt`. If a wheel doesn't build, maybe
  it needs to be removed from the pip cache with `pip cache remove cyipopt`
- To test the Helmholtz EOS external functions, the IDAES_HELMHOLTZ_DATA_PATH environment
  variable needs to be set to $HOME/.idaes/bin/helm_data/ (Note that the trailing '/'
  appears to be necessary)

"""


# TODO: This directory will be wherever we download the binaries that we
# want to test, likely just in the current working directory.
if "IDAES_DIR" in os.environ:
    IDAES_DIR = os.environ["IDAES_DIR"]
else:
    # Note that this directory is specific to mac and linux
    IDAES_DIR = os.path.join(os.environ["HOME"], ".idaes")
ipopts_to_test = [
    ("ipopt", os.path.join(IDAES_DIR, "bin", "ipopt")),
    ("ipopt_l1", os.path.join(IDAES_DIR, "bin", "ipopt_l1")),
    ("cyipopt", None),
]
ipopt_options_to_test = [
    ("default", {}),
    ("mumps", {"print_user_options": "yes", "linear_solver": "mumps"}),
    ("ma27", {"print_user_options": "yes", "linear_solver": "ma27"}),
    ("ma57", {"print_user_options": "yes", "linear_solver": "ma57"}),
    ("ma57_metis", {"print_user_options": "yes", "linear_solver": "ma57", "ma57_pivot_order": 4}),
]
sensitivity_solvers = [
    ("ipopt", "k_aug", "dot_sens"),
    ("ipopt_sens", "ipopt_sens", None),
    ("ipopt_sens_l1", "ipopt_sense_l1", None),
]
TEE = True
ipopt_test_data = list(itertools.product(ipopts_to_test, ipopt_options_to_test))
ipopt_test_data = [
    (f"{ipoptname}_{optname}", ipoptname, ipoptexe, options)
    for (ipoptname, ipoptexe), (optname, options) in ipopt_test_data
]


def _test_ipopt_with_options(name, exe, options):
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1,2], initialize=1.5)
    m.con = pyo.Constraint(expr=m.x[1]*m.x[2] == 0.5)
    m.obj = pyo.Objective(expr=m.x[1]**2 + 2*m.x[2]**2)

    if exe is None:
        solver = pyo.SolverFactory(name, options=options)
    else:
        solver = pyo.SolverFactory(name, executable=exe, options=options)

    if "ipopt_l1" in name:
        # Run this in a temp dir so we don't pollute the working directory with
        # ipopt_l1's files. See https://github.com/IDAES/idaes-ext/issues/275
        context = InTempDir()
    else:
        context = nullcontext()
    with context:
        solver.solve(m, tee=TEE)

    target_sol = [("x[1]", 0.840896415), ("x[2]", 0.594603557)]
    assert all(
        math.isclose(m.find_component(name).value, val, abs_tol=1e-7)
        for name, val in target_sol
    )


class TestIpopt(unittest.TestCase):

    @parameterized.parameterized.expand(ipopt_test_data)
    def test_ipopt(self, test_name, solver_name, exe, options):
        _test_ipopt_with_options(solver_name, exe, options)


class TestBonmin:

    exe = os.path.join(IDAES_DIR, "bin", "bonmin")

    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1.5)
        m.y = pyo.Var(domain=pyo.PositiveIntegers)
        m.con = pyo.Constraint(expr=m.x[1]*m.x[2] == m.y)
        m.obj = pyo.Objective(expr=m.x[1]**2 + 2*m.x[2]**2)
        return m

    def test_bonmin_default(self):
        m = self._make_model()
        solver = pyo.SolverFactory("bonmin", executable=self.exe)
        solver.solve(m, tee=TEE)

        assert math.isclose(m.y.value, 1.0, abs_tol=1e-7)
        assert math.isclose(m.x[1].value, 1.18920710, abs_tol=1e-7)
        assert math.isclose(m.x[2].value, 0.84089641, abs_tol=1e-7)


class TestCouenne:

    exe = os.path.join(IDAES_DIR, "bin", "couenne")

    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1.5)
        m.y = pyo.Var(domain=pyo.PositiveIntegers)
        m.con = pyo.Constraint(expr=m.x[1]*m.x[2] == m.y)
        m.obj = pyo.Objective(expr=(m.x[1] + 0.01)**2 + 2*(m.x[2] + 0.01)**2)
        return m

    def test_couenne_default(self):
        m = self._make_model()
        solver = pyo.SolverFactory("couenne", executable=self.exe)
        solver.solve(m, tee=TEE)

        assert math.isclose(m.y.value, 1.0, abs_tol=1e-7)
        assert math.isclose(m.x[1].value, -1.18816674, abs_tol=1e-7)
        assert math.isclose(m.x[2].value, -0.84163270, abs_tol=1e-7)


def _test_sensitivity(
    solver_name,
    solver_exe,
    sens_name,
    sens_exe,
    update_name,
    update_exe,
):
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1,2], initialize=1.5)
    m.p = pyo.Param(mutable=True, initialize=0.5)
    m.con = pyo.Constraint(expr=m.x[1]*m.x[2] == m.p)
    m.obj = pyo.Objective(expr=m.x[1]**2 + 2*m.x[2]**2)

    if solver_exe is None:
        solver = pyo.SolverFactory(solver_name)
    else:
        solver = pyo.SolverFactory(solver_name, executable=solver_exe)
    solver.solve(m, tee=TEE, keepfiles=True)
    if sens_name == "k_aug":
        sensitivity_executable = (sens_exe, update_exe)
    else:
        sensitivity_executable = sens_exe
    sensitivity_calculation(
        sens_name,
        m,
        [m.p],
        [0.7],
        cloneModel=False,
        tee=TEE,
        #sensitivity_executable=sensitivity_executable,
        #solver_executable=solver_exe,
    )
    solution = {"x[1]": 0.95, "x[2]": 0.75}
    if sens_name == "sipopt":
        # sipopt puts the perturbed solution in suffixes
        for var, val in solution.items():
            # Use a loose tolerance because methods seem to give different solutions...
            assert math.isclose(
                m.sens_sol_state_1[m.find_component(var)], val, abs_tol=1e-1
            )
    elif sens_name == "k_aug":
        # K_aug puts the perturbed solution back in the model
        for var, val in solution.items():
            # Use a loose tolerance because methods seem to give different solutions...
            assert math.isclose(m.find_component(var).value, val, abs_tol=1e-1)


class TestSensitivity:

    ipopt_exe = os.path.join(IDAES_DIR, "bin", "ipopt")
    sipopt_exe = os.path.join(IDAES_DIR, "bin", "ipopt_sens")
    k_aug_exe = os.path.join(IDAES_DIR, "bin", "k_aug")
    dot_sens_exe = os.path.join(IDAES_DIR, "bin", "dot_sens")

    def test_k_aug(self):
        _test_sensitivity(
            "ipopt",
            self.ipopt_exe,
            "k_aug",
            self.k_aug_exe,
            "dot_sens",
            self.dot_sens_exe,
        )

    def test_sipopt(self):
        _test_sensitivity(
            "ipopt", self.ipopt_exe, "sipopt", self.sipopt_exe, None, None
        )


class TestPyNumeroASL:

    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=0.0, bounds=(0, 10))
        m.ineq = pyo.Constraint(pyo.PositiveIntegers)
        m.ineq[1] = m.x[1] + 3*m.x[2] <= 7
        m.ineq[2] = 5*m.x[1] + m.x[2] <= 10
        m.obj = pyo.Objective(expr=m.x[1] + m.x[2], sense=pyo.maximize)
        return m

    def test_asl(self):
        m = self._make_model()
        # Note that PyomoNLP relies on Pyomo being able to find libpynumero_ASL.
        # Here, these should be on LD_LIBRARY_PATH or DYLD_LIBRARY_PATH.
        nlp = PyomoNLP(
            # Explicitly sort variables and constraints so we know where in Jacobian
            # each coefficient will be.
            m, nl_file_options=dict(file_determinism=FileDeterminism.SORT_SYMBOLS)
        )
        jac = nlp.evaluate_jacobian()
        row = [0, 0, 1, 1]
        col = [0, 1, 0, 1]
        data = [1.0, 3.0, 5.0, 1.0]
        expected = set((r, c, d) for r, c, d in zip(row, col, data))
        actual = set((r, c, d) for r, c, d in zip(jac.row, jac.col, jac.data))
        assert expected == actual


class TestCLP:

    exe = os.path.join(IDAES_DIR, "bin", "clp")

    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=0.0, bounds=(0, 10))
        m.ineq = pyo.Constraint(pyo.PositiveIntegers)
        m.ineq[1] = m.x[1] + 3*m.x[2] <= 7
        m.ineq[2] = 5*m.x[1] + m.x[2] <= 10
        m.obj = pyo.Objective(expr=m.x[1] + m.x[2], sense=pyo.maximize)
        return m

    def test_clp(self):
        solution = {"x[1]": 7-15/2.8, "x[2]": 5/2.8}
        m = self._make_model()
        solver = pyo.SolverFactory("clp", executable=self.exe)
        solver.solve(m, tee=TEE)
        for name, val in solution.items():
            assert math.isclose(m.find_component(name).value, val, abs_tol=1e-5)


class TestCBC:

    exe = os.path.join(IDAES_DIR, "bin", "cbc")

    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=0.0, bounds=(0, 10))
        m.y = pyo.Var(domain=pyo.Binary)
        m.ineq = pyo.Constraint(pyo.PositiveIntegers)
        m.ineq[1] = m.x[1] + 3*m.x[2] <= 7
        m.ineq[2] = 5*m.x[1] + m.x[2] <= 10
        m.ineq[3] = m.y + 1 >= m.x[1]
        m.obj = pyo.Objective(expr=m.x[1] + m.x[2], sense=pyo.maximize)
        return m

    def test_cbc(self):
        solution = {"x[1]": 7-15/2.8, "x[2]": 5/2.8, "y": 1.0}
        m = self._make_model()
        solver = pyo.SolverFactory("cbc", executable=self.exe)
        solver.solve(m, tee=TEE)
        for name, val in solution.items():
            assert math.isclose(m.find_component(name).value, val, abs_tol=1e-5)


class TestExternalFunctions:

    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1.0, bounds=(0, 10))
        m.cubic_root_low = pyo.ExternalFunction(
            library="cubic_roots", function="cubic_root_l"
        )
        m.x_over_exp_x_minus_one = pyo.ExternalFunction(
            library="functions", function="x_over_exp_x_minus_one"
        )
        m.isothermal_compressibility_vap_up = pyo.ExternalFunction(
            library="general_helmholtz_external",
            function="isothermal_compressibility_vap_up",
        )
        m.eq = pyo.Constraint(pyo.PositiveIntegers)
        m.eq[1] = m.x[1] == m.cubic_root_low(-1, 2, m.x[2])
        m.obj = pyo.Objective(expr=m.x[1]**2 + m.x[2]**2)
        return m

    def test_cubic_roots(self):
        m = self._make_model()
        cbrt = pyo.value(m.cubic_root_low(-1.0, 2.0, -3.0))
        assert math.isclose(cbrt**3 - cbrt**2 + 2*cbrt - 3.0, 0.0, abs_tol=1e-8)

    def test_functions(self):
        m = self._make_model()
        x = 1.2
        y = pyo.value(m.x_over_exp_x_minus_one(x))
        assert math.isclose(y, x / (math.exp(x) - 1.0), abs_tol=1e-8)

    def test_helmholtz(self):
        m = self._make_model()
        # delta is P/Pc ration; tau is T/Tc ratio
        delta = 2.0
        tau = 3.0
        # This is a little more difficult to test...
        # I just chose some random point, and now am asserting that we keep getting
        # the same answer. This is not based on any first-principles calculation.
        z = pyo.value(m.isothermal_compressibility_vap_up("co2", delta, tau))
        assert math.isclose(z, 333.47, abs_tol=0.1)


class TestPetsc:

    exe = os.path.join(IDAES_DIR, "bin", "petsc")

    def test_petsc_snes(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1.0)
        m.eq = pyo.Constraint(pyo.PositiveIntegers)
        m.eq[1] = m.x[1] - 2*pyo.log(m.x[2]) == 3
        m.eq[2] = m.x[1] * m.x[2]**1.5 == 1

        solver = pyo.SolverFactory("petsc", executable=self.exe)
        res = solver.solve(m, tee=True)
        pyo.assert_optimal_termination(res)
        # I just pulled these numbers from a successful solve.
        assert math.isclose(m.x[1].value, 2.045, abs_tol=0.01)
        assert math.isclose(m.x[2].value, 0.621, abs_tol=0.01)
        # Just checking that we actually ran SNES here (and not some
        # other petsc solver)
        assert res.solver.message == "SNES_CONVERGED_FNORM_ABS"

    def test_petsc_ts(self):
        m = pyo.ConcreteModel()
        tf = 1.0
        m.time = ContinuousSet(bounds=[0.0, tf])
        m.direction = pyo.Set(initialize=["in", "out"])
        m.area = pyo.Param(initialize=2.5)
        m.flow_coef = pyo.Param(initialize=4.0)
        m.height = pyo.Var(m.time, initialize=1.0)
        m.flow = pyo.Var(m.time, m.direction, initialize=1.0)
        m.dheight_dt = DerivativeVar(m.height, wrt=m.time, initialize=0.0)

        @m.Constraint(m.time)
        def flow_out_eqn(m, t):
            return m.flow[t, "out"] == m.flow_coef * pyo.sqrt(m.height[t])

        @m.Constraint(m.time)
        def height_diff_eqn(m, t):
            return m.area * m.dheight_dt[t] == m.flow[t, "in"] - m.flow[t, "out"]

        # Fix initial condition
        m.height[0].fix()
        # Fix input variable
        m.flow[:, "in"].fix()
        for t in m.time:
            if t < tf / 2:
                m.flow[t, "in"].fix(1.0)
            else:
                m.flow[t, "in"].fix(1.2)

        ALGEBRAIC = 0
        DIFFERENTIAL = 1
        DERIVATIVE = 2

        m.dae_suffix = pyo.Suffix(
            direction=pyo.Suffix.IMPORT_EXPORT,
            datatype=pyo.Suffix.INT,
        )
        m.dae_link = pyo.Suffix(
            direction=pyo.Suffix.IMPORT_EXPORT,
            datatype=pyo.Suffix.INT,
        )

        diffvar_index = 0
        t = tf
        m.dae_suffix[m.height[t]] = DIFFERENTIAL
        m.dae_suffix[m.flow[t, "out"]] = ALGEBRAIC
        m.dae_suffix[m.dheight_dt[t]] = DERIVATIVE
        m.dae_link[m.height[t]] = 0
        m.dae_link[m.dheight_dt[t]] = 0

        options = {"--dae_solve": "", "--ts_init_time": 0.0, "--ts_max_time": tf}
        solver = pyo.SolverFactory("petsc", executable=self.exe, options=options)

        res = solver.solve(m, tee=True)
        pyo.assert_optimal_termination(res)
        m.height.pprint()

        assert math.isclose(m.height[tf].value, 0.43, abs_tol=0.01)


if __name__ == "__main__":
    pytest.main([__file__])

# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver
import copy

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestNewSubroutines(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl", debug=False)
        self.avl_solver.add_constraint("alpha", 25.0)

    def test_residual(self):
        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.VRTX_R.RES[:]
        rhs = self.avl_solver.avl.VRTX_R.RHS[:]
        np.testing.assert_allclose(
            res,
            -1 * rhs,
            atol=1e-15,
        )

        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.VRTX_R.RES[:]
        res_d = self.avl_solver.avl.VRTX_R.RES_D[:]

        np.testing.assert_allclose(
            res,
            np.zeros_like(res),
            atol=5e-14,
        )

        np.testing.assert_allclose(
            res_d,
            np.zeros_like(res_d),
            atol=5e-14,
        )

    def test_new_solve(self):
        self.avl_solver.add_constraint("Elevator", 10.00)
        self.avl_solver.add_constraint("alpha", 10.00)
        self.avl_solver.add_constraint("beta", 10.00)
        
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.velsum()
        wv_new = self.avl_solver.avl.SOLV_R.WV
        self.avl_solver.avl.aero()
        gam_new = self.avl_solver.avl.VRTX_R.GAM[:]
        coef_data_new = self.avl_solver.get_case_total_data()
        coef_derivs_new = self.avl_solver.get_case_coef_derivs()

        self.avl_solver.execute_run()
        gam = self.avl_solver.avl.VRTX_R.GAM[:]
        wv = self.avl_solver.avl.SOLV_R.WV
        coef_data = self.avl_solver.get_case_total_data()
        coef_derivs = self.avl_solver.get_case_coef_derivs()

        np.testing.assert_allclose(
            wv,
            wv_new,
            atol=1e-15,
        )
        np.testing.assert_allclose(
            gam,
            gam_new,
            atol=1e-15,
        )
        for func_key in coef_data:
            np.testing.assert_allclose(
                coef_data[func_key],
                coef_data_new[func_key],
                atol=1e-14,
                err_msg=f"func_key {func_key}",
            )

        for func_key in coef_derivs:
            for consurf_key in coef_derivs[func_key]:
                np.testing.assert_allclose(
                    coef_derivs[func_key][consurf_key],
                    coef_derivs_new[func_key][consurf_key],
                    err_msg=f"deriv of func_key {func_key} wrt {consurf_key}",
                    atol=1e-14,
                )
                

if __name__ == "__main__":
    unittest.main()

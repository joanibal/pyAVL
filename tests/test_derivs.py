# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver

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


class TestPartials(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_fwd_cltot_alpha_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.case_r_d.alfad = 1.0
        self.avl_solver.avl.aero_d()
        dcl_dalfa = self.avl_solver.avl.case_r_d.cltotd

        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-8
        alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)

        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dalfa_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h

        np.testing.assert_allclose(
            dcl_dalfa,
            dcl_dalfa_fd,
            rtol=1e-6,
        )

    def test_rev_cltot_alpha_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.case_r_d.cltotd = 1.0
        self.avl_solver.avl.aero_b()
        dcl_dalfa = self.avl_solver.avl.case_r_d.alfad

        print(dcl_dalfa)
        # coef_data = self.avl_solver.get_case_total_data()

        # # use finite difference
        # h = 1e-8
        # alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        # self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)

        # self.avl_solver.avl.aero()
        # coef_data_peturb = self.avl_solver.get_case_total_data()
        # dcl_dalfa_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h

        # np.testing.assert_allclose(
        #     dcl_dalfa,
        #     dcl_dalfa_fd,
        #     rtol=1e-6,
        # )

    def test_fwd_cltot_gam_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.vrtx_r_d.gamd[1] = 1.0
        self.avl_solver.avl.aero_d()
        dcl_dgam = self.avl_solver.avl.case_r_d.cltotd

        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-1

        gam = self.avl_solver.get_avl_fort_var("VRTX_R", "GAM")
        gam[1] = gam[1] + h
        self.avl_solver.set_avl_fort_var("VRTX_R", "GAM", gam)

        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dgam_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h

        np.testing.assert_allclose(
            dcl_dgam,
            dcl_dgam_fd,
            rtol=1e-14,
        )

    # there is limited need to test this
    # def test_jacobian(self):
    #     """Test the Jacobian $\partial{res}/\partial{gam} of ADVL"""
    #     self.avl_solver.avl.get_res()


class TestTotals(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver.add_constraint("alpha", 5.0)
        self.avl_solver.add_constraint("beta", 9.0)
        self.avl_solver.add_constraint("roll rate", 1.2)
        self.avl_solver.add_constraint("pitch rate", 0.1)
        self.avl_solver.add_constraint("yaw rate", 0.8)

        # base line CL

    def test_new_solve(self):

        self.avl_solver.execute_run()
        gam = self.avl_solver.avl.vrtx_r.gam[:]
        self.avl_solver.avl.exec_rhs()
        gam_new = self.avl_solver.avl.vrtx_r.gam[:]

        np.testing.assert_allclose(
            gam,
            gam_new,
            atol=1e-15,
        )

    def test_residual(self):

        # self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.vrtx_r.res[:]
        rhs = self.avl_solver.avl.vrtx_r.rhs[:]

        np.testing.assert_allclose(
            res,
            -1 * rhs,
            atol=1e-15,
        )

        print(np.linalg.norm(rhs))
        print(np.linalg.norm(res))

        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.vrtx_r.res[:]
        print(res[:10])

        np.testing.assert_allclose(
            res,
            np.zeros_like(res),
            atol=5e-14,
        )

    # def test_cltot_alpha_derivs(self):
    #     pass


if __name__ == "__main__":
    unittest.main()

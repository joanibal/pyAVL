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


class TestFunctionPartials(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_fwd_cltot_alpha_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.CASE_R_diff.ALFA_diff = 1.0
        self.avl_solver.avl.aero_d()
        dcl_dalfa = self.avl_solver.avl.CASE_R_diff.CLTOT_diff

        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)

        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dalfa_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h

        np.testing.assert_allclose(
            dcl_dalfa,
            dcl_dalfa_fd,
            rtol=1e-5,
        )

    def test_rev_cltot_alpha_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.CASE_R_diff.CLTOT_diff = 1.0
        self.avl_solver.avl.aero_b()
        dcl_dalfa = self.avl_solver.avl.CASE_R_diff.ALFA_diff
        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)

        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dalfa_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h

        np.testing.assert_allclose(
            dcl_dalfa,
            dcl_dalfa_fd,
            rtol=1e-5,
        )

    def test_fwd_cltot_gam_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.VRTX_R_diff.GAM_diff[1] = 1.0
        self.avl_solver.avl.aero_d()
        dcl_dgam = self.avl_solver.avl.CASE_R_diff.CLTOT_diff

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

    def test_rev_cltot_gam_derivs(self):
        # base line CL
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        
        self.avl_solver.avl.CASE_R_diff.CLTOT_diff = 1
        self.avl_solver.avl.aero_b()
        dcl_dgam = self.avl_solver.avl.VRTX_R_diff.GAM_diff

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
            dcl_dgam[1],
            dcl_dgam_fd,
            rtol=1e-14,
        )

class TestResidualPartials(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_residual(self):

        # self.avl_solver.avl.exec_rhs()
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

        np.testing.assert_allclose(
            res,
            np.zeros_like(res),
            atol=5e-14,
        )
    
    def test_fwd_res_alpha_deriv(self):
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.CASE_R_diff.ALFA_diff = 1.0
        self.avl_solver.avl.get_res_d()
        dres_dalfa = self.avl_solver.avl.VRTX_R_diff.RES_diff
        
        res_0 = copy.deepcopy(self.avl_solver.avl.VRTX_R.RES[:])

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)
        self.avl_solver.avl.get_res()
        res_1 = self.avl_solver.avl.VRTX_R.RES[:]
        dres_dalfa_fd = (res_1 - res_0) / h

        np.testing.assert_allclose(
            dres_dalfa,
            dres_dalfa_fd,
            rtol=1e-7,
        )
        
    def test_rev_res_alpha_deriv(self):
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()

        self.avl_solver.avl.VRTX_R_diff.RES_diff[1] = 1.0
        self.avl_solver.avl.get_res_b()
        dres_dalfa = self.avl_solver.avl.CASE_R_diff.ALFA_diff
        
        self.avl_solver.avl.get_res()
        res_0 = self.avl_solver.avl.VRTX_R.RES[1]
        
        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)
        self.avl_solver.avl.get_res()
        res_1 = self.avl_solver.avl.VRTX_R.RES[1]
        dres_dalfa_fd = (res_1 - res_0) / h

        np.testing.assert_allclose(
            dres_dalfa,
            dres_dalfa_fd,
            rtol=1e-7,
        )
        
    # TODO: add dot product test
    # there is limited need to test this since it is linear (although it would be nice)
    # def test_jacobian(self):
    #     """Test the Jacobian $\partial{res}/\partial{gam} of ADVL"""
    #     self.avl_solver.avl.get_res()


class TestTotals(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver.add_constraint("alpha", 5.0)
        # self.avl_solver.add_constraint("beta", 9.0)
        # self.avl_solver.add_constraint("roll rate", 1.2)
        # self.avl_solver.add_constraint("pitch rate", 0.1)
        # self.avl_solver.add_constraint("yaw rate", 0.8)

        # base line CL

    def test_new_solve(self):

        self.avl_solver.execute_run()
        gam = self.avl_solver.avl.VRTX_R.GAM[:]
        self.avl_solver.avl.exec_rhs()
        gam_new = self.avl_solver.avl.VRTX_R.GAM[:]

        np.testing.assert_allclose(
            gam,
            gam_new,
            atol=1e-15,
        )
        
    def test_cltot_alpha_derivs(self):
        # get df_du and df_dx
        self.avl_solver.execute_run()
        
        # use reverse mode here
        self.avl_solver.avl.CASE_R_diff.CLTOT_diff = 1
        self.avl_solver.avl.aero_b()
        # dcl_dgam = self.avl_solver.avl.VRTX_R_diff.GAM_diff
        dcl_dalfa = self.avl_solver.avl.CASE_R_diff.ALFA_diff[()]

        # solve for the adjoint variable
        self.avl_solver.avl.solve_adjoint()
        dcl_dres = self.avl_solver.avl.VRTX_R_diff.RES_diff
        
        # combine adjoint with pr/px
        self.avl_solver.avl.get_res_b()
        dcl_dalfa += -1*self.avl_solver.avl.CASE_R_diff.ALFA_diff
        
        # use finite difference
        coef_data = self.avl_solver.get_case_total_data()
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_var("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_var("CASE_R", "ALFA", alpha + h)
        
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.aero()
        
        
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dalfa_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h
        
        np.testing.assert_allclose(
            dcl_dalfa,
            dcl_dalfa_fd,
            rtol=1e-8,
        )        


if __name__ == "__main__":
    unittest.main()

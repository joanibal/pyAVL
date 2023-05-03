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
        # self.avl_solver = AVLSolver(geo_file='aircraft_L1.avl')
        # print(self.avl_solver.avl.STRP_I.IJFRST[:20])
        # print(self.avl_solver.avl.SURF_I.IFRST[:10])
        self.avl_solver.add_constraint("alpha", 1.0)
        self.avl_solver.execute_run()
        

    def test_fwd_cltot_alpha_derivs(self):
        # base line CL

        self.avl_solver.avl.CASE_R_diff.ALFA_diff = 1.0
        self.avl_solver.avl.aero_d()
        dcl_dalfa = self.avl_solver.avl.CASE_R_diff.CLTOT_diff

        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_arr("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_arr("CASE_R", "ALFA", alpha + h)

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

        self.avl_solver.avl.CASE_R_diff.CLTOT_diff = 1.0
        self.avl_solver.avl.aero_b()
        dcl_dalfa = self.avl_solver.avl.CASE_R_diff.ALFA_diff
        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_arr("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_arr("CASE_R", "ALFA", alpha + h)

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

        self.avl_solver.avl.VRTX_R_diff.GAM_diff[1] = 1.0
        self.avl_solver.avl.aero_d()
        dcl_dgam = self.avl_solver.avl.CASE_R_diff.CLTOT_diff

        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-1

        gam = self.avl_solver.get_avl_fort_arr("VRTX_R", "GAM")
        gam[1] = gam[1] + h
        self.avl_solver.set_avl_fort_arr("VRTX_R", "GAM", gam)

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

        self.avl_solver.avl.CASE_R_diff.CLTOT_diff = 1
        self.avl_solver.avl.aero_b()
        dcl_dgam = self.avl_solver.avl.VRTX_R_diff.GAM_diff

        coef_data = self.avl_solver.get_case_total_data()

        # use finite difference
        h = 1e-1

        gam = self.avl_solver.get_avl_fort_arr("VRTX_R", "GAM")
        gam[1] = gam[1] + h
        self.avl_solver.set_avl_fort_arr("VRTX_R", "GAM", gam)

        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dgam_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h
        np.testing.assert_allclose(
            dcl_dgam[1],
            dcl_dgam_fd,
            rtol=1e-14,
        )

    def test_fwd_cltot_geom_derivs(self):
        # base line CL
        idx_surf = 0
        num_sec = self.avl_solver.avl.SURF_GEOM_I.NSEC[idx_surf]
        self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff[:num_sec, idx_surf] = 1.0
        self.avl_solver.avl.update_surfaces_d()
        self.avl_solver.avl.aero_d()
        dcl_daincs = self.avl_solver.avl.CASE_R_diff.CLTOT_diff
        # print("dcl_daincs", dcl_daincs)

        # # use finite difference
        coef_data = self.avl_solver.get_case_total_data()
        h = 1e-1
        self.avl_solver.avl.SURF_GEOM_R.AINCS[:num_sec, idx_surf] += h
        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_daincs_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h
        # print("dcl_daincs", dcl_daincs_fd)

        np.testing.assert_allclose(
            dcl_daincs,
            dcl_daincs_fd,
            rtol=1e-14,
        )

    def test_rev_cltot_geom_derivs(self):
        # base line CL
        idx_surf = 0
        idx_res = 1
        num_sec = self.avl_solver.avl.SURF_GEOM_I.NSEC[idx_surf]
        coef_data = self.avl_solver.get_case_total_data()
        
        # print('cl', coef_data["CL"])
        gam = self.avl_solver.avl.VRTX_R.GAM[:]
        # print('gam', np.linalg.norm(gam))


        self.avl_solver.avl.aero_b()
        self.avl_solver.avl.update_surfaces_b()

        dcl_daincs = self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff[:num_sec, idx_surf]

        # print("dcl_daincs", dcl_daincs)

        # # use finite difference
        self.avl_solver.avl.update_surfaces()


        h = 1e-1
        self.avl_solver.avl.SURF_GEOM_R.AINCS[:num_sec, idx_surf] += h
        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_daincs_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h
        # print(coef_data_peturb["CL"] , coef_data["CL"])
        # print("dcl_daincs_fd", dcl_daincs_fd)

        np.testing.assert_allclose(
            np.sum(dcl_daincs),
            dcl_daincs_fd,
            rtol=1e-14,
        )


class TestNewSubroutines(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 10.0)

    def test_residual(self):

        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.VRTX_R.RES[:]
        rhs = self.avl_solver.avl.VRTX_R.RHS[:]
        # print("res", np.linalg.norm(res))
        # print("rhs", np.linalg.norm(rhs))
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


class TestResidualPartials(unittest.TestCase):
    # TODO: add dot product tests
    # TODO: add test of jacobian dr/dgam
    # there is limited need to test this since it is linear (although it would be nice)
    # def test_jacobian(self):
    #     """Test the Jacobian $\partial{res}/\partial{gam} of ADVL"""
    #     self.avl_solver.avl.get_res()

    def setUp(self):
        self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 10.0)
        self.avl_solver.execute_run()
        # self.avl_solver.avl.exec_rhs()

    def test_fwd_res_alpha_deriv(self):

        self.avl_solver.avl.CASE_R_diff.ALFA_diff = 1.0
        self.avl_solver.avl.get_res_d()
        dres_dalfa = self.avl_solver.avl.VRTX_R_diff.RES_diff

        res_0 = copy.deepcopy(self.avl_solver.avl.VRTX_R.RES[:])

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_arr("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_arr("CASE_R", "ALFA", alpha + h)
        self.avl_solver.avl.get_res()
        res_1 = self.avl_solver.avl.VRTX_R.RES[:]
        dres_dalfa_fd = (res_1 - res_0) / h

        np.testing.assert_allclose(
            dres_dalfa,
            dres_dalfa_fd,
            rtol=1e-7,
        )

    def test_rev_res_alpha_deriv(self):

        self.avl_solver.avl.VRTX_R_diff.RES_diff[1] = 1.0
        self.avl_solver.avl.get_res_b()
        dres_dalfa = self.avl_solver.avl.CASE_R_diff.ALFA_diff

        self.avl_solver.avl.get_res()
        res_0 = self.avl_solver.avl.VRTX_R.RES[1]

        # use finite difference
        h = 1e-7
        alpha = self.avl_solver.get_avl_fort_arr("CASE_R", "ALFA")
        self.avl_solver.set_avl_fort_arr("CASE_R", "ALFA", alpha + h)
        self.avl_solver.avl.get_res()
        res_1 = self.avl_solver.avl.VRTX_R.RES[1]
        dres_dalfa_fd = (res_1 - res_0) / h

        np.testing.assert_allclose(
            dres_dalfa,
            dres_dalfa_fd,
            rtol=1e-7,
        )

    def test_fwd_res_geom_derivs(self):
        # base line CL
        idx_surf = 0
        num_sec = self.avl_solver.avl.SURF_GEOM_I.NSEC[idx_surf]
        self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff[:num_sec, idx_surf] = 1.0
        self.avl_solver.avl.update_surfaces_d()
        enc_diff = self.avl_solver.avl.VRTX_R_diff.ENC_diff

        self.avl_solver.avl.get_res_d()
        dres_aincs = self.avl_solver.avl.VRTX_R_diff.RES_diff

        # # use finite difference
        res_0 = copy.deepcopy(self.avl_solver.avl.VRTX_R.RES[:])
        h = 1e-8
        self.avl_solver.avl.SURF_GEOM_R.AINCS[:num_sec, idx_surf] += h
        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        res_1 = copy.deepcopy(self.avl_solver.avl.VRTX_R.RES[:])

        dres_aincs_fd = (res_1 - res_0) / h

        np.testing.assert_allclose(
            dres_aincs,
            dres_aincs_fd,
            rtol=1e-7,
        )

    def test_rev_res_geom_derivs(self):
        idx_surf = 0
        idx_res = 0
        idx_sec = 1
        num_sec = self.avl_solver.avl.SURF_GEOM_I.NSEC[idx_surf]

        self.avl_solver.avl.VRTX_R_diff.RES_diff[idx_res] = 1.0
        self.avl_solver.avl.get_res_b()
        self.avl_solver.avl.update_surfaces_b()

        dres_aincs = self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff[idx_sec, idx_surf]
        denc_aincs = self.avl_solver.avl.VRTX_R_diff.ENC_diff[:, idx_res]
        

        self.avl_solver.avl.get_res()
        res_0 = copy.deepcopy(self.avl_solver.avl.VRTX_R.RES[idx_res])
        enc_0 = copy.deepcopy(self.avl_solver.avl.VRTX_R.ENC[:, idx_res])
        # # use finite difference
        h = 1e-8
        self.avl_solver.avl.SURF_GEOM_R.AINCS[idx_sec, idx_surf] += h
        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        res_1 = copy.deepcopy(self.avl_solver.avl.VRTX_R.RES[idx_res])
        enc_1 = copy.deepcopy(self.avl_solver.avl.VRTX_R.ENC[:, idx_res])

        dres_aincs_fd = (res_1 - res_0) / h
        denc_aincs_fd = (enc_1 - enc_0) / h
        
        # print('res ad', dres_aincs)
        # print('res fd', dres_aincs_fd)
        
        # print('enc ad', denc_aincs)
        # print('enc fd', denc_aincs_fd)
        
        
        np.testing.assert_allclose(
            np.sum(dres_aincs),
            np.sum(dres_aincs_fd),
            rtol=1e-8,
        )


class TestTotals(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver = AVLSolver(geo_file="aircraft.avl")
        self.avl_solver.add_constraint("alpha", 5.0)
        # self.avl_solver.add_constraint("beta", 9.0)
        # self.avl_solver.add_constraint("roll rate", 1.2)
        # self.avl_solver.add_constraint("pitch rate", 0.1)
        # self.avl_solver.add_constraint("yaw rate", 0.8)

        # base line CL
        self.avl_solver.execute_run()

    def test_cltot_alpha_derivs(self):

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
        dcl_dalfa += -1 * self.avl_solver.avl.CASE_R_diff.ALFA_diff

        # use finite difference
        coef_data = self.avl_solver.get_case_total_data()
        h = 1e-6
        alpha = self.avl_solver.get_avl_fort_arr("CASE_R", "ALFA")
        # self.avl_solver.set_avl_fort_arr("CASE_R", "ALFA", alpha + h)
        self.avl_solver.add_constraint("alpha", np.rad2deg(alpha) + h)
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        dcl_dalfa_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h
        dcl_dalfa_fd *= 180 / np.pi

        np.testing.assert_allclose(
            dcl_dalfa,
            dcl_dalfa_fd,
            rtol=1e-8,
        )

    def test_cltot_geom_derivs(self):
        # base line CL
        idx_surf = 0
        num_sec = self.avl_solver.avl.SURF_GEOM_I.NSEC[idx_surf]
        self.avl_solver.avl.CASE_R_diff.CLTOT_diff = 1
        self.avl_solver.avl.aero_b()
        self.avl_solver.avl.update_surfaces_b()
        # print('python nvor', self.avl_solver.avl.CASE_I.NVOR)
        # quit()
        # dcl_geom = self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff
        # print("dcl_geom", dcl_geom[:1], np.linalg.norm(dcl_geom))


        dcl_daincs = copy.deepcopy(self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff[:num_sec, idx_surf])
        
        rhs =   copy.deepcopy(self.avl_solver.avl.VRTX_R_diff.GAM_diff[:num_sec])
        # solve for the adjoint variable
        self.avl_solver.avl.solve_adjoint()
        dcl_dres = self.avl_solver.avl.VRTX_R_diff.RES_diff

        # combine adjoint with pr/px
        
        self.avl_solver.avl.get_res_b()
        self.avl_solver.avl.update_surfaces_b()
        dcl_daincs += -1 * self.avl_solver.avl.SURF_GEOM_R_diff.AINCS_diff[:num_sec, idx_surf]

        # # use finite difference
        coef_data = self.avl_solver.get_case_total_data()

        h = 1e-8
        self.avl_solver.avl.SURF_GEOM_R.AINCS[:num_sec, idx_surf] += h
        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.aero()
        coef_data_peturb = self.avl_solver.get_case_total_data()

        dcl_daincs_fd = (coef_data_peturb["CL"] - coef_data["CL"]) / h
        
        
        np.testing.assert_allclose(
            np.sum(dcl_daincs),
            dcl_daincs_fd,
            rtol=1e-8,
        )


if __name__ == "__main__":
    unittest.main()

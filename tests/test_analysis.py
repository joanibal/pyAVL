# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver

# =============================================================================
# Standard Python Modules
# =============================================================================
import os
import copy
# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")


class TestAnalysisSweep(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_constrained_alpha_sweep(self):
        self.avl_solver.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
        self.avl_solver.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")

        alpha_array = np.arange(0, 10)
        cl_ref_arr = np.array(
            [
                1.2349061238204873,
                1.3282779560514677,
                1.42113603863798,
                1.5134418218893897,
                1.6051572516385355,
                1.6962448080698307,
                1.7866675434128751,
                1.8763891184181036,
                1.9653738375342211,
                2.053586682710641,
            ]
        )
        cd_ref_arr = np.array(
            [
                0.030753634421301752,
                0.033676747268002544,
                0.03679211508827116,
                0.040092804096450746,
                0.04357129507051221,
                0.047219503496354835,
                0.05102880128527856,
                0.0549900400052339,
                0.05909357556230493,
                0.06332929426488881,
            ]
        )
        for idx_alpha, alpha in enumerate(alpha_array):
            self.avl_solver.add_constraint("alpha", alpha)

            self.avl_solver.execute_run()
            run_data = self.avl_solver.get_case_total_data()

            np.testing.assert_allclose(
                cl_ref_arr[idx_alpha],
                run_data["CL"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                cd_ref_arr[idx_alpha],
                run_data["CD"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                run_data["CM"],
                0.0,
                atol=1e-8,
            )

    def test_constrained_cl_sweep(self):
        self.avl_solver.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
        self.avl_solver.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")

        cd_ref_arr = np.array(
            [
                0.01654814255244833,
                0.018124778345018383,
                0.01994896285331091,
                0.02202067604738557,
                0.024339528639367926,
                0.026904749658315106,
                0.029715171166770832,
                0.03276920986323026,
                0.03606484526169707,
                0.03959959407831794,
                0.043370480383046583,
            ]
        )
        cl_arr = np.arange(0.6, 1.7, 0.1)
        for idx_cl, cl in enumerate(cl_arr):
            self.avl_solver.add_trim_condition("CL", cl)
            self.avl_solver.execute_run()
            run_data = self.avl_solver.get_case_total_data()

            np.testing.assert_allclose(
                cl,
                run_data["CL"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                cd_ref_arr[idx_cl],
                run_data["CD"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                run_data["CM"],
                0.0,
                atol=1e-8,
            )


class TestHingeMom(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_con_surf_mom(self):
        self.avl_solver.add_constraint("Elevator", 10.00)
        self.avl_solver.add_constraint("Rudder", 0.00)
        self.avl_solver.execute_run()

        mom_data = self.avl_solver.get_hinge_moments()

        np.testing.assert_allclose(mom_data["Elevator"], -0.04381216304, rtol=1e-8)
        np.testing.assert_allclose(mom_data["Rudder"], 0.0, atol=1e-8)

        self.avl_solver.add_constraint("Elevator", 0.00)
        self.avl_solver.add_constraint("Rudder", 10.00)
        self.avl_solver.execute_run()
        mom_data = self.avl_solver.get_hinge_moments()

        np.testing.assert_allclose(mom_data["Rudder"], -1.0957068091302286e-3, rtol=1e-8)
        np.testing.assert_allclose(mom_data["Elevator"], -1.065702741142345e-2, rtol=1e-8)

class TestCaseDerivs(unittest.TestCase):
    def setUp(self) -> None:
        # self.avl_solver = AVLSolver(geo_file=geom_file)
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver = AVLSolver(geo_file="rect.avl", debug=False)
        
    def test_coefs_wrt_con_surfs(self):
        var = "CD"
        cs = "Elevator"
        
        self.avl_solver.add_constraint("alpha", 10.00)
        self.avl_solver.add_constraint("Elevator", 0.00)
        # self.avl_solver.execute_run()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()

        
        
        run_data = self.avl_solver.get_case_total_data()
        coef_derivs = self.avl_solver.get_case_coef_derivs()
        var_base = run_data[var]
        
        dvar_dcs = coef_derivs[var][cs]
        
        h = 1e-0
        self.avl_solver.add_constraint("Elevator", h)
        # self.avl_solver.execute_run( )
        self.avl_solver.avl.exec_rhs( )
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()

        run_data = self.avl_solver.get_case_total_data()
        var_perturb = run_data[var]
        
        dvar_dcs_fd = (var_perturb - var_base)/h
        
        print(f"h:{h} AD:{dvar_dcs} FD:{dvar_dcs_fd}")
        
class TestGamD(unittest.TestCase):
    def setUp(self) -> None:
        # self.avl_solver = AVLSolver(geo_file=geom_file)
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver = AVLSolver(geo_file="rect.avl", debug=False)

    def test_gam_d(self):
        var = "CD"
        cs = "Elevator"
        
        self.avl_solver.add_constraint("alpha", 10.00)
        self.avl_solver.add_constraint("Elevator", 0.00)
        # self.avl_solver.execute_run()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        
        run_data = self.avl_solver.get_case_total_data()
        coef_derivs = self.avl_solver.get_case_coef_derivs()
        var_base = run_data[var]
        
        dvar_dcs = coef_derivs[var][cs]

        
        num_gam = self.avl_solver.get_mesh_size()
        slicer_gam = (slice(0, self.avl_solver.get_mesh_size()))
        num_cs = self.avl_solver.get_num_control_surfs()
        slicer_gam_d = (slice(0, self.avl_solver.get_num_control_surfs()), slice(0, self.avl_solver.get_mesh_size()))
        
        gam_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("VRTX_R", "GAM", slicer=slicer_gam))
        dcp_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("VRTX_R", "DCP", slicer=slicer_gam))
        gam_d = self.avl_solver.get_avl_fort_arr("VRTX_R", "GAM_D", slicer=slicer_gam_d)
        dcp_d = self.avl_solver.get_avl_fort_arr("VRTX_R", "DCP_D", slicer=slicer_gam_d)
        cdstrp_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CDSTRP"))
        cxstrp_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CXSTRP"))
        czstrp_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CZSTRP"))
        cdstrp_d = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CDST_D"))
        cxstrp_d = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CXST_D"))
        czstrp_d = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CZST_D"))
        
        wv_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("SOLV_R", "wv", slicer=slicer_gam))
        wc_base = copy.deepcopy(self.avl_solver.get_avl_fort_arr("SOLV_R", "wc", slicer=slicer_gam))
        
        wv_d = self.avl_solver.get_avl_fort_arr("SOLV_R", "WV_D", slicer=slicer_gam_d)
        wc_d = self.avl_solver.get_avl_fort_arr("SOLV_R", "WC_D", slicer=slicer_gam_d)
        
        
        h = 1e-3
        self.avl_solver.add_constraint("Elevator", h)
        # self.avl_solver.execute_run( )
        self.avl_solver.avl.exec_rhs( )
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        
        gam_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("VRTX_R", "GAM", slicer=slicer_gam))
        dcp_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("VRTX_R", "DCP", slicer=slicer_gam))
        wv_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("SOLV_R", "wv", slicer=slicer_gam))
        wc_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("SOLV_R", "wc", slicer=slicer_gam))
        cdstrp_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CDSTRP"))
        cxstrp_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CXSTRP"))
        czstrp_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CZSTRP"))
        # cdstrp_perturb = copy.deepcopy(self.avl_solver.get_avl_fort_arr("STRP_R", "CDSTRP"))
        
        
        gam_d_fd = (gam_perturb - gam_base)/h
        wv_d_fd = (wv_perturb - wv_base)/h
        wc_d_fd = (wc_perturb - wc_base)/h
        dcp_d_fd = (dcp_perturb - dcp_base)/h
        run_data = self.avl_solver.get_case_total_data()
        var_perturb = run_data[var]
        
        dvar_dcs_fd = (var_perturb - var_base)/h
        
        dcdstrp_d_fd = (cdstrp_perturb - cdstrp_base)/h
        dcxstrp_d_fd = (cxstrp_perturb - cxstrp_base)/h
        dczstrp_d_fd = (czstrp_perturb - czstrp_base)/h
        
        
        print(f"gam AD:{gam_d} FD:{gam_d_fd}")
        
        print(f"wv AD:{wv_d} FD:{wv_d_fd}")
        print(f"wc AD:{wc_d} FD:{wc_d_fd}")
        
        print(f"dcp AD:{dcp_d} FD:{dcp_d_fd}")
        print(f"cdstrp AD:{cdstrp_d[0,0]} FD:{dcdstrp_d_fd[0]}" )
        print(f"cxstrp AD:{cxstrp_d[0,0]} FD:{dcxstrp_d_fd[0]}" )
        print(f"czstrp AD:{czstrp_d[0,0]} FD:{dczstrp_d_fd[0]}" )
        
        print(f"{var} AD:{dvar_dcs} FD:{dvar_dcs_fd}")



class TestVariableSetting(unittest.TestCase):
    def setUp(self) -> None:
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_alpha_set(self):
        """
        Test that setting the alpha works
        """
        
        self.avl_solver.add_constraint("alpha", 10.00)
        self.avl_solver.execute_run()
        alfa_new = self.avl_solver.get_case_parameter('alpha')
        np.testing.assert_allclose(alfa_new, 10.00, rtol=1e-15)
        

    def test_con_surf_set(self):
        """
        Test that setting the control surface works
        """
        
        self.avl_solver.add_constraint("Elevator", 10.00)
        self.avl_solver.add_constraint("alpha", 10.00)
        self.avl_solver.execute_run()
        alfa_new = self.avl_solver.get_case_parameter('alpha')
        np.testing.assert_allclose(alfa_new, 10.00, rtol=1e-15)
        def_dict = self.avl_solver.get_control_deflections()
        np.testing.assert_allclose(def_dict['Elevator'], 10.00, rtol=1e-15)
        
        

if __name__ == "__main__":
    unittest.main()

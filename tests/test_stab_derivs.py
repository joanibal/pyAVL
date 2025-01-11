# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver
import copy
import psutil

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


class TestStabDerivs(unittest.TestCase):
    # TODO: beta derivatives likely wrong

    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 5.0)
        # self.avl_solver.add_constraint("beta", 0.0)
        self.avl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def finite_dif(self, con_list, geom_seeds, step=1e-7):
        con_seeds = {}

        for con in con_list:
            con_seeds[con] = 1.0

        self.avl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
        self.avl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)

        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        # self.avl_solver.execute_run()
        coef_data_peturb = self.avl_solver.get_case_total_data()
        consurf_derivs_peturb = self.avl_solver.get_case_coef_derivs()

        self.avl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-step)
        self.avl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-step)

        self.avl_solver.execute_run()

        coef_data = self.avl_solver.get_case_total_data()
        consurf_derivs = self.avl_solver.get_case_coef_derivs()

        func_seeds = {}
        for func_key in coef_data:
            func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

        consurf_derivs_seeds = {}
        for func_key in consurf_derivs:
            consurf_derivs_seeds[func_key] = {}
            for surf_key in consurf_derivs[func_key]:
                consurf_derivs_seeds[func_key][surf_key] = (
                    consurf_derivs_peturb[func_key][surf_key] - consurf_derivs[func_key][surf_key]
                ) / step

        return func_seeds, consurf_derivs_seeds

    def test_deriv_values(self):
        # compare the analytical gradients with finite difference for each constraint and function
        base_data = self.avl_solver.get_case_total_data()
        stab_derivs = self.avl_solver.get_case_stab_derivs()

        con_keys =  ["alpha", "beta", "roll rate", "pitch rate", "yaw rate"]
        # con_keys =  ["alpha", "beta", 'roll rate']
        for con_key in con_keys:
            h = 1e-8
            val = self.avl_solver.get_case_constraint(con_key)
            self.avl_solver.add_constraint(con_key, val + h)
            self.avl_solver.execute_run()
            perb_data = self.avl_solver.get_case_total_data()
            self.avl_solver.add_constraint(con_key, val)
            
            for func_key in stab_derivs:
                ad_dot = stab_derivs[func_key][con_key] 
                fd_dot = (perb_data[func_key] - base_data[func_key]) / h 
                
                if con_key in ['alpha', 'beta']:
                   fd_dot *= 180/np.pi 
                   # convert to radians from degrees!
                
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))
                # print(f"{func_key:5} wrt {con_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-13
                if np.abs(ad_dot) < tol or np.abs(fd_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        atol=1e-10,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        rtol=5e-5,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )   
            
            


if __name__ == "__main__":
    unittest.main()
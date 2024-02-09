# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver, AVLGroup

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np
import openmdao.api as om
import warnings
import numpy as np

# Set DeprecationWarning to be treated as an error
warnings.simplefilter('error', DeprecationWarning)

# Optionally, you can also set NumPy options to raise errors for specific conditions
np.seterr(all='raise')

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")


class TestOMAnalaysis(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
    
        model = om.Group()
        model.add_subsystem("avlsolver", AVLGroup(geom_file=geom_file, mass_file=mass_file))

        self.prob = om.Problem(model)

    def test_aero_coef(self):
        
        self.avl_solver.execute_run()
        run_data = self.avl_solver.get_case_total_data()

        prob = self.prob
        prob.setup(mode='rev')
        prob.run_model()
        
        for func in run_data:        
            om_val = prob.get_val(f"avlsolver.{func}")
            assert om_val == run_data[func]
    
    def test_surface_param_setting(self):
        prob = self.prob
        prob.setup(mode='rev')
                
       
        surf_data = self.avl_solver.get_surface_params(include_geom=True, include_panneling=True, include_con_surf=True)
        np.random.seed(111)

        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                arr += np.random.rand(*arr.shape)*0.1
                # #set surface data
                self.avl_solver.set_surface_param(surf_key, geom_key, arr)
                self.avl_solver.avl.update_surfaces()
                self.avl_solver.execute_run()
                run_data = self.avl_solver.get_case_total_data()
                # set om surface data
                om_arr = prob.set_val(f"avlsolver.{surf_key}:{geom_key}", arr)
                prob.run_model()
                
                for func in run_data:        
                    om_val = prob.get_val(f"avlsolver.{func}")
                    assert om_val == run_data[func]
        
    def test_CL_solve(self):
        prob = self.prob
        cl_star = 1.5
        prob.model.add_design_var("avlsolver.alpha", lower=-10, upper=10)
        prob.model.add_constraint("avlsolver.CL", equals=cl_star)
        prob.model.add_objective("avlsolver.CD", scaler=1e3)
        prob.setup(mode='rev')
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['optimizer'] = 'SLSQP'
        prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
        prob.driver.options['tol'] = 1e-6
        prob.driver.options['disp'] = True
    
        prob.setup(mode='rev')
        om.n2(prob, show_browser=False, outfile="vlm_opt.html")
        prob.run_driver()
        om_val = prob.get_val(f"avlsolver.alpha")
        
        
        self.avl_solver.add_trim_condition("CL", cl_star)
        self.avl_solver.execute_run()
        alpha = self.avl_solver.get_case_parameter("alpha")
        
        np.testing.assert_allclose(om_val,
                            alpha,
                            rtol=1e-5,
                            err_msg=f"solved alpha",
                        )
        
    def test_CM_solve(self):
        prob = self.prob
        prob.model.add_design_var("avlsolver.alpha", lower=-10, upper=10)
        prob.model.add_constraint("avlsolver.CM", equals=0.0, scaler=1e3)
        prob.model.add_objective("avlsolver.CD", scaler=1e3)
        prob.setup(mode='rev')
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['optimizer'] = 'SLSQP'
        prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
        prob.driver.options['tol'] = 1e-6
        prob.driver.options['disp'] = True
    
        prob.setup(mode='rev')
        om.n2(prob, show_browser=False, outfile="vlm_opt.html")
        prob.run_driver()
        om_val = prob.get_val(f"avlsolver.alpha")
        
        
        self.avl_solver.add_constraint("alpha", 0.00, con_var="Cm pitch moment")
        self.avl_solver.execute_run()
        alpha = self.avl_solver.get_case_parameter("alpha")
        
        np.testing.assert_allclose(om_val,
                            alpha,
                            rtol=1e-5,
                            err_msg=f"solved alpha",
                        )

    def test_OM_total_derivs(self):
        prob = self.prob
        cl_star = 1.5
        prob.model.add_design_var("avlsolver.Wing:xyzles")
        prob.model.add_design_var("avlsolver.Wing:chords")
        prob.model.add_design_var("avlsolver.Wing:aincs")
        prob.model.add_design_var("avlsolver.Elevator", lower=-10, upper=10)
        prob.model.add_design_var("avlsolver.alpha", lower=-10, upper=10)
        prob.model.add_constraint("avlsolver.CL", equals=cl_star)
        prob.model.add_objective("avlsolver.CD", scaler=1e3)
        prob.setup(mode='rev')
        om.n2(prob, show_browser=False, outfile="vlm_opt.html")
        prob.run_model()
        deriv_err = prob.check_totals()
        rtol = 1e-4
        for key, data in deriv_err.items():
                np.testing.assert_allclose(
                    data['J_fd'],
                    data['J_rev'],
                    rtol=rtol,
                    err_msg=f"deriv of {key[0]} wrt {key[1]} does not agree with FD to rtol={rtol}"
                )
        
if __name__ == "__main__":
    unittest.main()
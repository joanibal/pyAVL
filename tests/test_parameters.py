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


class TestParameterAPI(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.params_baseline = {
            # "bank": 0.0,
            # "elevation": 0,
            "Mach": 0.1,
            "velocity": 0.0,
            "density": 1.225,
            "grav.acc.": 9.81,
            "mass": 5.63,
            "Ixx": 0.196,
            "Iyy": 1.54,
            "Izz": 1.74,
            "X cg": 0.0654,
            "Y cg": 0,
            "Z cg": 0,
            "CD0": 0.0116,
        }

    def test_get_parameters(self):

        for key in self.params_baseline:
            param = self.avl_solver.get_case_parameter(key)
            self.assertEqual(param, self.params_baseline[key], msg=key)

    def test_set_parameters(self):

        for key in self.params_baseline:
            # add each key to the update dict one at a time

            self.avl_solver.set_case_parameter(key, self.params_baseline[key] + 0.1)

            # check that the parameter was updated
            param_updated = self.avl_solver.get_case_parameter(key)
            self.assertEqual(param_updated, self.params_baseline[key] + 0.1, msg=key)

        # make sure the parameters stay set after an update
        self.avl_solver.execute_run()
        for key in self.params_baseline:
            param_updated = self.avl_solver.get_case_parameter(key)
            self.assertEqual(param_updated, self.params_baseline[key] + 0.1, msg=key)

    def test_set_fort_var(self):
        """
            test that the parameter changes effect the correct fortran variables
        """
        
        for key in self.params_baseline:
            # add each key to the update dict one at a time

            self.avl_solver.set_case_parameter(key, self.params_baseline[key] + 0.1)

        # other parameters only get updated in exec subroutine
        self.avl_solver.execute_run()
        
        self.assertEqual(self.avl_solver.get_avl_fort_arr("CASE_R", "MACH"), 
                         self.params_baseline["Mach"] + 0.1)
        self.assertEqual(self.avl_solver.get_avl_fort_arr("CASE_R", "CDREF"), 
                         self.params_baseline["CD0"] + 0.1)
        
        xyz_ref = self.avl_solver.get_avl_fort_arr("CASE_R", "XYZREF")
        self.assertEqual(xyz_ref[0], self.params_baseline["X cg"] + 0.1)
        self.assertEqual(xyz_ref[1], self.params_baseline["Y cg"] + 0.1)
        self.assertEqual(xyz_ref[1], self.params_baseline["Z cg"] + 0.1)
        
              

if __name__ == "__main__":
    unittest.main()
    
 
# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver
import copy
import resource

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


avl_solver = AVLSolver(geo_file=geom_file)
class TestCaseDerivs(unittest.TestCase):
    def setUp(self) -> None:
        # self.avl_solver = AVLSolver(geo_file=geom_file)
        self.avl_solver = avl_solver
        pass
    
    def tearDown(self):
        # Without the following line a copy of large_list will be kept in
        # memory for each test that runs, uncomment the line to allow the
        # large_list to be garbage collected.
        # self.large_list = None
        mb_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
        print("Memory usage: %s MB" % mb_memory)

    def test_1(self):
        self.avl_solver.execute_run()

    def test_2(self):
        self.avl_solver.execute_run()



if __name__ == "__main__":
    unittest.main()
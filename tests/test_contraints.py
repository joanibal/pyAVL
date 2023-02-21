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


class TestConstraints(unittest.TestCase):
    # TODO: add reference values for comparison by running avl binnary by hand or with wrappper

    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    # TODO: add test for roll, pitch, and yaw rate
    # def test_rates(self):

    def test_angles(self):
        self.avl_solver.addConstraint("alpha", 6.00)
        self.avl_solver.addConstraint("beta", 2.00)
        self.avl_solver.executeRun()

    def test_control_surfaces(self):
        self.avl_solver.addConstraint("D1", 0.00)
        self.avl_solver.addConstraint("D2", 0.00)
        self.avl_solver.executeRun()

    def test_trim(self):
        # TODO: parametrize for the other options
        self.avl_solver.addTrimCondition("CL", 1.0)
        self.avl_solver.executeRun()


if __name__ == "__main__":
    unittest.main()

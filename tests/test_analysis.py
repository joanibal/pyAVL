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


class TestAnalysisSweep(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_constrained_alpha_sweep(self):
        self.avl_solver.addConstraint('Elevator',  0.00, con_var='Cm pitch moment')
        self.avl_solver.addConstraint('Rudder', 0.00, con_var="Cn yaw moment")
        self.avl_solver.alphaSweep(0, 10)

        np.testing.assert_allclose(
            self.avl_solver.CL,
            np.array(
                [
                    1.229212,
                    1.322154,
                    1.414586,
                    1.506469,
                    1.597764,
                    1.688435,
                    1.778445,
                    1.867757,
                    1.956337,
                    2.044148,
                    2.131159,
                ]
            ),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CD,
            np.array(
                [
                    0.03059,
                    0.033488,
                    0.036577,
                    0.039849,
                    0.043298,
                    0.046915,
                    0.050691,
                    0.054619,
                    0.058687,
                    0.062887,
                    0.067207,
                ]
            ),
            rtol=1e-4,
        )
        np.testing.assert_allclose(self.avl_solver.CM, np.zeros_like(self.avl_solver.CM), atol=1e-8)

    def test_constrained_cl_sweep(self):
        self.avl_solver.addConstraint('Elevator',  0.00, con_var='Cm pitch moment')
        self.avl_solver.addConstraint('Rudder', 0.00, con_var="Cn yaw moment")
        start_cl = 0.6
        end_cl = 1.6
        increment = 0.1
        self.avl_solver.CLSweep(start_cl, end_cl, increment)

        np.testing.assert_allclose(self.avl_solver.CL, np.arange(start_cl, end_cl + increment, increment))
        np.testing.assert_allclose(
            self.avl_solver.CD,
            np.array(
                [
                    0.016546,
                    0.018124,
                    0.019949,
                    0.022023,
                    0.024343,
                    0.02691,
                    0.029722,
                    0.032778,
                    0.036076,
                    0.039612,
                    0.043385,
                ]
            ),
            rtol=1e-4,
        )
        np.testing.assert_allclose(self.avl_solver.CM, np.zeros_like(self.avl_solver.CM), atol=1e-8)

if __name__ == "__main__":
    unittest.main()
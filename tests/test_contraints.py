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

        np.testing.assert_allclose(
            np.rad2deg(self.avl_solver.alpha),
            np.array([6.0]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            np.rad2deg(self.avl_solver.beta),
            np.array([2.0]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CL,
            np.array([1.821967]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CD,
            np.array([0.053871]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CM,
            np.array([-0.195797]),
            rtol=1e-4,
        )

    def test_control_surfaces(self):
        self.avl_solver.addConstraint("D1", 10.00)
        self.avl_solver.addConstraint("D2", 5.00)
        self.avl_solver.executeRun()
        np.testing.assert_allclose(
            self.avl_solver.CL,
            np.array([1.005553]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CD,
            np.array([0.038557]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CM,
            np.array([0.956897]),
            rtol=1e-4,
        )
        
    def test_control_surfaces_names(self):
        """test that the control surface names are can be used as well"""
        self.avl_solver.addConstraint("Elevator", 10.00)
        self.avl_solver.addConstraint("Rudder", 5.00)
        self.avl_solver.executeRun()
        np.testing.assert_allclose(
            self.avl_solver.CL,
            np.array([1.005553]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CD,
            np.array([0.038557]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CM,
            np.array([0.956897]),
            rtol=1e-4,
        )

    def test_trim(self):
        # TODO: parametrize for the other options
        self.avl_solver.addTrimCondition("CL", 1.0)
        self.avl_solver.executeRun()


        np.testing.assert_allclose(
            self.avl_solver.alpha,
            np.array([-0.031734]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.beta,
            np.array([0.0]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CL,
            np.array([1.0]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CD,
            np.array([0.025359]),
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            self.avl_solver.CM,
            np.array([0.252352]),
            rtol=1e-4,
        )


if __name__ == "__main__":
    unittest.main()

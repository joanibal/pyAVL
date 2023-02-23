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
geom_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestGeom(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file)

    def test_surface_params(self):
        reference_data = {
            "Wing": {
                "nchordwise": 7,
                "cspace": 1.0,
                "nspan": 20,
                "sspace": -2.0,
                "yduplicate": 0.0,
                "scale": np.array([1.1, 1.2, 1.3]),
                "translate": np.array([0.1, 0.2, 0.3]),
                "angle": 1.23,
                "nspans": np.array([5, 4, 3, 2, 1]),
                "sspaces": np.array([-1.0, 0.0, 1.0, 2.0, 3.0]),
                "aincs": np.deg2rad(np.array([0.5, 0.4, 0.3, 0.2, 0.1])),
                "chords": np.array([0.5,0.4,0.3,0.2,0.1]),
                "xyzles": np.array([[0,   0,     0],
                                    [0.1, 1.0, 0.01],
                                    [0.2, 2.0, 0.02],
                                    [0.3, 3.0, 0.03],
                                    [0.4, 4.0, 0.04]]),
            },
        }

        data = self.avl_solver.get_surface_params()

        for surf in reference_data:
            for key in reference_data[surf]:
                # print(key, data[surf][key], reference_data[surf][key])
                np.testing.assert_allclose(
                    data[surf][key],
                    reference_data[surf][key],
                    atol=1e-8,
                    err_msg=f"Surface `{surf}` key `{key}` does not match reference data",
                )
        
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
            np.array([5.407351]),
            rtol=1e-4,
        )
        
        self.avl_solver.resetData()
        self.avl_solver.set_surface_params(data)                
        
        assert self.avl_solver.get_num_surfaces() == 5
        assert self.avl_solver.get_num_strips() == 90
        assert self.avl_solver.get_mesh_size() == 780
        
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
            np.array([5.407351]),
            rtol=1e-4,
        )

        # data
        # np.testing.assert_allclose(self.avl_solver.CM, np.zeros_like(self.avl_solver.CM), atol=1e-8)


if __name__ == "__main__":
    unittest.main()
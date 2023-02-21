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


class TestGetGEOM(unittest.TestCase):
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
                "scale": np.array([1.1,  1.2,  1.3]),
                "translate": np.array([0.1,  0.2,  0.3]),
                "angle": 1.23,
            },
        }
        
        data = self.avl_solver.get_surface_params()
        
        for surf in reference_data:
            for key in data[surf]:
                # print(key, data[surf][key], reference_data[surf][key])
                np.testing.assert_allclose(data[surf][key], reference_data[surf][key], atol=1e-8, 
                                           err_msg=f"Surface `{surf}` key `{key}` does not match reference data")
        
        # data
        # np.testing.assert_allclose(self.avl_solver.CM, np.zeros_like(self.avl_solver.CM), atol=1e-8)



if __name__ == "__main__":
    unittest.main()

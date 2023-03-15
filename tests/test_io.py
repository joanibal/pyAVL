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
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestInput(unittest.TestCase):
    # TODO: add test for expected input output errors
    def test_read_geom(self):
        avl_solver = AVLSolver(geo_file=geom_file)
        assert avl_solver.get_num_surfaces() == 5
        assert avl_solver.get_num_strips() == 90
        assert avl_solver.get_mesh_size() == 780

    def test_read_geom_and_mass(self):
        avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        assert avl_solver.get_avl_fort_var("CASE_L", "LMASS")
# TODO: add test for output (write out geom, mass, and stabilty)


class TestFortranLevelAPI(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_mod_file, mass_file=mass_file)

    def test_get_scalar(self):
        avl_version = 3.35

        version = self.avl_solver.get_avl_fort_var("CASE_R", "VERSION")
        self.assertEqual(version, avl_version)

        # test that this works with lower case
        version = self.avl_solver.get_avl_fort_var("case_r", "version")
        self.assertEqual(version, avl_version)

    def test_get_array(self):
        chords = self.avl_solver.get_avl_fort_var("SURF_GEOM_R", "CHORDS")

        self.assertEqual(chords.shape, (30, 400))
        np.testing.assert_array_equal(chords[0, :5], np.array([0.5, 0.4, 0.3, 0.2, 0.1]))


if __name__ == "__main__":
    unittest.main()

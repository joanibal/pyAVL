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
geom_output_file = os.path.join(base_dir, "aircraft_out.avl")

# TODO: add test for expected input output errors


class TestInput(unittest.TestCase):
    def test_read_geom(self):
        avl_solver = AVLSolver(geo_file=geom_file)
        assert avl_solver.get_num_surfaces() == 5
        assert avl_solver.get_num_strips() == 90
        assert avl_solver.get_mesh_size() == 780

    def test_read_geom_and_mass(self):
        avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        assert avl_solver.get_avl_fort_var("CASE_L", "LMASS")


class TestOutput(unittest.TestCase):
    def test_write_geom(self):
        """check that the file written by pyavl is the same as the original file"""
        avl_solver = AVLSolver(geo_file=geom_file)
        avl_solver.write_geom_file(geom_output_file)
        baseline_data = avl_solver.get_surface_params()

        del avl_solver
        avl_solver = AVLSolver(geo_file=geom_output_file)
        new_data = baseline_data = avl_solver.get_surface_params()

        for surf in baseline_data:
            for key in baseline_data[surf]:
                print(new_data[surf][key])
                data = new_data[surf][key]
                # check if it is a list of strings
                if isinstance(data, list) and isinstance(data[0], str):
                    for a, b in zip(data, baseline_data[surf][key]):
                        assert a == b
                else:
                    np.testing.assert_allclose(
                        new_data[surf][key],
                        baseline_data[surf][key],
                        atol=1e-8,
                        err_msg=f"Surface `{surf}` key `{key}` does not match reference data",
                    )


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

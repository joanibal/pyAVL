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

geom_file = os.path.join(base_dir, 'aircraft.txt')
mass_file = os.path.join(base_dir, 'aircraft.mass')

class TestInput(unittest.TestCase):
    # TODO: add test for expected input output errors
    def test_read_geom(self):
        AVLSolver(geo_file=geom_file)
    
    def test_read_geom_and_mass(self):
        AVLSolver(geo_file=geom_file, mass_file=mass_file)
        
    
# TODO: add test for output (write out geom, mass, and stabilty) 


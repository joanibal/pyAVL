from pyavl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

avl_solver = AVLSolver(geo_file="aircraft.avl")
avl_solver.plot_geom()

# # set the angle of attack
# avl_solver.add_constraint("alpha", 5.00)
# avl_solver.set_case_parameter("Mach", 0.3)

# # This is the method that acutally runs the analysis
# avl_solver.execute_run()

# avl_solver.plot_cp()

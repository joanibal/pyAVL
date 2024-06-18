from pyavl import AVLSolver
import numpy as np


avl_solver = AVLSolver(geo_file="aircraft.avl", mass_file="aircraft.mass",  debug=True)

# set the angle of attack
avl_solver.add_constraint("alpha", 5.00)

avl_solver.set_case_parameter("velocity", 10)

# This is the method that acutally runs the analysis
avl_solver.set_avl_fort_arr('CASE_R', 'EXEC_TOL', 0.00002)
avl_solver.execute_run()
avl_solver.avl.mode()


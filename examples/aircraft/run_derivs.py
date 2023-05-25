from pyavl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)
alpha = 3.0
avl_solver.add_constraint("alpha", alpha)
avl_solver.execute_run()
run_data = avl_solver.get_case_total_data()
avl_solver.get_deriv()

step = 1e-6
avl_solver.add_constraint("alpha", alpha + step)
avl_solver.execute_run()
run_data_p = avl_solver.get_case_total_data()
print("fd", (run_data_p["CL"] - run_data["CL"]) / step)

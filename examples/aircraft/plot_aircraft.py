""" This scipt demos the use of the ways to vizualie the geometry and the solution in OptVL"""
from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)
avl_solver.plot_geom()

avl_solver.add_constraint("alpha", 5.00)
avl_solver.execute_run()
avl_solver.avl.cpoml()

avl_solver.plot_cp()



# ax.axis('equal')
# plt.show()

# # set the angle of attack
# avl_solver.add_constraint("alpha", 0.00)

# # control surface names from geometry file
# avl_solver.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
# avl_solver.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")

# avl_solver.set_case_parameter("Mach", 0.3)

# # This is the method that acutally runs the analysis
# avl_solver.execute_run()
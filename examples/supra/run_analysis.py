from optvl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="supra.avl", debug=False)
# set the angle of attack
avl_solver.set_case_parameter("Mach", 0.0)

print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    avl_solver.add_constraint("alpha", alpha)
    avl_solver.execute_run()
    run_data = avl_solver.get_case_total_data()
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )




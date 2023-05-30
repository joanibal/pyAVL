from pyavl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)
avl_solver.add_constraint("alpha", 0.00)

# control surface names from geometry file
avl_solver.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
avl_solver.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")

avl_solver.set_case_parameter("Mach", 0.3)

# # This is the method that acutally runs the analysis
# avl_solver.execute_run()

# print("----------------- alpha sweep ----------------")
# print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
# for alpha in range(10):
#     avl_solver.add_constraint("alpha", alpha)
#     avl_solver.execute_run()
#     run_data = avl_solver.get_case_total_data()
#     print(
#         f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
#     )


# avl_solver.CLSweep(0.6, 1.6)

print("----------------- CL sweep ----------------")
print("   Angle        Cl           Cd          Cdff          Cdv          Cm    CN")
for cl in np.arange(0.6, 1.7, 0.001):
    avl_solver.add_trim_condition("CL", cl)
    avl_solver.execute_run()
    run_data = avl_solver.get_case_total_data()
    alpha = avl_solver.get_case_parameter("alpha")
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f} {run_data["CN SA"]:10.6f} '
    )
    
avl_solver.add_trim_condition("CL", 1.6)
avl_solver.execute_run()
run_data = avl_solver.get_case_total_data()
alpha = avl_solver.get_case_parameter("alpha")
print(
    f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f} {run_data["CN SA"]:10.6f} '
)

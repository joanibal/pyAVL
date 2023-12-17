from pyavl import AVLSolver
import numpy as np
import time
start_time = time.time()
last_time = time.time()
avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False, timing=True)
alpha = 3.0
avl_solver.add_constraint("alpha", alpha)
print('--- analysis ---')
avl_solver.execute_run()


run_data = avl_solver.get_case_total_data()
print(f"analysis time: {time.time()-last_time}")
last_time = time.time()
np.testing.assert_allclose(
                1.5031447347753768,
                run_data["CL"],
                rtol=3e-14,
            )
# print(run_data["CL"])
# print('--- sens ---')

# sens = avl_solver.execute_run_sensitivies(['CL'])

# step = 1e-6
# avl_solver.add_constraint("alpha", alpha + step)
# print(f"sens time: {time.time()-last_time}")
# last_time = time.time()
# print('--- analysis ---')
# avl_solver.execute_run()
# run_data_p = avl_solver.get_case_total_data()
# # print("fd     ", (run_data_p["CL"] - run_data["CL"]) / step)
# # print("adjoint", sens['CL']['alpha'])
# print(f"analysis time: {time.time()-last_time}")
# last_time = time.time()
# print(f"runscript time: {time.time()-start_time}")
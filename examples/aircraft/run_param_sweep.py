from pyavl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False, timing=False)

# set the angle of attack
avl_solver.add_constraint("alpha", 5.00)

for y_scale in np.linspace(0.5, 1.5, 5):
    avl_solver.set_surface_params({"Wing":{"scale":np.array([1, y_scale, 1])}})

    avl_solver.execute_run()
    stab_derivs = avl_solver.get_case_stab_derivs()

    print(f"----------------- y_scale: {y_scale} ----------------")
    for func in stab_derivs:
        for con_key in stab_derivs[func]:
            print(f"d{func}/d{con_key:5}: {stab_derivs[func][con_key]:.6f}")
    
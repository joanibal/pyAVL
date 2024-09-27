from pyavl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=True)
avl_solver.avl.cpoml()

xyz_lo, xyz_up, cp_lo, cp_up = avl_solver.get_cp_data()
# import pdb; pdb.set_trace()
# print(xyz_lo.shape)
# for isurf in range(xyz_lo.shape[2]):
#     print('isurf', isurf)
#     for i, xyz in enumerate(xyz_lo[isurf]):
#         print(i, 'XYZ LO', xyz[0], xyz[1], xyz[2])
num_surf = avl_solver.get_num_surfaces()
surf_names = avl_solver.get_surface_names()

for idx_surf in range(num_surf):
    print(surf_names[idx_surf], xyz_lo[idx_surf].shape, 'zeros', np.sum(xyz_lo[idx_surf].flatten() == 0.0))
    xyzs= xyz_lo[idx_surf]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for xyz in xyzs:
        ax.scatter(xyz[0], xyz[1] , xyz[2])

    xyzs= xyz_up[idx_surf]
    for xyz in xyzs:
        ax.scatter(xyz[0], xyz[1] , xyz[2])


    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    # ax.set_box_aspect([1,1,1])

    plt.show()

# fig, ax = plt.subplots(1, 1, sharex=True)
# avl_solver.add_mesh_plot(ax)


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
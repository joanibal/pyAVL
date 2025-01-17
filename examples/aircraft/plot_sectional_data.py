from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

ovl = AVLSolver(geo_file="aircraft.avl", debug=False)
ovl.add_constraint('alpha', 5.0)
ovl.execute_run()

strip_data = ovl.get_strip_data()
for surf_key in strip_data:
    span_distance = strip_data[surf_key]['XYZ LE'][:,1]
    plt.plot(span_distance, strip_data[surf_key]['lift dist'])

plt.ylabel('lift distribution')
plt.xlabel('spanwise position')
plt.show()


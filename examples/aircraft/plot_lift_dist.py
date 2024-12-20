from pyavl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pprint
avl = AVLSolver(geo_file="aircraft.avl")

# set the angle of attack
avl.add_constraint("alpha", 5.00)

avl.set_case_parameter("Mach", 0.3)

# This is the method that acutally runs the analysis
avl.execute_run()

data_dict = avl.get_strip_data()

for surf in data_dict:
    
    # plot just the unduplicated surfaces so we can see a little easier
    if "DUP" in surf:
        continue
    
    xyz_le = data_dict[surf]['XYZ LE']
    
    plt.plot(xyz_le[:,1], data_dict[surf]['lift dist'])
    plt.plot(xyz_le[:,1], data_dict[surf]['CL perp'])
    plt.plot(xyz_le[:,1], data_dict[surf]['CL'],'--')
    
plt.show()
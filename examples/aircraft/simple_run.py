from pyavl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)
avl_solver.addConstraint('alpha', 0.00)

# surface names form geometry file
avl_solver.addConstraint('Elevator',  0.00, con_var='Cm pitch moment')
avl_solver.addConstraint('Rudder', 0.00, con_var="Cn yaw moment")

avl_solver.executeRun()

print('----------------- Neutral Point ----------------')
avl_solver.calcNP()
print(avl_solver.NP, ' X Np')
avl_solver.resetData()


avl_solver.alphaSweep(0, 10)

print('----------------- alpha sweep ----------------')
print('   Angle        Cl           Cd          Cdff          Cdv          Cm')
for i in range(len(avl_solver.alpha)):
    print(' %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   '%(avl_solver.alpha[i]*(180/np.pi),avl_solver.CL[i],avl_solver.CD[i], avl_solver.CDFF[i], avl_solver.CDV[i],  avl_solver.CM[i]))


avl_solver.resetData()

avl_solver.CLSweep(0.6, 1.6)

print('----------------- CL sweep ----------------')
print('   Angle        Cl           Cd          Cdff          Cdv          Cm')
for i in range(len(avl_solver.alpha)):
    print(' %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   '%(avl_solver.alpha[i]*(180/np.pi),avl_solver.CL[i],avl_solver.CD[i], avl_solver.CDFF[i], avl_solver.CDV[i],  avl_solver.CM[i]))



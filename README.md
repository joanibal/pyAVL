# pyAVL

pyAVL is a stripped down version of Mark Drela and Harold Youngren's famous AVL code wrapped in python with f2py.
This allows one to more easily conduct large parameter sweeps in AVL or to include AVL into a larger model. 
Additionally, this wrapper provides access to more data than is available through traditional file output. 

# Installation
To compile the avl library use 
```
make
```
This code has only been tested with gfortran and gnu95 compilers. 
If you want to use something besides gfortran you will have to modify the Makefile


and to install the pyavl package on your python path use 
```
pip install . 
```
or 
```
pip install . -e 
```
to install in development mode 

# Basic usage
The API of pyAVL was made to mirror the usage of AVL through its text interface. 

The AVL wrapper is implemented in the `AVLSolver` class. 
To use this wrapper, first one must initialize the `AVLSoilver` object with a geometry file and optionally a mass file. 
After, the user can add constraints and then execute the run to generate data. 
Below is a basic example of this workflow. 
```
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

```



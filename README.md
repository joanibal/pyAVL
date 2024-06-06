# pyAVL
[![Downloads](https://static.pepy.tech/badge/pyavl-wrapper)](https://pepy.tech/project/pyavl-wrapper)
<!--- [![Downloads](https://static.pepy.tech/badge/pyavl-wrapper/month)](https://pepy.tech/project/pyavl-wrapper) --->
[Documentation](https://joanibal.github.io/pyAVL/)

pyAVL is a stripped down version of Mark Drela and Harold Youngren's famous AVL code wrapped in python with f2py.
This allows one to more easily conduct large parameter sweeps in AVL or to include AVL into a larger model. 
Additionally, this wrapper provides access to more data than is available through traditional file output. 
Unlike in the output files which is limit to about 4 digits, the user has access to the full double precision data. 

# Installation
The best way to get pyAVL is to install it through pip
```
pip install pyavl-wrapper
```
This version even comes packaged with OpenBLAS for faster analysis. 


Currently, only Linux and macOS are supported. 
The process of building on Windows still has issues. 
For now Windows users will have to use pyAVL through Windows subsystem for Linux (WSL).


## building locally
If you want to make pyAVL locally then you have to clone the repository and use the following process.

In the root directory run
```
pip install . 
```

<!-- ## building step by step

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
to install in development mode  -->

# Basic usage
The API of pyAVL was made to mirror the usage of AVL through its text interface. 
The user loads in a geometry file, adds constraints, and then executes analysis runs.

The AVL wrapper is implemented in the `AVLSolver` class. 
To use this wrapper, first one must initialize the `AVLSolver` object with a geometry file and optionally a mass file. 
After, the user can add constraints and then execute the run to generate data. 
Below is a basic example of this workflow. 
```python
from pyavl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="aircraft.avl")
avl_solver.add_constraint("alpha", 0.00)

# control surface names from geometry file
avl_solver.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
avl_solver.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")

avl_solver.set_case_parameter("Mach", 0.3)

# This is the method that acutally runs the analysis
avl_solver.execute_run()

print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    avl_solver.add_constraint("alpha", alpha)
    avl_solver.execute_run()
    run_data = avl_solver.get_case_total_data()
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )

print("----------------- CL sweep ----------------")
print("   Angle        Cl           Cd          Cdff          Cdv          Cm")
for cl in np.arange(0.6,1.6,0.1):
    avl_solver.add_trim_condition("CL", cl)
    avl_solver.execute_run()
    run_data = avl_solver.get_case_total_data()
    alpha = avl_solver.get_case_parameter("alpha")
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )
```



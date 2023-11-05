# Basic Usage Overview

`pyAVL` offers a simple and intuitive API, mirroring AVL's text interface.
The typical workflow involves loading a geometry file, adding constraints, and executing analysis runs.

## Initializing and Setting up AVL Solver
To begin with `pyAVL`, start by initializing the `AVLSolver` class:

```python
from pyavl import AVLSolver
```
avl_solver = AVLSolver(geo_file="aircraft.avl")

After initializing, you can set up various constraints and parameters:

```python
avl_solver.add_constraint("alpha", 0.00)
avl_solver.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
avl_solver.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")
avl_solver.set_case_parameter("Mach", 0.3)
```
## Running Analysis

Once you've set up the solver, running the analysis is straightforward:

```python
avl_solver.execute_run()
```

For a more detailed example and advanced use cases, see the analysis guide.
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



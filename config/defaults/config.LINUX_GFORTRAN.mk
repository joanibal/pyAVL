# Config File for LINUX and GFORTRAN Compiler
AR       = ar
AR_FLAGS = -rvs

FF90       = gfortran
FF90_FLAGS = -O2 -fdefault-real-8 -fPIC -Wno-align-commons -std=legacy

C_FLAGS = -O2 -fPIC -g

F2PY = f2py
F2PY_FF90 = gfortran

PYTHON = python

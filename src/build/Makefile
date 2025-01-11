# Include user supplied configuration file
include ../../config/config.mk

# Group Fortran compiler flags
FF90_ALL_FLAGS = $(FF90_FLAGS)

CC_ALL_FLAGS     = $(C_FLAGS)
# Include full list of files
include fileList

# Include full list of directories
include directoryList

# Include the list of rules
include rules

# Set the make VPATH variable to the "dirs" variable from
# directorylist. We must first append the '../'
dirs:=$(addprefix ../,$(dirs))
VPATH:=$(dirs)

# Add ../ prefix to all source files, remove directories
fortranFiles:=$(addprefix ../,$(fortranFiles))
cFiles:=$(addprefix ../,$(cFiles))
f77Files:=$(addprefix ../,$(f77Files))

fortranFilesNoDir=$(notdir $(fortranFiles))
cFilesNoDir=$(notdir $(cFiles))
f77FilesNoDir=$(notdir $(f77Files))

# Generate two separate list of .F90 and .f90 files using the filter command
f90Files=$(filter %.f90,$(fortranFilesNoDir))
F90Files=$(filter %.F90,$(fortranFilesNoDir))

# Finally convert all source files to .o
OFILES=$(f90Files:%.f90=%.o) $(F90Files:%.F90=%.o) $(f77FilesNoDir:%.f=%.o) $(cFilesNoDir:%.c=%.o)

OS := $(shell uname)

# Define platform string based on the OS
ifeq ($(OS), Linux)
    OSNVMAX := 5450
else ifeq ($(OS), Darwin)
    OSNVMAX := 5450
else ifeq ($(OS), Windows_NT)
    OSNVMAX := 3500
else
    $(error Unsupported operating system)
endif

# Compile sources

# .PHONY: ../include/AVL.INC sources lib

default: lib ../f2py/libavl.pyf
	$(F2PY) --fcompiler=$(F2PY_FF90) --f90flags='$(FF90_ALL_FLAGS)' -c -m  libavl ../f2py/libavl.pyf libavl.a
	python importTest.py libavl
	mv libavl.*so ../../pyavl/libavl.so


# # Rule to generate AVL.INC
../includes/AVL.INC: ../includes/AVL.INC.in
	python ../create_avl_inc.py $(OSNVMAX)


# Ensure all object files are recompiled if AVL.INC changes
# $(OFILES): ../includes/AVL.INC

sources: $(OFILES) 

lib: sources
	$(AR) $(AR_FLAGS) libavl.a $(OFILES) 

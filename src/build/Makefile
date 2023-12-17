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

.PHONY: sources lib

default: lib ../f2py/libavl.pyf
	$(F2PY) -llapack -lblas -lgomp --fcompiler=$(F2PY_FF90) --f90flags='$(FF90_ALL_FLAGS)' -c -m  libavl ../f2py/libavl.pyf libavl.a 
	python importTest.py libavl
	mv libavl.*so ../../pyavl/libavl.so

sources: $(OFILES)

lib: sources
	$(AR) $(AR_FLAGS) libavl.a $(OFILES) 

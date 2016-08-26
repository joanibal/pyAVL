#*********************************************************
# Makefile for XFOIL V6.93 programs
# H.Youngren 4/24/01
# M.Drela
#*********************************************************
INCLUDE_COMMON = ../common.mk
include ${INCLUDE_COMMON}

AVL_LIB = libavl.a
OBJ =	aero.o\
	 aic.o\
	 ainput.o\
	 airutil.o\
	 amake.o\
	 amass.o\
	 amode.o\
	 aoper.o\
	 aoutput.o\
	 asetup.o\
	 atpforc.o\
	 atrim.o\
	 autil.o\
	 avl.o\
	 cdcl.o\
	 getvm.o\
	 hidden.o\
	 matrix.o\
	 matrix-lapack.o\
	 matrix-linpack.o\
	 matrix-lapack.o\
	 second.o\
	 second_g77.o\
	 second_ifc.o\
	 sgutil.o\
	 spline.o\
	 userio.o\
	 eispack.o\




default: libavl

clean:
	$(RM) *.o libavl.a

libavl: $(OBJ) 
	$(AR) $(AR_FLAGS) $(AVL_LIB) $(OBJ)
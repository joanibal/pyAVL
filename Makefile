#*********************************************************
# Makefile for pyXLIGHT and pyXLIGHT Complex
# G. Kenway Feb 7, 2010
#*********************************************************

default: 
	@echo "Please select one of the following configurations"
	@echo "Also ensure that config.tar.gz is extracted in config/"
	@echo "make intel       -> for linux pcs with intel compiler"
	@echo "make gfortran    -> for linux pcs with gfortran compiler"

# intel:
# 	@echo "Linux - Intel"
# 	-rm common.mk
# 	ln -s ./config/config.LINUX_INTEL.mk common.mk
# 	( cd src && make) || exit 1; 
# 	( cd src_cs && make) || exit 1;
# 	f2py  --fcompiler=intel --f90flags=-r8 -c -m pyxlight src/pyxlight.pyf src/libxfoil.a
# 	f2py  --fcompiler=intel --f90flags=-r8 -c -m pyxlight_cs src_cs/pyxlight_cs.pyf src_cs/libxfoil_cs.a
# 	mv pyxlight.so ./python
# 	mv pyxlight_cs.so ./python
# 	-rm common.mk

default:
	@echo "Linux - Gfortran"
	-rm common.mk
	rm src/avl.o || :
	# rm src/aoper.o || :
	# rm src/amode.o || :
	# rm src/eispack.o || :

	rm python/pyavl.so || :
	ln -s ./config/config.LINUX_GFORTRAN.mk common.mk
	( cd src && make) || exit 1; 
	f2py  --fcompiler=gnu95 --f77flags='-O2 -fdefault-real-8 -fPIC -g -fcheck=all' -c -m  libavl src/libavl.pyf src/libavl.a 
	mv libavl.*so ./pyavl/libavl.so
	-rm common.mk


clean:
	rm -fr src/*.mod
	rm -fr src/*.o
	rm -fr src/*.a
	rm -fr src/*.so
	# rm -f *~ config.mk;


gfortran_test:
	@echo "Linux - Gfortran"
	-rm common.mk

	rm src_test/avl_test.o || :
	rm python/pyavl_test.so || :

	ln -s ./config/config.LINUX_GFORTRAN.mk common.mk
	( cd src_test && make) || exit 1; 
	f2py  --fcompiler=gfortran --f77flags=-fdefault-real8 -c -m  pyavl_test src_test/pyavl_test.pyf src_test/libavl_test.a 
	mv pyavl_test.so ./python
	-rm common.mk
# clean:
# 	-rm common.mk
# 	ln -s ./config/config.LINUX_INTEL.mk common.mk || exit 1;# Doesn't matter which 
# 	-(cd python && rm *.so)
# 	(cd src && $(MAKE) $@) || exit 1;
# 	(cd src_cs && $(MAKE) $@) || exit 1;
# 	-rm common.mk

# Note we have to copy a dummy config file to common.mk so it doesn't complain
# about the file not existing
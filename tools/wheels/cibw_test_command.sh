set -xe

PROJECT_DIR="$1"

cd $PROJECT_DIR/tests
# python -m unittest -v
#HACK: if the tests are not split up the CI runs out of memory...

# test package built and installed correctly
python -m unittest -v test_import.py
python -m unittest -v test_io.py

# test basic avl functionality
python -m unittest -v test_parameters.py
python -m unittest -v test_analysis.py
python -m unittest -v test_surf_geom.py
python -m unittest -v test_contraints.py

# tests for adjoint
python -m unittest -v test_new_subroutines.py
python -m unittest -v test_partial_derivs.py
python -m unittest -v test_total_derivs.py

# test mem ussage of pyavl and test framework
python -m unittest -v test_tear_down.py


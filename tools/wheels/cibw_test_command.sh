set -xe

PROJECT_DIR="$1"

cd $PROJECT_DIR/tests
python -m unittest -v

python -m unittest -v test_analysis.py  &&python -m unittest -v test_io.py     &&
python -m unittest -v test_surf_geom.py &&
python -m unittest -v test_contraints.py  &&python -m unittest -v test_import.py      &&
python -m unittest -v test_parameters.py
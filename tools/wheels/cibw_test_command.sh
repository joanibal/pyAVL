set -xe

PROJECT_DIR="$1"

cd $PROJECT_DIR/tests
python -m unittest -v


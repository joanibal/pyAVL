#!/bin/bash

# Find the .ini file in the specified directory pattern
file_path=$(find /Users/runner/work/pyAVL/pyAVL -name 'meson-python-native-file.ini' -print -quit)

# Output the found file path
echo $file_path

# Installation Guide for pyAVL

## Installing with pip
The recommended method to install `pyAVL` is via pip:

```shell
pip install pyavl-wrapper
```
This package is optimized with OpenBLAS for quicker analysis.

### Supported Platforms
Currently, `pyAVL` supports:
- Linux
- macOS

!!! Note
    Windows support is under development. For now, Windows users can utilize `pyAVL` through the Windows Subsystem for Linux (WSL).

## Building Locally
If you'd like to build `pyAVL` manually, follow the steps below:

1. Clone the repository to your local machine.
2. Install OpenBlas and make sure its libraries can be seen by the compiler. The following script should be able to install OpenBlas for you.
   ```shell
   sh -c 'bash /path/to/pyavl/tools/wheels/cibw_before_build_linux.sh /path/to/pyavl
   '```
3. Navigate to the root directory and run:
   ```
   pip install .
   ```

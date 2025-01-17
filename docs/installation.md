# Installation Guide for optvl

## Installing with pip
The recommended method to install `optvl` is via pip:

```shell
pip install optvl
```
This package is optimized with OpenBLAS for quicker analysis.

### Supported Platforms
Currently, `optvl` supports:
- Linux
- macOS (Apple Silicon and Intel)

!!! Note
    Windows support is under development. For now, Windows users can use `optvl` through the Windows Subsystem for Linux (WSL).

## Building Locally
If you'd like to build `optvl` manually, follow the steps below:

1. Clone the repository to your local machine.
2. Install OpenBlas and make sure its libraries can be seen by the compiler. The following script should be able to install OpenBlas for you.
   ```shell
   sh -c 'bash /path/to/optvl/tools/wheels/cibw_before_build_linux.sh /path/to/optvl
   '```
3. Navigate to the root directory and run:
   ```
   pip install .
   ```

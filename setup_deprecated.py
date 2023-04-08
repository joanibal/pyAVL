from numpy.distutils.core import setup
import re
import os
import subprocess
from numpy.distutils.core import Extension

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("pyavl/__init__.py").read(),
)[0]

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


def make_avl_lib():
    cwd = os.getcwd()
    # os.chdir("pyavl")
    # subprocess.run(["make", "clean"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    config_file = "config/defaults/config.LINUX_GFORTRAN.mk"

    build_dir = os.path.join(cwd, "src/build")
    subprocess.run(["cp", config_file, "config/config.mk"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    subprocess.run(["make", "clean"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p2 = subprocess.run(["make"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    compile_log = os.path.join(build_dir, "compile.log")
    with open(compile_log, "wb") as f:
        f.write(p2.stdout)
    if p2.returncode != 0:
        with open(compile_log, "r") as f:
            print(f.read())
        raise OSError("make", f"The compile command failed! Check the log at {compile_log} for more information.")


# open src/build/fileList makefile and return list of f77Files


ext_avl = Extension(
    "libavl",
    sources=[
        "src/AINDEX.INC"
        "src/AVL.INC"
        "src/aero.f",
        "src/aic.f",
        "src/ainput.f",
        "src/airutil.f",
        "src/amake.f",
        # "src/amass.f",
        # "src/amode.f",
        # "src/aoper.f",
        # "src/aoutput.f",
        # "src/asetup.f",
        # "src/atpforc.f",
        # "src/atrim.f",
        # "src/autil.f",
        # "src/avl.f",
        # "src/cdcl.f",
        # "src/getvm.f",
        # "src/hidden.f",
        # "src/matrix-lapackdp.f",
        # "src/second.f",
        # "src/second_g77.f",
        # "src/second_ifc.f",
        # "src/sgutil.f",
        # "src/spline.f",
        # "src/userio.f",
        # "src/eispack.f",
        "src/f2py/libavl.pyf",
    ],
    #    **lapack)
    f2py_options=[],
    extra_f90_compile_args=[
        "-O2",
        "-fdefault-real-8",
        "-fPIC",
        "-Wno-align-commons",
        "-std=legacy",
        "-lblas",
        "-llapack",
    ],
)


if __name__ == "__main__":
    # make_avl_lib()
    # cwd = os.getcwd()

    setup(
        name="pyavl",
        version=__version__,
        description="A direct Python interface for Mark Drela and Harold Youngren's famous AVL code.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/joanibal/pyavl",
        license="GPL-2.0",
        packages=[
            "pyavl",
        ],
        # package_data={"pyavl": ["*.so"]},
        ext_modules=[ext_avl],
        install_requires=["numpy"],
        extras_require={"plotting": ["matplotlib"], "testing": ["testflo>=1.4.5"]},
        classifiers=["Programming Language :: Python, Fortran"],
    )

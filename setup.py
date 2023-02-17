from setuptools import setup
import re
import os

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("pyavl/__init__.py").read(),
)[0]

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pyavl",
    version=__version__,
    description="A Python wrapped version of Mark Drela's AVL code with the GUI and eigen features.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/joanibal/pyavl",
    license="",
    packages=[
        "pyavl",
    ],
    package_data={"pyavl": ["*.so"]},
    install_requires=["numpy"],
    extras_require={
        "plotting": ["matplotlib"],
    },
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)

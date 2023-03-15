project(
  'pyavl',
  'c',
  version: '1.0.0',
  license: 'GPL-2.0',
  meson_version: '>= 0.64.0',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'fortran_std=legacy',
    'blas=openblas',
    'lapack=openblas'
  ],
)

cc = meson.get_compiler('c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
# if ff.has_argument('-Wno-conversion')
#   add_project_arguments('-Wno-conversion', language: 'fortran')
# endif


# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()


incdir_numpy = run_command(py3,
    [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
    ],
    check: true
).stdout().strip()
inc_np = include_directories(incdir_numpy)

# # TODO: pyoptsparse supports numpy>=1.16 but numpy.f2py.get_include() wasnt added until later, raise numpy version?
#incdir_f2py = run_command(py3_target,
#  [
#    '-c',
#    'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'
#  ],
#  check : true
#).stdout().strip()
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)

numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

fortranobject_c = incdir_f2py / 'fortranobject.c'
# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: numpy_nodepr_api,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)


# TODO: 64-bit BLAS and LAPACK
#
# Note that this works as long as BLAS and LAPACK are detected properly via
# pkg-config. By default we look for OpenBLAS, other libraries can be configured via
# `meson configure -Dblas=blas -Dlapack=lapack` (example to build with Netlib
# BLAS and LAPACK).
# For MKL and for auto-detecting one of multiple libs, we'll need a custom
# dependency in Meson (like is done for scalapack) - see
# https://github.com/mesonbuild/meson/issues/2835
blas_name = get_option('blas')
lapack_name = get_option('lapack')
# pkg-config uses a lower-case name while CMake uses a capitalized name, so try
# that too to make the fallback detection with CMake work
if blas_name == 'openblas'
  blas_name = ['openblas', 'OpenBLAS']
endif
if lapack_name == 'openblas'
  lapack_name = ['openblas', 'OpenBLAS']
endif
blas = dependency(blas_name)
lapack = dependency(lapack_name)

# TODO: Add `pybind11` when available as a dependency
dependency_map = {
  'BLAS': blas,
  'LAPACK': lapack,
}


# c = run_command('source' / 'grab-all-fortran-files.py', check: true)
# snopt_source_files = c.stdout().strip().split('\n')
avl_source_files = files(
        # 'src/AINDEX.INC',
        # 'src/AVL.INC',
        'src/aero.f',
        'src/aic.f',
        'src/ainput.f',
        'src/airutil.f',
        'src/amake.f',
        'src/amass.f',
        'src/amode.f',
        'src/aoper.f',
        'src/aoutput.f',
        'src/asetup.f',
        'src/atpforc.f',
        'src/atrim.f',
        'src/autil.f',
        'src/avl.f',
        'src/cdcl.f',
        'src/getvm.f',
        'src/hidden.f',
        'src/matrix-lapackdp.f',
        'src/second.f',
        # 'src/second_g77.f', #TODO: patch back in as needed
        # 'src/second_ifc.f', #TODO: patch back in as needed
        'src/sgutil.f',
        'src/spline.f',
        'src/userio.f',
        'src/eispack.f',
)

avl_c_wrapper = custom_target('libavlmodule.c',
input : ['src/f2py/libavl.pyf'],
output : ['libavlmodule.c', 'libavl-f2pywrappers.f'],
command: [py3, '-m', 'numpy.f2py',  '@INPUT@',
'--lower', ]
)

py3.extension_module('libavl',
    avl_source_files,
    avl_c_wrapper,
    fortranobject_c,
    include_directories: [inc_np, inc_f2py],
    dependencies: [lapack, blas, fortranobject_dep],
    subdir: 'pyavl',
    install : true,
    link_language: 'fortran',
    fortran_args: ['-ffixed-line-length-80', 
                   '-std=legacy',
                   '-Wno-align-commons',
                   '-Wmaybe-uninitialized', 
                   '-O2',
                   '-fdefault-real-8',
                   '-fPIC']
                   
)

py3.install_sources([
  'pyavl/__init__.py',
  'pyavl/pyAVL.py',
  'pyavl/MExt.py',
],
  subdir: 'pyavl'
)


        
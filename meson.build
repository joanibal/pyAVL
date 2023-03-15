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

is_windows = host_machine.system() == 'windows'

is_mingw = is_windows and cc.get_id() == 'gcc'

if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_project_link_arguments(gcc_link_args, language: ['c', 'cpp'])
    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_project_arguments('-mlong-double-64', language: 'c')
    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_project_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command('_build_utils/gcc_build_bitness.py').stdout().strip()
    if bitness == '64'
      add_project_arguments('-DMS_WIN64', language: ['c', 'cpp', 'fortran'])
    endif
  endif
  if meson.get_compiler('fortran').get_id() == 'gcc'
    add_project_link_arguments(gcc_link_args, language: ['fortran'])
    # Flag needed to work around BLAS and LAPACK Gfortran dependence on
    # undocumented C feature when passing single character string
    # arguments.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
    # https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
    add_project_arguments('-fno-optimize-sibling-calls', language: ['fortran'])
  endif
endif


_linker_script = meson.project_source_root() / '_build_utils/link-version-pyinit.map'
version_link_args = ['-Wl,--version-script=' + _linker_script]

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

# Fortran warning flags

_fflag_Wno_maybe_uninitialized = ff.get_supported_arguments('-Wno-maybe-uninitialized')
_fflag
# The default list of warnings to ignore from Fortran code. There is a lot of
# old, vendored code that is very bad and we want to compile it silently (at
# least with GCC and Clang)
ff_args = ff.get_supported_arguments(
 _fflag_Wno_maybe_uninitialized,
)

if is_windows
  ff_args += ff.get_supported_arguments('-fomit-frame-pointer')
endif

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
    link_args: version_link_args,
    fortran_args: ff_args,
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


        
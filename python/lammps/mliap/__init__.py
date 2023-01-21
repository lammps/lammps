
# Check compatiblity of this build with the python shared library.
# If this fails, lammps will segfault because its library will
# try to improperly start up a new interpreter.
import sysconfig
import ctypes
import platform
import warnings

py_ver = sysconfig.get_config_vars('VERSION')[0]
OS_name = platform.system()

if OS_name == "Darwin":
    SHLIB_SUFFIX = '.dylib'
    library = 'libpython' + py_ver + SHLIB_SUFFIX
elif OS_name == "Windows":
    SHLIB_SUFFIX = '.dll'
    library = 'python' + py_ver + SHLIB_SUFFIX
else:
    SHLIB_SUFFIX = '.so'
    library = 'libpython' + py_ver + SHLIB_SUFFIX

try:
    pylib = ctypes.CDLL(library)
except Exception as e:
    raise OSError("Unable to locate python shared library") from e

if not pylib.Py_IsInitialized():
    warnings.warn("This interpreter is not compatible with python-based MLIAP for LAMMPS. "
                  "Attempting to activate the MLIAP-python coupling from python may result "
                  "in undefined behavior.")
else:
    from .loader import load_model, load_unified, activate_mliappy
    try:
         from .loader import  load_model_kokkos,  activate_mliappy_kokkos
    except Exception as ee:
        # ignore import error, it means that the KOKKOS package was not included in LAMMPS
        pass
del sysconfig, ctypes, library, pylib

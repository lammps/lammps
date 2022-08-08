
# Check compatiblity of this build with the python shared library.
# If this fails, lammps will segfault because its library will
# try to improperly start up a new interpreter.
import sysconfig
import ctypes
import platform

library_dir = sysconfig.get_config_vars('LIBDIR')[0]
library_name = sysconfig.get_config_vars('LIBRARY')[0]
library = library_dir + "/" + library_name

OS_name = platform.system()

if OS_name == "Linux":
    SHLIB_SUFFIX = '.so'
elif OS_name == "Darwin":
    SHLIB_SUFFIX = '.dylib'
elif OS_name == "Windows":
    SHLIB_SUFFIX = '.dll'
else:
    SHLIB_SUFFIX = sysconfig.get_config_vars('SHLIB_SUFFIX')

try:
    pylib = ctypes.CDLL(library)
except OSError as e:
    if library.endswith(".a"):
        library = library.strip(".a") + ".so"
        pylib = ctypes.CDLL(library)
    else:
        raise e

if not pylib.Py_IsInitialized():
    raise RuntimeError("This interpreter is not compatible with python-based mliap for LAMMPS.")

del sysconfig, ctypes, library, pylib

from .loader import load_model, activate_mliappy

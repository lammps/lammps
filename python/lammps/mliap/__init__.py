
# Check compatiblity of this build with the python shared library.
# If this fails, lammps will segfault because its library will
# try to improperly start up a new interpreter.
import sysconfig
import ctypes
import platform

library_dir = sysconfig.get_config_vars('prefix')[0]
py_ver = sysconfig.get_config_vars('VERSION')[0]
OS_name = platform.system()

if OS_name == "Linux":
    SHLIB_SUFFIX = '.so'
    library = library_dir + '/lib/' + 'libpython' + py_ver + SHLIB_SUFFIX
elif OS_name == "Darwin":
    SHLIB_SUFFIX = '.dylib'
    library = library_dir + '/lib/' + 'libpython' + py_ver + SHLIB_SUFFIX
elif OS_name == "Windows":
    SHLIB_SUFFIX = '.dll'
    library = library_dir + '\\' + 'python' + py_ver + SHLIB_SUFFIX

try:
    pylib = ctypes.CDLL(library)
except:
    OSError

if not pylib.Py_IsInitialized():
    raise RuntimeError("This interpreter is not compatible with python-based mliap for LAMMPS.")

del sysconfig, ctypes, library, pylib

from .loader import load_model, activate_mliappy

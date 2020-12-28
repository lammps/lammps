
# Check compatiblity of this build with the python shared library.
# If this fails, lammps will segfault because its library will
# try to improperly start up a new interpreter.
import sysconfig
import ctypes
library = sysconfig.get_config_vars('INSTSONAME')[0]
pylib = ctypes.CDLL(library)
if not pylib.Py_IsInitialized():
    raise RuntimeError("This interpreter is not compatible with python-based mliap for LAMMPS.")
del sysconfig, ctypes, library, pylib

from .loader import load_model, activate_mliappy

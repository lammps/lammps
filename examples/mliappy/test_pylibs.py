import sysconfig, os,ctypes

library = sysconfig.get_config_vars('INSTSONAME')[0]

pylib = ctypes.CDLL(library)

connected = pylib.Py_IsInitialized()

if not connected:
      print("FAILURE: This interpreter is not compatible with python-driven mliappy.")
else:
      print("SUCCESS: This interpreter is compatible with python-driven MLIAPPY")

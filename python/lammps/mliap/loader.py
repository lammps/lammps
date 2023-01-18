# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

# ----------------------------------------------------------------------
#   Contributing author: Nicholas Lubbers (LANL)
# -------------------------------------------------------------------------


import sys
import importlib
import importlib.util
import importlib.machinery
import importlib.abc

from ctypes import pythonapi, c_int, c_void_p, py_object

# This dynamic loader imports a python module embedded in a shared library.
# The default value of api_version is 1013 because it has been stable since 2006.
class DynamicLoader(importlib.abc.Loader):
    def __init__(self,module_name,library,api_version=1013):
        self.api_version = api_version

        attr = "PyInit_"+module_name
        initfunc = getattr(library,attr)
        # c_void_p is standin for PyModuleDef *
        initfunc.restype = c_void_p
        initfunc.argtypes = ()
        self.module_def = initfunc()

    def create_module(self, spec):
        createfunc = pythonapi.PyModule_FromDefAndSpec2
        # c_void_p is standin for PyModuleDef *
        createfunc.argtypes = c_void_p, py_object, c_int
        createfunc.restype = py_object
        module = createfunc(self.module_def, spec, self.api_version)
        return module

    def exec_module(self, module):
        execfunc = pythonapi.PyModule_ExecDef
        # c_void_p is standin for PyModuleDef *
        execfunc.argtypes = py_object, c_void_p
        execfunc.restype = c_int
        result = execfunc(module, self.module_def)
        if result<0:
           raise ImportError()

def activate_mliappy(lmp):
    try:
        library = lmp.lib
        module_names = ["mliap_model_python_couple", "mliap_unified_couple"]
        api_version = library.lammps_python_api_version()

        for module_name in module_names:
            # Make Machinery
            loader = DynamicLoader(module_name,library,api_version)
            spec = importlib.util.spec_from_loader(module_name,loader)

            # Do the import
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
    except Exception as ee:
        raise ImportError("Could not load ML-IAP python coupling module.") from ee

def activate_mliappy_kokkos(lmp):
    try:
        library = lmp.lib
        module_names = ["mliap_model_python_couple_kokkos"]
        api_version = library.lammps_python_api_version()

        for module_name in module_names:
            # Make Machinery
            loader = DynamicLoader(module_name,library,api_version)
            spec = importlib.util.spec_from_loader(module_name,loader)

            # Do the import
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
    except Exception as ee:
        raise ImportError("Could not load ML-IAP python coupling module.") from ee

def load_model(model):
    try:
        import mliap_model_python_couple
    except ImportError as ie:
        raise ImportError("ML-IAP python module must be activated before loading\n"
                          "the pair style. Call lammps.mliap.activate_mliappy(lmp)."
                          ) from ie
    mliap_model_python_couple.load_from_python(model)

def load_model_kokkos(model):
    try:
        import mliap_model_python_couple_kokkos
    except ImportError as ie:
        raise ImportError("ML-IAP python module must be activated before loading\n"
                          "the pair style. Call lammps.mliap.activate_mliappy(lmp)."
                          ) from ie
    mliap_model_python_couple_kokkos.load_from_python(model)


def load_unified(model):
    try:
        import mliap_unified_couple
    except ImportError as ie:
        raise ImportError("ML-IAP python module must be activated before loading\n"
                          "the pair style. Call lammps.mliap.activate_mliappy(lmp)."
                          ) from ie
    mliap_unified_couple.load_from_python(model)


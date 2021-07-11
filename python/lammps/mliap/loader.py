# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
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
import importlib.util
import importlib.machinery

def activate_mliappy(lmp):
    try:
        # Begin Importlib magic to find the embedded python module
        # This is needed because the filename for liblammps does not
        # match the spec for normal python modules, wherein
        # file names match with PyInit function names.
        # Also, python normally doesn't look for extensions besides '.so'
        # We fix both of these problems by providing an explict
        # path to the extension module 'mliap_model_python_couple' in

        path = lmp.lib._name
        loader = importlib.machinery.ExtensionFileLoader('mliap_model_python_couple', path)
        spec = importlib.util.spec_from_loader('mliap_model_python_couple', loader)
        module = importlib.util.module_from_spec(spec)
        sys.modules['mliap_model_python_couple'] = module
        spec.loader.exec_module(module)
        # End Importlib magic to find the embedded python module

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


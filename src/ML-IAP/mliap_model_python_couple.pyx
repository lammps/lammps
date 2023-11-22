# cython: language_level=3
# distutils: language = c++

cimport cython

import pickle

# For converting C arrays to numpy arrays
import numpy as np

# For converting void * to integer for tracking object identity
from libc.stdint cimport uintptr_t

from libcpp.string cimport string


cdef extern from "mliap_data.h" namespace "LAMMPS_NS":
    cdef cppclass MLIAPData:
        # Array shapes
        int nlistatoms
        int nlocal
        int ndescriptors

        # Input data
        int * ielems                # types for all atoms in list
        double ** descriptors       # descriptors for all atoms in list

        # Output data to write to
        double ** betas             # betas for all atoms in list
        double * eatoms             # energy for all atoms in list
        double energy

cdef extern from "mliap_model_python.h" namespace "LAMMPS_NS":
    cdef cppclass MLIAPModelPython:
        void connect_param_counts()


class MLIAPPYModelNotLinked(Exception): pass


LOADED_MODELS = {}

cdef object c_id(MLIAPModelPython * c_model):
    """
    Use python-style id of object to keep track of identity.
    Note, this is probably not a perfect general strategy but it should work fine with LAMMPS pair styles.
    """
    return int(<uintptr_t> c_model)

cdef object retrieve(MLIAPModelPython * c_model) with gil:
    try:
        model = LOADED_MODELS[c_id(c_model)]
    except KeyError as ke:
        raise KeyError("Model has not been loaded.") from ke
    if model is None:
        raise MLIAPPYModelNotLinked("Model not linked, connect the model from the python side.")
    return model

cdef public int MLIAPPY_load_model(MLIAPModelPython * c_model, char* fname) with gil:
    str_fname = fname.decode('utf-8') # Python 3 only; not Python 2 not supported.
    if str_fname == "LATER":
        model = None
        returnval = 0
    else:
        if str_fname.endswith(".pt") or str_fname.endswith('.pth'):
            import torch
            model = torch.load(str_fname)
        else:
            with open(str_fname,'rb') as pfile:
                model = pickle.load(pfile)
        returnval = 1
    LOADED_MODELS[c_id(c_model)] = model
    return returnval

def load_from_python(model):
    unloaded_models = [k for k, v in LOADED_MODELS.items() if v is None]
    num_models = len(unloaded_models)
    cdef MLIAPModelPython * lmp_model

    if num_models == 0:
            raise ValueError("No model in the waiting area.")
    elif num_models > 1:
            raise ValueError("Model is amibguous, more than one model in waiting area.")
    else:
        c_id = unloaded_models[0]
        LOADED_MODELS[c_id]=model
        lmp_model = <MLIAPModelPython *> <uintptr_t> c_id
        lmp_model.connect_param_counts()


cdef public void MLIAPPY_unload_model(MLIAPModelPython * c_model) with gil:
    del LOADED_MODELS[c_id(c_model)]

cdef public int MLIAPPY_nparams(MLIAPModelPython * c_model) with gil:
    return int(retrieve(c_model).n_params)

cdef public int MLIAPPY_nelements(MLIAPModelPython * c_model) with gil:
    return int(retrieve(c_model).n_elements)

cdef public int MLIAPPY_ndescriptors(MLIAPModelPython * c_model) with gil:
    return int(retrieve(c_model).n_descriptors)

cdef public void MLIAPPY_compute_gradients(MLIAPModelPython * c_model, MLIAPData * data) with gil:
    model = retrieve(c_model)

    n_d = data.ndescriptors
    n_a = data.nlistatoms

    # Make numpy arrays from pointers
    beta_np = np.asarray(<double[:n_a,:n_d] > &data.betas[0][0])
    desc_np = np.asarray(<double[:n_a,:n_d]> &data.descriptors[0][0])
    elem_np = np.asarray(<int[:n_a]> &data.ielems[0])
    en_np = np.asarray(<double[:n_a]> &data.eatoms[0])

    # Invoke python model on numpy arrays.
    model(elem_np,desc_np,beta_np,en_np)

    # Get the total energy from the atom energy.
    energy = np.sum(en_np)
    data.energy = <double> energy
    return

# cython: language_level=3
# distutils: language = c++

cimport cython

import pickle

# For converting C arrays to numpy arrays
import numpy as np
cimport numpy as cnp

# For converting void * to integer for tracking object identity
from libc.stdint cimport uintptr_t

from libcpp.string cimport string


cdef extern from "mliap_data.h" namespace "LAMMPS_NS":
    cdef cppclass MLIAPData:
        # Array shapes
        int natoms
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
        ctypedef void (*CBPtr)(void * , MLIAPData);
        void set_model(CBPtr, void *);
     

LOADED_MODELS = {}
cdef public int MLIAPPY_load_model(MLIAPModelPython * c_model, char* fname) except 1 with gil:
    str_fname = fname.decode('utf-8') # Python 3 only; not Python 2 not supported.
    
    with open(str_fname,'rb') as pfile:
            model = pickle.load(pfile)

    LOADED_MODELS[int(<uintptr_t> c_model)] = model
    return 0

cdef public void MLIAPPY_unload_model(MLIAPModelPython * c_model) with gil:
    del LOADED_MODELS[int(<uintptr_t> c_model)]
    
cdef public int MLIAPPY_nparams(MLIAPModelPython * c_model) with gil:
    model = LOADED_MODELS[int(<uintptr_t> c_model)]
    n_params = int(model.n_params)
    return <int> n_params

cdef public int MLIAPPY_nelements(MLIAPModelPython * c_model) with gil:
    model = LOADED_MODELS[int(<uintptr_t> c_model)]
    n_elements = int(model.n_elements)
    return <int> n_elements

cdef public int MLIAPPY_ndescriptors(MLIAPModelPython * c_model) with gil:
    model = LOADED_MODELS[int(<uintptr_t> c_model)]
    n_descriptors = int(model.n_descriptors)
    return <int> n_descriptors

cdef public MLIAPPY_model_callback(MLIAPModelPython * c_model, MLIAPData * data) with gil:
    model = LOADED_MODELS[int(<uintptr_t> c_model)]

    n_d = data.ndescriptors
    n_a = data.natoms
    
    # Make numpy arrays from pointers
    beta_np = np.asarray(<double[:n_a,:n_d] > &data.betas[0][0])
    desc_np = np.asarray(<double[:n_a,:n_d]> &data.descriptors[0][0])
    type_np = np.asarray(<int[:n_a]> &data.ielems[0])
    en_np = np.asarray(<double[:n_a]> &data.eatoms[0])

    # Invoke python model on numpy arrays.
    model(type_np,desc_np,beta_np,en_np)
    
    # Get the total energy from the atom energy.
    energy = np.sum(en_np)
    data.energy = <double> energy
    return

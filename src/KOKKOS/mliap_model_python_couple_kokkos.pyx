# cython: language_level=3
# distutils: language = c++

cimport cython

from libc.stdint cimport uintptr_t

import pickle

# For converting C arrays to numpy arrays
import numpy
import torch
try:
    import cupy
except ImportError:
    pass

# For converting void * to integer for tracking object identity
from libc.stdint cimport uintptr_t

from libcpp.string cimport string


cdef extern from "mliap_data_kokkos.h" namespace "LAMMPS_NS":
    cdef cppclass MLIAPDataKokkosDevice:
        # Array shapes
        int nlistatoms
        int nlocal
        int ndescriptors

        # Input data
        int * ielems                # types for all atoms in list
        double * descriptors       # descriptors for all atoms in list

        # Output data to write to
        double * betas             # betas for all atoms in list
        double * eatoms             # energy for all atoms in list
        double *energy
        int dev
cdef extern from "mliap_model_python_kokkos.h" namespace "LAMMPS_NS":
    cdef cppclass MLIAPModelPythonKokkosDevice:
        void connect_param_counts()


class MLIAPPYKokkosModelNotLinked(Exception): pass


LOADED_MODELS = {}

cdef object c_id(MLIAPModelPythonKokkosDevice * c_model):
    """
    Use python-style id of object to keep track of identity.
    Note, this is probably not a perfect general strategy but it should work fine with LAMMPS pair styles.
    """
    return int(<uintptr_t> c_model)

cdef object retrieve(MLIAPModelPythonKokkosDevice * c_model) with gil:
    try:
        model = LOADED_MODELS[c_id(c_model)]
    except KeyError as ke:
        raise KeyError("Model has not been loaded.") from ke
    if model is None:
        raise MLIAPPYKokkosModelNotLinked("Model not linked, connect the model from the python side.")
    return model

cdef public int MLIAPPYKokkos_load_model(MLIAPModelPythonKokkosDevice * c_model, char* fname) with gil:
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
    cdef MLIAPModelPythonKokkosDevice * lmp_model

    if num_models == 0:
            raise ValueError("No model in the waiting area.")
    elif num_models > 1:
            raise ValueError("Model is amibguous, more than one model in waiting area.")
    else:
        c_id = unloaded_models[0]
        LOADED_MODELS[c_id]=model
        lmp_model = <MLIAPModelPythonKokkosDevice *> <uintptr_t> c_id
        lmp_model.connect_param_counts()


cdef public void MLIAPPYKokkos_unload_model(MLIAPModelPythonKokkosDevice * c_model) with gil:
    del LOADED_MODELS[c_id(c_model)]

cdef public int MLIAPPYKokkos_nparams(MLIAPModelPythonKokkosDevice * c_model) with gil:
    return int(retrieve(c_model).n_params)

cdef public int MLIAPPYKokkos_nelements(MLIAPModelPythonKokkosDevice * c_model) with gil:
    return int(retrieve(c_model).n_elements)

cdef public int MLIAPPYKokkos_ndescriptors(MLIAPModelPythonKokkosDevice * c_model) with gil:
    return int(retrieve(c_model).n_descriptors)

cdef create_array(device, void *pointer, shape,is_int):

    size=1
    for i in shape:
        size = size*i

    if ( device == 1):
        mem = cupy.cuda.UnownedMemory(ptr=int( <uintptr_t> pointer), owner=None, size=size)
        memptr = cupy.cuda.MemoryPointer(mem, 0)
        type=cupy.double
        if (is_int):
            type=cupy.int32
        return cupy.ndarray(shape, type, memptr=memptr)
    else:
        if (len(shape) == 1 ):
            if (is_int):
                return numpy.asarray(<int[:shape[0]]>pointer)
            else:
                return numpy.asarray(<double[:shape[0]]>pointer)
        else:
            if (is_int):
                return numpy.asarray(<int[:shape[0],:shape[1]]>pointer)
            else:
                return numpy.asarray(<double[:shape[0],:shape[1]]>pointer)


cdef public void MLIAPPYKokkos_compute_gradients(MLIAPModelPythonKokkosDevice * c_model, MLIAPDataKokkosDevice * data) with gil:

    dev=data.dev

    torch.cuda.nvtx.range_push("set data fields")
    model = retrieve(c_model)
    n_d = data.ndescriptors
    n_a = data.nlistatoms

    cdef void* ptr = data.ielems
    # Make numpy arrays from pointers
    elem_cp = create_array(dev, data.ielems, (n_d,), True)
    en_cp   = create_array(dev, data.eatoms, (n_a,), False)
    beta_cp = create_array(dev, data.betas, (n_a, n_d), False)
    desc_cp = create_array(dev, data.descriptors, (n_a, n_d), False)
    torch.cuda.nvtx.range_pop()

    # Invoke python model on numpy arrays.
    torch.cuda.nvtx.range_push("call model")
    model(elem_cp,desc_cp,beta_cp,en_cp,dev==1)
    torch.cuda.nvtx.range_pop()

    # Get the total energy from the atom energy.
    energy = cupy.sum(en_cp)
    data.energy[0] = <double> energy
    return

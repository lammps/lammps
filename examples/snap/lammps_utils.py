import numpy as np
import ctypes

def set_cmdlinevars(cmdargs, argdict):
    for key in argdict.keys():
        cmdargs += ["-var",key,f"{argdict[key]}"]
    return cmdargs

def extract_commands(string):
    return [x for x in string.splitlines() if x.strip() != '']

def extract_compute_np(lmp,name,compute_type,result_type,array_shape):
    """                                                                                                               
    Convert a lammps compute to a numpy array.                                                                        
    Assumes the compute returns a floating point numbers.                                                             
    Note that the result is a view into the original memory.                                                          
    If the result type is 0 (scalar) then conversion to numpy is skipped and a python float is returned.              
    """
    ptr = lmp.extract_compute(name, compute_type, result_type)  # 1,2: Style (1) is per-atom compute, returns array type (2).
    if result_type == 0: return ptr # No casting needed, lammps.py already works                                      
    if result_type == 2: ptr = ptr.contents
    total_size = np.prod(array_shape)
    buffer_ptr = ctypes.cast(ptr, ctypes.POINTER(ctypes.c_double * total_size))
    array_np = np.frombuffer(buffer_ptr.contents, dtype=float)
    array_np.shape = array_shape
    return array_np

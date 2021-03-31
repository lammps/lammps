
from __future__ import print_function

def square(val):
    return val*val

def printnum():
    print("2.25")

def printtxt():
    print("sometext")

def getidxvar(lmpptr):
    from lammps import lammps, LMP_VAR_EQUAL
    lmp = lammps(ptr=lmpptr)

    val = lmp.extract_variable("idx",None,LMP_VAR_EQUAL)
    print(val)

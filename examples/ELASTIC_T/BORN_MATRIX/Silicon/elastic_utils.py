import numpy as np
from lammps import lammps, LMP_VAR_EQUAL

# method for rotating elastic constants

def rotate_cij(cij, r):

    # K_1
    # R_11^2 R_12^2 R_13^2
    # R_21^2 R_22^2 R_23^2
    # R_31^2 R_32^2 R_33^2

    k1 = r*r

    # K_2
    # R_12.R_13 R_13.R_11 R_11.R_12
    # R_22.R_23 R_23.R_21 R_21.R_22
    # R_32.R_33 R_33.R_31 R_31.R_32

    k2 = np.array([
        [r[0][1]*r[0][2], r[0][2]*r[0][0], r[0][0]*r[0][1]],
        [r[1][1]*r[1][2], r[1][2]*r[1][0], r[1][0]*r[1][1]],
        [r[2][1]*r[2][2], r[2][2]*r[2][0], r[2][0]*r[2][1]],
    ])

    # K_3
    # R_21.R_31 R_22.R_32 R_23.R_33
    # R_31.R_11 R_32.R_12 R_33.R_13
    # R_11.R_21 R_12.R_22 R_13.R_23

    k3 = np.array([
        [r[1][0]*r[2][0], r[1][1]*r[2][1], r[1][2]*r[2][2]],
        [r[2][0]*r[0][0], r[2][1]*r[0][1], r[2][2]*r[0][2]],
        [r[0][0]*r[1][0], r[0][1]*r[1][1], r[0][2]*r[1][2]],
    ])

    # K_4a
    # R_22.R_33 R_23.R_31 R_21.R_32
    # R_32.R_13 R_33.R_11 R_31.R_12
    # R_12.R_23 R_13.R_21 R_11.R_22

    k4a = np.array([
        [r[1][1]*r[2][2], r[1][2]*r[2][0], r[1][0]*r[2][1]],
        [r[2][1]*r[0][2], r[2][2]*r[0][0], r[2][0]*r[0][1]],
        [r[0][1]*r[1][2], r[0][2]*r[1][0], r[0][0]*r[1][1]],
    ])

    # K_4b
    # R_23.R_32 R_21.R_33 R_22.R_31
    # R_33.R_12 R_31.R_13 R_32.R_11
    # R_13.R_22 R_11.R_23 R_12.R_21

    k4b = np.array([
        [r[1][2]*r[2][1], r[1][0]*r[2][2], r[1][1]*r[2][0]],
        [r[2][2]*r[0][1], r[2][0]*r[0][2], r[2][1]*r[0][0]],
        [r[0][2]*r[1][1], r[0][0]*r[1][2], r[0][1]*r[1][0]],
    ])
    
    k = np.block([[k1, 2*k2],[k3, k4a+k4b]])
    cijrot = k.dot(cij.dot(k.T))
    return cijrot

def calculate_cij(cmdargs):
    lmp = lammps(cmdargs=cmdargs)
    lmp.file("in.elastic")

    C11 = lmp.extract_variable("C11",None, LMP_VAR_EQUAL)
    C22 = lmp.extract_variable("C22",None, LMP_VAR_EQUAL)
    C33 = lmp.extract_variable("C33",None, LMP_VAR_EQUAL)
    C44 = lmp.extract_variable("C44",None, LMP_VAR_EQUAL)
    C55 = lmp.extract_variable("C55",None, LMP_VAR_EQUAL)
    C66 = lmp.extract_variable("C66",None, LMP_VAR_EQUAL)
    
    C12 = lmp.extract_variable("C12",None, LMP_VAR_EQUAL)
    C13 = lmp.extract_variable("C13",None, LMP_VAR_EQUAL)
    C14 = lmp.extract_variable("C14",None, LMP_VAR_EQUAL)
    C15 = lmp.extract_variable("C15",None, LMP_VAR_EQUAL)
    C16 = lmp.extract_variable("C16",None, LMP_VAR_EQUAL)
    
    C23 = lmp.extract_variable("C23",None, LMP_VAR_EQUAL)
    C24 = lmp.extract_variable("C24",None, LMP_VAR_EQUAL)
    C25 = lmp.extract_variable("C25",None, LMP_VAR_EQUAL)
    C26 = lmp.extract_variable("C26",None, LMP_VAR_EQUAL)
    
    C34 = lmp.extract_variable("C34",None, LMP_VAR_EQUAL)
    C35 = lmp.extract_variable("C35",None, LMP_VAR_EQUAL)
    C36 = lmp.extract_variable("C36",None, LMP_VAR_EQUAL)
    
    C45 = lmp.extract_variable("C45",None, LMP_VAR_EQUAL)
    C46 = lmp.extract_variable("C46",None, LMP_VAR_EQUAL)
    
    C56 = lmp.extract_variable("C56",None, LMP_VAR_EQUAL)

    cij = np.array([
      [C11, C12, C13, C14, C15, C16],
      [  0, C22, C23, C24, C25, C26],
      [  0,   0, C33, C34, C35, C36],
      [  0,   0,   0, C44, C45, C46],
      [  0,   0,   0,   0, C55, C56],
      [  0,   0,   0,   0,   0, C66],
    ])
    cij = np.triu(cij) + np.tril(cij.T, -1)            

    return cij

def gen_varargs(varlist):
    cmdargs = []
    for key in varlist:
        cmdargs += ["-var",key,str(varlist[key])]
    return cmdargs

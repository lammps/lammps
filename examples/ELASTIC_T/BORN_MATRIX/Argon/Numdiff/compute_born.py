#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def reduce_Born(Cf):
    C = np.zeros((6,6), dtype=np.float64)
    C[0,0] = Cf[0]
    C[1,1] = Cf[1]
    C[2,2] = Cf[2]
    C[3,3] = Cf[3]
    C[4,4] = Cf[4]
    C[5,5] = Cf[5]
    C[0,1] = Cf[6]
    C[0,2] = Cf[7]
    C[0,3] = Cf[8]
    C[0,4] = Cf[9]
    C[0,5] = Cf[10]
    C[1,2] = Cf[11]
    C[1,3] = Cf[12]
    C[1,4] = Cf[13]
    C[1,5] = Cf[14]
    C[2,3] = Cf[15]
    C[2,4] = Cf[16]
    C[2,5] = Cf[17]
    C[3,4] = Cf[18]
    C[3,5] = Cf[19]
    C[4,5] = Cf[20]
    C = np.where(C,C,C.T)
    return C

def compute_delta():
    D = np.zeros((3,3,3,3))
    for a in range(3):
        for b in range(3):
            for m in range(3):
                for n in range(3):
                    D[a,b,m,n] = (a==m)*(b==n) + (a==n)*(b==m)
    d = np.zeros((6,6))
    d[0,0] = D[0,0,0,0]
    d[1,1] = D[1,1,1,1]
    d[2,2] = D[2,2,2,2]
    d[3,3] = D[1,2,1,2]
    d[4,4] = D[0,2,0,2]
    d[5,5] = D[0,1,0,1]
    d[0,1] = D[0,0,1,1]
    d[0,2] = D[0,0,2,2]
    d[0,3] = D[0,0,1,2]
    d[0,4] = D[0,0,0,2]
    d[0,5] = D[0,0,0,1]
    d[1,2] = D[1,1,2,2]
    d[1,3] = D[1,1,1,2]
    d[1,4] = D[1,1,0,2]
    d[1,5] = D[1,1,0,1]
    d[2,3] = D[2,2,1,2]
    d[2,4] = D[2,2,0,2]
    d[2,5] = D[2,2,0,1]
    d[3,4] = D[1,2,0,2]
    d[3,5] = D[1,2,0,1]
    d[4,5] = D[0,2,0,1]
    d = np.where(d,d,d.T)
    return d


def write_matrix(C, filename):
    with open(filename, 'w') as f:
        f.write("# Cij Matrix from post process computation\n")
        for i in C:
            f.write("{:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f}\n".format(
                i[0]*10**-9, i[1]*10**-9, i[2]*10**-9, i[3]*10**-9, i[4]*10**-9, i[5]*10**-9,
                )
            )
    return

def main():

    N = 500
    vol = 27.047271**3 * 10**-30  # m^3
    T = 60  # K
    kb = 1.380649 * 10**-23  # J/K
    kbT = T*kb  # J
    kcalmol2J = 4183.9954/(6.022*10**23)

    born = np.loadtxt('born.out')
    stre = np.loadtxt('vir.out')
    stre[:, 1:] = -stre[:, 1:]*101325  # -> Pa
    try:
        mean_born = np.mean(born[:, 1:], axis=0)
    except IndexError:
        mean_born = born[1:]

    CB = kcalmol2J/vol*reduce_Born(mean_born)  # -> J/m^3=Pa
    Cs = vol/kbT*np.cov(stre[:,1:].T)
    Ct = N*kbT/vol * compute_delta()
    C = CB - Cs + Ct
    write_matrix(CB, 'born_matrix.out')
    write_matrix(Cs, 'stre_matrix.out')
    write_matrix(Ct, 'temp_matrix.out')
    write_matrix(C, 'full_matrix.out')
    C11 = np.mean([C[0,0], C[1,1], C[2,2]]) * 10**-9
    C12 = np.mean([C[0,1], C[0,2], C[1,2]]) * 10**-9
    C44 = np.mean([C[3,3], C[4,4], C[5,5]]) * 10**-9
    eC11 = np.std([C[0,0], C[1,1], C[2,2]]) * 10**-9
    eC12 = np.std([C[0,1], C[0,2], C[1,2]]) * 10**-9
    eC44 = np.std([C[3,3], C[4,4], C[5,5]]) * 10**-9
    print(C*10**-9)
    print("C11 = {:f} ± {:f}; C12 = {:f} ± {:f}; C44 = {:f} ± {:f}".format(C11, eC11, C12, eC12, C44, eC44))

    return


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit("User interruption.")


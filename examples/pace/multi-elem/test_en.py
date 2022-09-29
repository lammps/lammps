from __future__ import print_function
import sys, os
import ctypes
import numpy as np
from ase.io import read,write
from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

# get MPI settings from LAMMPS
def run_struct(f):
    file_prefix = f.split('.')[0]
    atoms = read(f)
    lmp = lammps()

    me = lmp.extract_setting("world_rank")
    nprocs = lmp.extract_setting("world_size")

    write('%s.data' % file_prefix,atoms,format='lammps-data')

    cmds = ["-screen", "none", "-log", "none"]
    lmp = lammps(cmdargs = cmds)

    print("Made LAMMPS instance")

    def run_lammps(dgradflag):

        # simulation settings
        fname = file_prefix
        lmp.command("clear")
        lmp.command("info all out log")
        lmp.command('units  metal')
        lmp.command('atom_style  atomic')
        lmp.command("boundary    p p p")
        lmp.command("atom_modify    map hash")
        lmp.command('neighbor  2.3 bin')
        # boundary
        lmp.command('boundary  p p p')
        # read atoms
        lmp.command('read_data  %s.data' % fname )
        lmp.command('mass  1 1.00')
        lmp.command('mass  2 14.00')
        lmp.command('mass  3 15.999')

        # potential settings

        lmp.command(f"pair_style     zero 5.7")
        lmp.command(f"pair_coeff     * *")


        if dgradflag:
            lmp.command(f"compute     pace all pace coupling_coefficients.yace 1 1")
        else:
            lmp.command(f"compute     pace all pace coupling_coefficients.yace 1 0 ")

        # run

        lmp.command(f"thermo         100")
        lmp.command(f"run {nsteps}")

    # declare simulation/structure variables

    nsteps = 0
    ntypes = 3

    # declare compute pace variables


    bikflag = 1

    # NUMBER of descriptors
    nd = 91

    dgradflag = 0
    
    run_lammps(dgradflag)
    
    lmp_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
    print ('global shape',np.shape(lmp_pace))
    np.save('%s_chi_i.npy' % file_prefix,lmp_pace)
    lmp.close()
    del lmp
    return None

import glob
for f in sorted(glob.glob('*.xyz')):
    print ('running %s' % f)
    run_struct(f)

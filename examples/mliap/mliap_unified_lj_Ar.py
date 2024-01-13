# Demonstrate how to load a unified model from python.
# This is essentially the same as in.mliap.unified.lj.Ar
# except that python is the driving program, and lammps
# is in library mode.

before_loading =\
"""# 3d Lennard-Jones melt

units       lj
atom_style  atomic

lattice     fcc 0.8442
region      box block 0 10 0 10 0 10
create_box  1 box
create_atoms    1 box
mass        1 1.0

velocity    all create 3.0 87287 loop geom
"""
after_loading =\
"""

pair_style  mliap unified EXISTS
pair_coeff  * * Ar

neighbor    0.3 bin
neigh_modify    every 20 delay 0 check no

fix     1 all nve

thermo      50
run     250
"""

import lammps

lmp = lammps.lammps(cmdargs=['-echo','both'])

# Before defining the pair style, you must do the following:
import lammps.mliap
lammps.mliap.activate_mliappy(lmp)
# Otherwise, when running lammps in library mode,
# you will get an error

# Setup the simulation to just before declaring
# the mliap unified pair style
lmp.commands_string(before_loading)

# Define the model however you like. In this example
# we simply import the unified L-J example from mliap
from lammps.mliap.mliap_unified_lj import MLIAPUnifiedLJ
unified = MLIAPUnifiedLJ(["Ar"])

# You can also load the model from a pickle file.
# import pickle
# with open('mliap_unified_lj_Ar.pkl', 'rb') as pfile:
#     unified = pickle.load(pfile)

# Connect the L-J model to the mliap unified pair style.
lammps.mliap.load_unified(unified)


# Run the simulation with the mliap unified pair style
# Use pre-loaded model by specifying model filename as "EXISTS"
lmp.commands_string(after_loading)
lmp.close()
lmp.finalize()

# Demonstrate how to load a model from the python side.
# This is essentially the same as in.mliap.pytorch.Ta06A
# except that python is the driving program, and lammps
# is in library mode.

before_loading =\
"""# Demonstrate MLIAP/PyTorch interface to linear SNAP potential

# Initialize simulation

variable nsteps index 100
variable nrep equal 4
variable a equal 3.316
units           metal

# generate the box and atom positions using a BCC lattice

variable nx equal ${nrep}
variable ny equal ${nrep}
variable nz equal ${nrep}

boundary        p p p

lattice         bcc $a
region          box block 0 ${nx} 0 ${ny} 0 ${nz}
create_box      1 box
create_atoms    1 box

mass 1 180.88

# choose potential

# DATE: 2014-09-05 UNITS: metal CONTRIBUTOR: Aidan Thompson athomps@sandia.gov CITATION: Thompson, Swiler, Trott, Foiles and Tucker, arxiv.org, 1409.3880 (2014)

# Definition of SNAP potential Ta_Cand06A
# Assumes 1 LAMMPS atom type

variable zblcutinner equal 4
variable zblcutouter equal 4.8
variable zblz equal 73

# Specify hybrid with SNAP, ZBL

pair_style hybrid/overlay &
zbl ${zblcutinner} ${zblcutouter} &
mliap model mliappy LATER &
descriptor sna Ta06A.mliap.descriptor
pair_coeff 1 1 zbl ${zblz} ${zblz}
pair_coeff * * mliap Ta
"""
after_loading =\
"""

# Setup output

compute  eatom all pe/atom
compute  energy all reduce sum c_eatom

compute  satom all stress/atom NULL
compute  str all reduce sum c_satom[1] c_satom[2] c_satom[3]
variable press equal (c_str[1]+c_str[2]+c_str[3])/(3*vol)

thermo_style    custom step temp epair c_energy etotal press v_press
thermo          10
thermo_modify norm yes

# Set up NVE run

timestep 0.5e-3
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes

# Run MD

velocity all create 300.0 4928459 loop geom
fix 1 all nve
run             ${nsteps}
"""

import lammps

lmp = lammps.lammps(cmdargs=['-echo','both'])

# Before defining the pair style, one must do the following:
import lammps.mliap
lammps.mliap.activate_mliappy(lmp)
# Otherwise, when running lammps in library mode,
# you will get an error:
# "ERROR: Loading MLIAPPY coupling module failure."

# Setup the simulation and declare an empty model
# by specifying model filename as "LATER"
lmp.commands_string(before_loading)

# Define the model however you like. In this example
# we load it from disk:
import os
import torch
torch_model = 'Ta06A.mliap.pytorch.model.pt'
if not os.path.exists(torch_model):
    raise FileNotFoundError(f"Generate {torch_model} with convert_mliap_Ta06A.py")
model = torch.load(torch_model)

# Connect the PyTorch model to the mliap pair style.
lammps.mliap.load_model(model)
  
# run the simulation with the mliap pair style
lmp.commands_string(after_loading)

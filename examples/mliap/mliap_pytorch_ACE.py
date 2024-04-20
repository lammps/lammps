# Demonstrate how to load a model from the python side.
# This is essentially the same as in.mliap.pytorch.MOF
# except that python is the driving program, and lammps
# is in library mode.
before_loading =\
"""# Demonstrate MLIAP/PyTorch interface to torch model

# Initialize simulation

variable        nsteps index 100
variable        nrep equal 4
variable        a equal 3.316
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

# define potential with LATER mliappy


pair_style  mliap model mliappy LATER descriptor ace ccs_single_element.yace
pair_coeff * * Ta

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
thermo_modify   norm yes

# Set up NVE run
dump            1 all cfg 10 ats.*.cfg mass type xs ys zs
dump_modify     1 element Ta

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
torch_model = 'ACE_NN_Pytorch.pt'
if not os.path.exists(torch_model):
    raise FileNotFoundError(f"Generate {torch_model} first")
model = torch.load(torch_model)

# Connect the PyTorch model to the mliap pair style.
lammps.mliap.load_model(model)

# run the simulation with the mliap pair style
lmp.commands_string(after_loading)
lmp.close()
lmp.finalize()

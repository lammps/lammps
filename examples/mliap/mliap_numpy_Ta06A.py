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

import numpy as np

class LinearModel():
  def __init__(self,file):
      coeffs = np.genfromtxt(file,skip_header=6)
      self.bias = coeffs[0]
      self.weights = coeffs[1:]
      self.n_params = len(coeffs)
      self.n_descriptors = len(self.weights)
      self.n_elements = 1
  def __call__(self,elems,bispectrum,beta,energy):
      energy[:] = bispectrum @ self.weights + self.bias
      beta[:] = self.weights

mymodel = LinearModel("Ta06A.mliap.model")

import lammps

lmp = lammps.lammps(cmdargs=['-echo','both'])

import lammps.mliap
lammps.mliap.activate_mliappy(lmp)
lmp.commands_string(before_loading)
lammps.mliap.load_model(mymodel)
lmp.commands_string(after_loading)
lmp.close()
lmp.finalize()

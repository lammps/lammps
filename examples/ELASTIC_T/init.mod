# NOTE: This script can be modified for different atomic structures, 
# units, etc. See in.elastic for more info.
#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
variable up equal 2.0e-2
 
# metal units, elastic constants in GPa
units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

# Define MD parameters
variable nevery equal 10
variable nrepeat equal 10
variable nfreq equal ${nevery}*${nrepeat}
variable nthermo equal ${nfreq}
variable nequil equal 10*${nthermo}
variable nrun equal 3*${nthermo}
variable temp equal 2000.0
variable timestep equal 0.001
variable mass1 equal 28.06
variable tdamp equal 0.01
variable seed equal 123457

# generate the box and atom positions using a diamond lattice
variable a equal 5.431

boundary	p p p

lattice         diamond $a
region		box prism 0 3.0 0 3.0 0 3.0 0.0 0.0 0.0
create_box	1 box
create_atoms	1 box
mass 1 ${mass1}


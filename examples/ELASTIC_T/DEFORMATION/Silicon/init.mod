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
variable nevery equal 10                  # sampling interval
variable nrepeat equal 10                 # number of samples
variable nfreq equal ${nevery}*${nrepeat} # length of one average
variable nthermo equal ${nfreq}           # interval for thermo output
variable nequil equal 10*${nthermo}       # length of equilibration run
variable nrun equal 3*${nthermo}          # length of equilibrated run
variable temp equal 2000.0                # temperature of initial sample
variable timestep equal 0.001             # timestep
variable mass1 equal 28.06                # mass
variable adiabatic equal 0                # adiabatic (1) or isothermal (2)
variable tdamp equal 0.01                 # time constant for thermostat
variable seed equal 123457                # seed for thermostat

# generate the box and atom positions using a diamond lattice
variable a equal 5.431

boundary	p p p

lattice         diamond $a
region		box prism 0 3.0 0 3.0 0 3.0 0.0 0.0 0.0
create_box	1 box
create_atoms	1 box
mass 1 ${mass1}
velocity	all create ${temp} 87287



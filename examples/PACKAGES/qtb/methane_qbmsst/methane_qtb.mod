## This script first constructs a liquid methane structure of a given size. It then uses fix qtb to equilibrate the computational cell to the specified temperature and pressure.


## This part defines units, methane structure, and atomic information
#General
units			real
dimension		3
boundary		p p p
atom_style		charge

#Lattice
lattice			custom 1.0 &
			a1	3.9783624 0 0 &
			a2	0 3.9783624 0 &
			a3	0 0 3.9783624 &
			&
			basis	0.5 0.5 0.5 & 
			basis	0.663 0.663 0.663 &
			basis	0.337 0.337 0.663 &
			basis	0.663 0.337 0.337 &
			basis	0.337 0.663 0.337

#Computational Cell
region			simbox block 0 3.9783624 0 3.9783624 0 3.9783624 units box
create_box		2 simbox
create_atoms		1 box & 
			basis	1 1 &
			basis	2 2 &
			basis	3 2 &
			basis	4 2 &
			basis	5 2
replicate		${x_rep} ${y_rep} ${z_rep}

#Atomic Information
mass			1 12.011150
mass			2  1.007970


## This part defines the reax pair potential in methane, force field coefficients are specified in "ffield.reax"
#Pair Potentials
pair_style		reax/c NULL
pair_coeff		* * ffield.reax C H
fix                     0 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

#Neighbor Style
neighbor		2.5 bin
neigh_modify		every 10 delay 0 check no


## This part equilibrates liquid methane to a temperature of ${temperature}(unit temperatureture) with quantum nuclear effects
#Initialization
velocity		all create ${temperature} 93 dist gaussian sum no mom yes rot yes loop all

#Setup output
thermo_style		custom step temp press etotal vol
thermo			100

#Colored thermal bath
fix			scapegoat_qtb all nve									#NVE does the time integration
fix			methane_qtb all qtb temp ${temperature} damp ${damp_qtb} seed 35082 f_max 0.3 N_f 50	#Change f_max if your Debye frequency is higher 
timestep		${delta_t}
run	        	2000											#500 fs
unfix			methane_qtb
unfix			scapegoat_qtb

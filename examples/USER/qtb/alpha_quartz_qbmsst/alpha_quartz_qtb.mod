## This script first constructs an alpha quartz structure of a given size. It then uses fix qtb to equilibrate the computational cell to the specified temperature and pressure.


## This part defines units, alpha-quartz crystal, and atomic information
#General
units			metal
dimension		3
boundary		p p p
atom_style		charge

#Lattice
lattice			custom 1.0 &
			a1	4.916000 0.000000 0.000000 &
			a2	-2.45800 4.257381 0.000000 &
			a3	0.000000 0.000000 5.405400 &
			&		
			basis	0.469700 0.000000 0.000000 &
			basis	0.000000 0.469700 0.666667 &
			basis	0.530300 0.530300 0.333333 &
			&
			basis	0.413500 0.266900 0.119100 &
			basis	0.266900 0.413500 0.547567 &
			basis	0.733100 0.146600 0.785767 &
			basis	0.586500 0.853400 0.214233 &
			basis	0.853400 0.586500 0.452433 &
			basis	0.146600 0.733100 0.880900							#American Mineralogist 65 920 1980 (Space Group 154)

#Computational Cell
region			orthorhombic_unit_cell block 0 4.916000 0 8.514762 0 5.405400 units box 
create_box		2 orthorhombic_unit_cell  
create_atoms		1 box &
			basis	1 1 &
			basis	2 1 &
			basis	3 1 &
			basis	4 2 &
			basis	5 2 &
			basis	6 2 &
			basis	7 2 &
			basis	8 2 &
			basis	9 2 
replicate		${x_rep} ${y_rep} ${z_rep}

#Atomic Information
mass			1 28.085500
mass			2 15.999400
set			type 1 charge +2.4 
set			type 2 charge -1.2


## This part implements the BKS pair potential with a cut-off distance for the Buckingham term. Long range Coulomb interactions are evaluated with the pppm method.
include			alpha_quartz_potential.mod


## This part equilibrates your crystal to a pressure of ${pressure}(unit pressure) and a temperature of ${temperature}(unit temperatureture) with quantum nuclear effects
variable		p_damp equal ${delta_t}*1000								#Recommended pressure damping parameter in fix nph
fix			scapegoat_qtb all nph iso ${pressure} ${pressure} ${p_damp}				#NPH does the time integration
fix			quartz_qtb all qtb temp ${temperature} damp ${damp_qtb} seed 35082 f_max 120.00 N_f 100	#Change f_max (THz) if your Debye frequency is higher
thermo_style		custom step temp press etotal vol lx ly lz pxx pyy pzz pxy pyz pxz
thermo			100
run			2000											# 2 ps
unfix			quartz_qtb
unfix			scapegoat_qtb

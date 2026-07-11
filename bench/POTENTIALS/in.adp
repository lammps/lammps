# bulk Ni in ADP

units		metal
atom_style	atomic

lattice		fcc 3.52
region		box block 0 20 0 20 0 20
create_box	1 box
create_atoms	1 box

pair_style	adp
pair_coeff	* * Ni.adp Ni

velocity	all create 1600.0 376847 loop geom

neighbor	1.0 bin
neigh_modify    delay 5 every 1

fix		1 all nve

timestep	0.005

run		100


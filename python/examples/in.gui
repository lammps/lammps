# 3d Lennard-Jones melt

units		lj
dimension	2
atom_style	atomic

lattice		sq2 0.8442
region		box block 0 30 0 15 -0.5 0.5
create_box	1 box
create_atoms	1 box
mass		1 1.0

velocity	all create 1.44 87287 loop geom

pair_style	lj/cut 2.5
pair_coeff	1 1 1.0 1.0 2.5

neighbor	0.3 bin
neigh_modify	delay 0 every 1 check yes

timestep	0.003

fix		1 all nve
fix		3 all enforce2d

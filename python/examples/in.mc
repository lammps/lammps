# problem setup for Monte Carlo relaxation of perturbed 2d hex lattice
# same as example/MC/in.mc

units		lj
atom_style	atomic
atom_modify     map array sort 0 0.0

dimension       2

lattice		hex 1.0
region		box block 0 10 0 5 -0.5 0.5

create_box	1 box
create_atoms	1 box
mass		1 1.0

pair_style	lj/cut 2.5
pair_coeff	1 1 1.0 1.0 2.5
pair_modify     shift yes

neighbor	0.3 bin
neigh_modify	delay 0 every 1 check yes

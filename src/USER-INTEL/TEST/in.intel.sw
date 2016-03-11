# bulk Si via Stillinger-Weber
package intel 1 mode mixed balance $b
package omp 0
suffix $s

variable	x index 2
variable	y index 2
variable	z index 4

variable	xx equal 20*$x
variable	yy equal 20*$y
variable	zz equal 10*$z

units		metal
atom_style	atomic

lattice		diamond 5.431
region		box block 0 ${xx} 0 ${yy} 0 ${zz}
create_box	1 box
create_atoms	1 box

pair_style	sw
pair_coeff	* * ../../../bench/POTENTIALS/Si.sw Si
mass            1 28.06

velocity	all create 1000.0 376847 loop geom

neighbor	1.0 bin
neigh_modify    delay 5 every 1

fix		1 all nve

timestep	0.001

run             10
run		100

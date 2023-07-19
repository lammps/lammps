# mixture example

units           real
atom_style      full

pair_style      lj/cut/coul/long 12
pair_modify     mix arithmetic

bond_style      harmonic
angle_style     harmonic
dihedral_style  none
improper_style  none

kspace_style    pppm 1e-5

read_data       data.mixture

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# MM dynamics

timestep        0.01

fix		1 all nve

thermo_style    custom step cpu temp ke evdwl ecoul epair emol elong &
                pe etotal press

thermo          10
run		20

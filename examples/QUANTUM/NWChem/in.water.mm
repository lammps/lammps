# MM for water dimer

units		real
atom_style	full

bond_style      harmonic
angle_style     harmonic

read_data       data.water.mm

group           mm molecule 1
group           qm molecule 2

# pair style must define stand-alone short-range Coulombics

pair_style      lj/cut/coul/cut 6.0
pair_coeff      1 1 0.13506 3.166         
pair_coeff      2 2 0.0 1.0         

velocity        all create 300.0 458732

neighbor	1.0 bin
neigh_modify	delay 0 every 1 check yes

fix		1 all nve

timestep        1.0

thermo_style    custom step cpu temp ke evdwl ecoul epair emol elong &
                pe etotal press

thermo          1

run		10

# MM for water dimer

units		real
atom_style	full

bond_style      harmonic
angle_style     harmonic

read_data       data.water.dimer

group           mm molecule 1
group           qm molecule 2

# pair style must define stand-alone short-range Coulombics

pair_style      lj/cut/coul/cut 6.0
pair_coeff      1 1 0.13506 3.166         
pair_coeff      2 2 0.0 1.0         

#velocity        all create 300.0 458732

neighbor	1.0 bin
neigh_modify	delay 0 every 1 check yes

fix		1 all nve

compute         1 all pair/local dist
compute         2 all reduce max c_1

variable        fxabs atom abs(fx)
variable        fyabs atom abs(fy)
variable        fzabs atom abs(fz)
variable        qabs atom abs(q)
compute         3 all reduce max v_fxabs v_fyabs v_fzabs v_qabs

dump            1 all custom 1 dump.water.dimer.mm id x y z q fx fy fz
dump_modify     1 sort id format float "%20.16g"

timestep        1.0

thermo_style    custom step cpu temp ke evdwl ecoul epair emol elong &
                pe etotal press c_2 c_3[*]

thermo          1

run		10

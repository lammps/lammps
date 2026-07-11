units           metal
atom_style      atomic
boundary        p p p

lattice diamond 5.4309
region box block 0 10 0 10 0 10
create_box 1 box
mass 1 28.0855
create_atoms 1 box basis 1 1 basis 2 1 basis 3 1 basis 4 1 basis 5 1 basis 6 1 basis 7 1 basis 8 1

pair_style      sw
pair_coeff      * * Si.sw Si

neighbor        2.0 bin
neigh_modify    every 5 delay 0 check yes

fix             1 all nve
fix             twotemp all ttm/mod 1354684 Si.ttm_mod 10 10 10 set 1000.0 # outfile 100 T_out.txt

compute         pe all pe/atom
compute         ke all ke/atom

timestep        0.0001
thermo          100

thermo_style    custom step temp etotal f_twotemp[1] f_twotemp[2]
                thermo_modify format float "%20.16g"

run             1000

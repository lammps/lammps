
variable newtn index on off
label newtn

variable sfx index nosfx
variable tst index 01

clear

include ${sfx}.inc

units lj
atom_style ${astyle}
special_bonds lj/coul 1.0 1.0 1.0
special_bonds ${sbonds}
newton ${newtn}

read_data ../data/data.${astyle}-${tst}
include ${pstyle}.inc

group left  id 1
group right id 2

velocity all ramp vz -0.1 0.1 z -5.0 5.0 units box

neighbor 0.3 bin

fix  1 all nve

compute 1 left group/group right

variable force atom fx*fx+fy*fy+fz*fz
compute 2 all pe/atom
compute 3 all reduce sum c_2
compute 4 left reduce sum v_force
compute 5 right reduce sum v_force

thermo 1
thermo_style custom step etotal pe c_3 ke evdwl ecoul epair etail &
             pxx pyy pzz pxy pxz pyz &
             c_1[1] c_1[2] c_4 c_5 fmax fnorm
thermo_modify format float %-16.10f
run 10 post no

next newtn
jump SELF newtn


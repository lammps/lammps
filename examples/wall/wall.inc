
clear
read_restart walltest.restart
# log ${name}.log

compute ea move pe/atom
compute pe move reduce sum c_ea
special_bonds lj/coul 1.0 1.0 1.0

thermo_style custom step temp press ke pe etotal c_pe
thermo 1
thermo_modify format float %14.10g

fix 0 move nve
# dump 1 all custom 1 ${name}.dump id x fx
# dump_modify 1 format float %20.15g

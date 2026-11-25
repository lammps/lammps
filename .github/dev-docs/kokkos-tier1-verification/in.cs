include         mol_head.mod
variable        ps index coul/long/cs
if "${ps} == coul/long/cs" then "pair_style ${ps} 9.0" &
elif "${ps} == coul/wolf/cs" "pair_style ${ps} 0.2 9.0" &
elif "${ps} == born/coul/wolf/cs" "pair_style ${ps} 0.2 9.0 9.0" &
elif "${ps} == born/coul/dsf/cs" "pair_style ${ps} 0.1 9.0 9.0" &
elif "${ps} == born/coul/dsf" "pair_style ${ps} 0.1 9.0 9.0" &
elif "${ps} == born/coul/long/cs" "pair_style ${ps} 9.0 9.0" &
elif "${ps} == buck/coul/long/cs" "pair_style ${ps} 9.0 9.0" &
elif "${ps} == lj/cut/coul/long/cs" "pair_style ${ps} 9.0 9.0" &
else "pair_style ${ps} 9.0 9.0"
if "${ps} == coul/long/cs" then "pair_coeff * *" &
elif "${ps} == coul/wolf/cs" "pair_coeff * *" &
elif "${ps} == born/coul/wolf/cs" "pair_coeff * * 1.0 0.3 0.0 1.0 0.5" &
elif "${ps} == born/coul/dsf/cs" "pair_coeff * * 1.0 0.3 0.0 1.0 0.5" &
elif "${ps} == born/coul/dsf" "pair_coeff * * 1.0 0.3 0.0 1.0 0.5" &
elif "${ps} == born/coul/long/cs" "pair_coeff * * 1.0 0.3 0.0 1.0 0.5" &
elif "${ps} == buck/coul/long/cs" "pair_coeff * * 100.0 0.3 10.0" &
elif "${ps} == lj/class2/coul/long/cs" "pair_coeff * * 0.02 3.0" &
else "pair_coeff * * 0.02 3.0"
if "${ps} == coul/long/cs" then "kspace_style ewald 1.0e-6" &
elif "${ps} == born/coul/long/cs" "kspace_style ewald 1.0e-6" &
elif "${ps} == buck/coul/long/cs" "kspace_style ewald 1.0e-6" &
elif "${ps} == lj/cut/coul/long/cs" "kspace_style ewald 1.0e-6" &
elif "${ps} == lj/class2/coul/long/cs" "kspace_style ewald 1.0e-6"
fix             1 all nve
thermo_style    custom step temp pe evdwl ecoul elong press
thermo          5
run             20 post no

# sample LAMMPS input script for thermal conductivity of liquid LJ
# Muller-Plathe method via fix thermal_conductivity

# settings

variable        x equal 10
variable        y equal 10
variable        z equal 20

variable        rho equal 0.6
variable        t equal 1.35
variable        rc equal 2.5

#variable       rho equal 0.85
#variable        t equal 0.7
#variable       rc equal 3.0

# setup problem

units           lj
atom_style      atomic

lattice         fcc ${rho}
region          box block 0 $x 0 $y 0 $z
create_box      1 box
create_atoms    1 box
mass            1 1.0

velocity        all create $t 87287

pair_style      lj/cut ${rc}
pair_coeff      1 1 1.0 1.0

neighbor        0.3 bin
neigh_modify    delay 0 every 1

# 1st equilibration run

fix             1 all nvt temp $t $t 0.5
thermo          100
run             1000

velocity        all scale $t

unfix           1

# 2nd equilibration run

compute         ke all ke/atom
variable        temp atom c_ke/1.5

fix             1 all nve

compute         layers all chunk/atom bin/1d z lower 0.05 units reduced
fix             2 all ave/chunk 10 100 1000 layers v_temp file profile.mp
fix             3 all thermal/conductivity 10 z 20

variable        tdiff equal f_2[11][3]-f_2[1][3]
thermo_style    custom step temp epair etotal f_3 v_tdiff
thermo_modify   colname f_3 E_delta colname v_tdiff dTemp_step

thermo          1000
run             20000

# thermal conductivity calculation
# reset fix thermal/conductivity to zero energy accumulation
fix             3 all thermal/conductivity 10 z 20

variable start_time equal time
variable kappa equal (f_3/(time-${start_time})/(lx*ly)/2.0)*(lz/2.0)/f_ave

fix             ave all ave/time 1 1 1000 v_tdiff ave running
thermo_style    custom step temp epair etotal f_3 v_tdiff f_ave
thermo_modify   colname f_3 E_delta colname v_tdiff dTemp_step colname f_ave dTemp

run             20000
print           "Running average thermal conductivity: $(v_kappa:%.2f)"

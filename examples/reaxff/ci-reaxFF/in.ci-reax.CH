#ci-reax potential for CH systems with tabulated ZBL correction
atom_style      charge
units           real

read_data       CH4.dat

pair_style      hybrid/overlay reaxff control checkqeq no table linear 11000
pair_coeff      * * reaxff ffield.ci-reax.CH C H
pair_coeff      1 1 table ci-reaxFF_ZBL.dat CC_cireaxFF 
pair_coeff      1 2 table ci-reaxFF_ZBL.dat CH_cireaxFF 
pair_coeff      2 2 table ci-reaxFF_ZBL.dat HH_cireaxFF

timestep        0.25
fix             1 all nve
fix             2 all temp/berendsen 500.0 500.0 100.0

#dump           1 all atom 30 dump.ci-reax.lammpstrj

run             3000

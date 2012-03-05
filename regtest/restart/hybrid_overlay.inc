pair_style hybrid/overlay lj/cut 2.5 coul/long 5.0
pair_coeff 1 1 lj/cut 1.0 1.0
pair_coeff 1 1 coul/long 
kspace_style ewald 0.001
variable coeff delete
variable coeff index ""

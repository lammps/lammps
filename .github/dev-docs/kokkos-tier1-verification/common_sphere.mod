units           lj
atom_style      sphere
atom_modify     map array
lattice         fcc 0.8442
region          box block 0 8 0 8 0 8
create_box      1 box
create_atoms    1 box
variable        dia index 1.0
variable        den index 1.0
set             group all diameter ${dia} density ${den}
velocity        all create 1.2 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes

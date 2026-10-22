# shared LJ melt setup used by the fix/compute/region decks
units           lj
atom_style      atomic
atom_modify     map array
variable        bnd index "p p p"
boundary        ${bnd}
lattice         fcc 0.8442
region          box block 0 8 0 8 0 8
create_box      2 box
create_atoms    1 box
mass            * 1.0
set             group all type/fraction 2 0.35 5738
velocity        all create 1.2 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes

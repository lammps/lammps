variable        tdir index ../../../unittest/force-styles/tests
units           real
atom_style      full
atom_modify     map array
pair_style      zero 10.0
bond_style      zero
angle_style     zero
dihedral_style  zero
improper_style  zero
read_data       ${tdir}/data.fourmol
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

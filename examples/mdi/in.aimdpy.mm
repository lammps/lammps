# 3d Lennard-Jones melt - LAMMPS as MM engine for aimd_driver.py

variable        x index 5
variable        y index 5
variable        z index 5

units           lj
atom_style      atomic

lattice         fcc 0.8442
region          box block 0 $x 0 $y 0 $z
create_box      1 box
create_atoms    1 box
mass            1 1.0

velocity        all create 1.44 87287 loop geom

neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes

fix             1 all nve

mdi             engine

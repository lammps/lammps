#!/bin/sh

# Aysun Itai's LAMMPS files contain two molecules:

# The CNAD molecule has molecule-id 1

ltemplify.py -name CNAD -molid "1" cnad-cnt.in cnad-cnt.data > cnad.lt

# The CNT (carbon nanotube) corresponds to molecule-id 2
ltemplify.py -name CNT -molid "2" cnad-cnt.in cnad-cnt.data > cnt.lt

# This will extract both molecules and save them as separate .LT files.
# (We can include these files later when preparing new simulations.)

This directory contains scripts used for converting AMBER parameter files
into moltemplate (.LT) format.  When a newer version of the AMBER parameters
is eventually published, you can use these scripts to convert the new files
again.  (Some tinkering may be necessary.)

The main bash script is a wrapper which simply splits up the parameter (".dat")
file into fragments which (it thinks) correspond to the mass, pair, bond,
angle, dihedral, and improper section of the original .dat file.
(However sometimes it gets this wrong and you have to split it up manually!)

Then this bash script invokes the relevant python script to convert 
each section into .LT format:
amberparm_to_mass.py
amberparm_to_pair.py
amberparm_to_bond.py
amberparm_to_angle.py
amberparm_to_dihedral.py
amberparm_to_improper.py
In case this goes wrong, you may have to run these scripts manaully.


Find out how to run this bash script by invoking it without any arguments:

./amberparm2lt.sh

-Andrew
2013-8-04

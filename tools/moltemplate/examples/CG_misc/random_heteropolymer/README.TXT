This is a very crude attempt to create a polymer melt
consisting of 150 copies of the same polymer.

THE FORCE FIELD PARAMETERS AND CONFORMATION FOR THIS SYSTEM ARE COMPLETELY WRONG

The purpose of this example is to demonstrate one way to create a 
large number of randomly generated polymers, and to use 
an NPT simulation to pack them all into a small box.
(Smaller than their initial outstretched length.)

Each polymer is a random heteropolymer of length 200 monomers.
This polymer is a polymoer of PVDF and PCTFE monomers 
selected randomly in a 1:3 ratio.
Again, this is not a realistic simulation of PVDV or PCTFE polymers,

Furthermore, even after the simulation is done, the arrangement of the 
polymers in the box is not characteristic of a truly random polymer melt.


---- A note on size ---

This is a large system with nearly 200000 atoms.
It takes several minutes to compile this example
and (currently requiring at least 4.0 Gb of memory).
(I'm working on reducing that requirement. -Andrew 2012-9-12)

-----------------------

Instructions on how to build LAMMPS input files and 
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh

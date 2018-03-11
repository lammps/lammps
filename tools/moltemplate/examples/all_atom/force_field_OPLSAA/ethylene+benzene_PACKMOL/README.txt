
This is an example of how to use the OPLSAA force-field in LAMMPS

This example also shows how to use moltemplate in combination with PACKMOL.
(PACKMOL is a useful program for generating atomic coordinates. In this example,
 moltemplate.sh is only used to create the topology, force-field and charges,
 and PACKMOL generates the coordinates, which moltemplate reads (in "step 1").
 Moltemplate can also be used for generating atomic coordinates, especially
 for mixing many small molecules together, as we do in this example.  However
 I wanted to demonstrate how to combine PACKMOL with moltemplate.sh.
 In some other scenarios, such as protein solvation, PACKMOL does a much
 better job than moltemplate.)

As of 2016-11-21, this code has not been tested for accuracy.
(See the WARNING.TXT file.)

step 1)
To build the files which LAMMPS needs, follow the instructions in:
README_setup.sh

step 2)
To run LAMMPS with these files, follow these instructions:
README_run.sh

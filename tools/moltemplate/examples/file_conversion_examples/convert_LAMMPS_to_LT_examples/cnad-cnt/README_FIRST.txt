###########################################################
# Interaction of a carbon nanotube with a pair of "mystery
# molecules" (extracted from the cnat-cnt.data/in files).
###########################################################
# Author: Aysun Itai and Andrew Jewett

This example uses "ltemplify.py" to create molecule templates out
of two different molecules in a pre-existing LAMMPS IN/DATA file.
Then I show how to use "moltemplate.sh" to make copies of these
molecules and to move and rotate them (creating new LAMMPS IN/DATA files).

   Disclaimer:
The molecules in this example are not physically realistic.
The purpose of this example is to demonstrate ltemplify usage.

   REQUIRED INPUT FILES

cnad-cnt.data cnad-cnt.in system.lt

 cnad-cnt.data
 This is a LAMMPS data file containing the coordinates and the topology
 for a system combining the two molecules together.  ltemplify will extract
 molecules from this file, one at a time.

 cnad-cnt.in
 This file contains force-field parameters and old run settings for the system.
 (We ignore the run settings in this file.)  The force-field parameters in
 the "cnad-cnt.in" file are only necessary because we are going to build
 a completely new set of simulation input files. (We are not only going to
 rotate them and duplicate the molecules.)  ltemplify.py will extract the
 force field parameters from this file.  This approach allows us to combine
 these molecules with other types of molecules later on.)

 system.lt
 The "system.lt" contains the instructions what we will do with these molecules
 after ltemplify.py has converted them into .LT format.  In this example
 it contains instructions for rotating and copying the two molecules,
 (It also defines the periodic boundary conditions.)

   OUTPUT FILES

cnad.lt
cnt.lt

These files are referenced in system.lt.
Running moltemplate.sh on system.lt  (using "moltemplate.sh system.lt")
creates new LAMMPS data and input files:
system.data, system.in, system.in.init, system.in.settings
(These files are referenced in run.in.nvt.)

You can run a simulation from the files created by moltemplate using

lmp_linux -i run.in.nvt

NOTE: BECAUSE ALL OF THE ORIGINAL FORCE FIELD PARAMETERS WERE INTENTIONALLY
      ALTERED, THE SYSTEM WILL MOVE IN A VERY UNREALISTIC WAY WHEN SIMULATED.
      (This was done to protect the original source of the files.)
      The goal of this example is only to demonstrate how to use
      "ltemplify.py" to convert lammps input and data files into
      LT format and back again.)

    -----------

Instructions:
Run the commands (follow the instructions) in these files:

step 1)
README_step1_run_ltemplify.sh

and then

step 2)
README_step2_run_moltemplate.sh

step 3)  OPTIONAL

To run a short LAMMPS simulation, you can use the "in.nvt" file, for example:

$LAMMPS_BINARY -i run.in.nvt

where "$LAMMPS_BINARY" is the name of the command you use to invoke lammps
(such as lmp_linux, lmp_g++, lmp_mac, lmp_ubuntu, lmp_cygwin, etc...).
    -----------

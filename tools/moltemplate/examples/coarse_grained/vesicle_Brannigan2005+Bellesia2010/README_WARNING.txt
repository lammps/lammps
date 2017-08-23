WARNING:

   This is not a realistic simulation of proteins in a lipid membrane.  This
example was only intented to be a technical demonstration to show how to
combine totally different kinds of coarse-grained molecules (with different
kinds of force-fields) together in the same simulation in LAMMPS.  Tuning the
force-field parameters to get realistic results was not the goal.  I did
not take the extra time to do this.  If you have suggestions for changes,
please email me (jewett.aij at gmail dot com).

   In addition, I have noticed that newer versions of PACKMOL do not
always succeed at generating a spherical vesicle in a reasonable amount of time.
(You may have to play with the .inp files in the packmol_files directory
 to get PACKMOL to produce any files at all.

(NOTE: This example also demonstrantes how to use an external program
 ("packmol") to generate the coordinates for the atoms in the system.
 PLEASE USE "packmol", NOT "ppackmol". -the parallel version of "packmol".
 This is because "ppackmol" is more likely to get caught in infinite loops.)


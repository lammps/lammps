This example using the electron Force Field (eFF) was created by
Andres Jaramillo-Botero and distributed with LAMMPS in the
"examples/USERS/eff/CH4" subdirectory.
The files from that example were converted into moltemplate format using
"ltemplify.py" and then edited by hand (to rename the atom types,
and replace all of the "pair_coeff ..." commands with "pair_coeff * *")


--- Original README: ---
Methane, valence electron ionization and full molecule tests (spe, dynamics).
Note: electron mass set to 1


-----
WARNING: Regarding the "run.in.ch4_ionized" file
         As of 2014-3-12, the "pair_style eff/cut 5000.0 0 0" command
         located in "orig_files/in.ch4_ionized.dynamics" (as well as the
         files "moltemplate_files/ch4_ionized.lt" and "system.in.settings",
         which are both derived from it) causes LAMMPS to hang.
         Running LAMMPS on Andres' original eFF example has the same behavior.
         This appears to be an eFF/LAMMPS issue (not a moltemplate issue).
         The "pair_style eff/cut 100" command works, so
         try reducing the cutoff (or ask Andres Jaramillo-Botero for help).
         Please let me know if you solve this issue (jewett.aij -at- g mail)


How to add a fix to LAMMPS:

1. Put the .h and .cpp files in the src/ folder, and make sure it follows the same style as other fix files already included in LAMMPS.

2. Edit style_fix.h and add an your fix_foo.h file to the list of other fixes.

3. Put fix_sicrack files and the other crack script files in separate folder "USER-SICRACK", and make an "Install.csh" script which will (re)move files to LAMMPS src directory when the command "make yes-user-sicrack" or "make no-user-sicrack" is issued before making LAMMPS

4. Edit "Makefile" in the src-directory, add "user-sicrack" in the list of package variables among the other alternative packages, in the "PACKUSER =" command

5. Make LAMMPS again.

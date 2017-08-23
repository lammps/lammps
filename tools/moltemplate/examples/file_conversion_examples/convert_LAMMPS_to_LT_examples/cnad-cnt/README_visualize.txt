
 ------- To view the trajectory in VMD --------


1) Build a PSF file for use in viewing with VMD.

This step works with VMD 1.9 and topotools 1.2.
(Older versions, like VMD 1.8.6, don't support this.)


a) Start VMD
b) Menu  Extensions->Tk Console
c) Enter:

(I assume that the the DATA file is called "system.data")

   topo readlammpsdata system.data full
   animate write psf system.psf


Later, to Load a trajectory in VMD:
  Start VMD
  Select menu: File->New Molecule
 -Browse to select the PSF file you created above, and load it.
  (Don't close the window yet.)
 -Browse to select the trajectory file.
  If necessary, for "file type" select: "LAMMPS Trajectory"
  Load it

-----  Wrap the coordinates to the unit cell

a) Start VMD
b) Load the trajectory in VMD (see above)
c) Menu  Extensions->Tk Console
d) Enter:

    DISCLAIMER: I'M NOT SURE THESE COMMANDS ARE CORRECT.
                LOOKUP "pbctools" FOR DETAILS.

    pbc wrap -compound res -all
    pbc box

3) Optional: If you like, change the atom types in the PSF file so
   that VMD recognizes the atom types, use something like:

sed -e 's/   1    1      /   C    C      /g' < system.psf > temp1.psf
sed -e 's/   2    2      /   H    H      /g' < temp1.psf  > temp2.psf
sed -e 's/   3    3      /   P    P      /g' < temp2.psf  > system.psf

(If you do this, it might effect step 2 above.)

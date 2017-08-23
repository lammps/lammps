 ------- Instructions to view a trajectory in VMD --------

 ------- Disclaimer -------

1) Build a PSF file for use in viewing with VMD.

This step works with VMD 1.9 and topotools 1.2.
(Older versions, like VMD 1.8.6, don't support this.)


a) Start VMD
b) Menu  Extensions->Tk Console
c) Enter:

(I assume that the the DATA file is called "system.data")

   topo readlammpsdata system.data full
   animate write psf system.psf

   (Note, at this point the image shown in the VMD graphics window may
   not appear correct or incomplete.  The coordinates of the atoms may
   overlap if you asked moltemplate.sh to load your coordinates from
   a PDB or XYZ file.
   However, later after you have run a simulation, the trajectories
   should appear reasonably correct when you load them in VMD using
   the PSF file you just generated.)


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

    # If you have a solute of type 1, then use this:
    #pbc wrap -sel type=1 -all -centersel type=2 -center com

"1" corresponds to the "O" atom type
"2" corresponds to the "H" atom type

3) Optional: If you like, change the atom types in the PSF file so
   that VMD recognizes the atom types:

sed -e 's/   1    1      /   O    O      /g' < system.psf > temp1.psf
sed -e 's/   2    2      /   H    H      /g' < temp1.psf > system.psf

(If you do this, I guess that you might have to use
 "type=O" and "type=H" in step 2 above.)

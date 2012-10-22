
 ------- To view a lammps trajectory in VMD --------


1) Build a PSF file for use in viewing with VMD.

This step works with VMD 1.9 and topotools 1.2.  
(Older versions, like VMD 1.8.6, don't support this.)


a) Start VMD
b) Menu  Extensions->Tk Console
c) Enter:

(I assume that the the DATA file is called "system.data")

   topo readlammpsdata system.data full
   animate write psf system.psf

2) 

Later, to Load a trajectory in VMD:

  Start VMD
  Select menu: File->New Molecule
 -Browse to select the PSF file you created above, and load it.
  (Don't close the window yet.)
 -Browse to select the trajectory file.
  If necessary, for "file type" select: "LAMMPS Trajectory"
  Load it.

   ----  A note on trajectory format: -----
If the trajectory is a DUMP file, then make sure the it contains the
information you need for pbctools (see below.  I've been using this 
command in my LAMMPS scripts to create the trajectories:

  dump 1 all custom 5000 DUMP_FILE.lammpstrj id mol type x y z ix iy iz

It's a good idea to use an atom_style which supports molecule-ID numbers 
so that you can assign a molecule-ID number to each atom.  (I think this 
is needed to wrap atom coordinates without breaking molecules in half.)

Of course, you don't have to save your trajectories in DUMP format, 
(other formats like DCD work fine)  I just mention dump files 
because these are the files I'm familiar with.

3) -----  Wrap the coordinates to the unit cell
          (without cutting the molecules in half)

a) Start VMD
b) Load the trajectory in VMD (see above)
c) Menu  Extensions->Tk Console
d) Try entering these commands:

    pbc wrap -compound res -all
    pbc box

    ----- Optional ----
    It can help to shift the location of the periodic boundary box 
    To shift the box in the y direction (for example) do this:

    pbc wrap -compound res -all -shiftcenterrel {0.0 0.15 0.0}
    pbc box -shiftcenterrel {0.0 0.15 0.0}

    Distances are measured in units of box-length fractions, not Angstroms.

    Alternately if you have a solute whose atoms are all of type 1, 
    then you can also try this to center the box around it:

    pbc wrap -sel type=1 -all -centersel type=2 -center com

4) 
    You should check if your periodic boundary conditions are too small.
    To do that:
       select Graphics->Representations menu option
       click on the "Periodic" tab, and 
       click on the "+x", "-x", "+y", "-y", "+z", "-z" checkboxes.

5) Optional: If you like, change the atom types in the PSF file so 
   that VMD recognizes the atom types, use something like:

sed -e 's/   1    1      /   C    C      /g' < system.psf > temp1.psf
sed -e 's/   2    2      /   H    H      /g' < temp1.psf  > temp2.psf
sed -e 's/   3    3      /   P    P      /g' < temp2.psf  > system.psf

(If you do this, it might effect step 2 above.)

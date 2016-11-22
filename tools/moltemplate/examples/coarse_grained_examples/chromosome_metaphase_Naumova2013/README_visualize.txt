NOTE: VMD DOES NOT ALLOW YOU TO VISUALIZE SYSTEMS WITH MANY BONDS ATTACHED
      TO EACH ATOM.  (IF IT DID, THE RESULTS WOULD BE UGLY ANWAY.)

HOWEVER THIS MODEL ATTACHES APPROXIMATELY 60 BONDS TO EACH CONDENSIN ATOM.
IN ORDER TO PULL THE CONDENSIN MONOMERS TOGETHER.  YOU MUST DELETE THOSE
BONDS (of type "1" or "2") FROM THE "system.data" FILE BEFORE YOU CARRY
OUT THE COMMANDS BELOW.  (...And backup your "system.data" file.  You'll need 
all the bonds when you run the simulations.)

-------------- COLORS ---------------
In order to show how the polymer is distributed along the length of the
cylinder, I recommend to select the
Graphics->Graphical Representations 
menu option, and select "Index" from the "Coloring Method" pull-down menu.

After doing this, you can switch from a red-white-blue scheme, to a 
rainbow ("jet") scheme, by selecting the Extensions->Tk Console menu option
and loading the "vmd_colorscale_jet.tcl" file located in the "images" directory.
-------------------------------------------

First, if you have not done so, download and install VMD:

http://www.ks.uiuc.edu/Research/vmd/
http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD


 ------- To view a lammps trajectory in VMD --------

The system coordinates are initialy stored in a LAMMPS' ".data" file.
(If that file was built with moltemplate, it will be named "system.data".)

The first step is to view that file.
Then you should create a ".psf" file
(The .psf file is necessary after you run the simulation 
 for viewing trajectories.)

1) Build a PSF file for use in viewing with VMD 

a) Start VMD
b) Menu  Extensions->Tk Console
c) Enter:

(I assume that the the DATA file is called "system.data")

   topo readlammpsdata system.data full
   animate write psf system.psf

You will see a snapshot of the system on the screen.
(presumably the initial conformation at t=0)

2) 

Later once you have run a simulation, 
to Load a trajectory in VMD:

  Start VMD
  Select menu: File->New Molecule
 -Browse to select the PSF file you created above, and load it.
  (Don't close the window yet.)
 -Browse to select the trajectory file
  (It usually has names like "traj.lammpstrj".  It depends on how you saved it.)
  If necessary, for "file type" select: "LAMMPS Trajectory".
  (However VMD should recognize the file type by the file extension.)
  Load it.





##################### PERIODIC BOUNDARY CONDITIONS #####################
  If you are only simulating a single molecule and you are not 
  using periodic boundary conditions, then ignore everything below.
########################################################################

   ----  A note on trajectory format: -----
If the trajectory is the standard LAMMPS format, (aka a "DUMP" file with 
a ".lammpstrj" extension), then it's a good idea when you run the simulation
to tell LAMMPS you want to store the information needed for showing periodic
boundary conditions.  (Even if you are not using periodic boundaries.  
It never hurts to include a tiny bit of extra information.) To do that,
I've been using this command in my LAMMPS scripts to create the trajectories:

  dump 1 all custom 5000 traj.lammpstrj id mol type x y z ix iy iz

(Also: it's a good idea to use an atom_style which supports molecule-ID numbers 
so that you can assign a molecule-ID number to each atom. I think this is needed
to wrap atom coordinates visually without breaking molecules in half. Again
you don't need to worry about this if you are not using periodic boundaries.)


3) -----  Wrap the coordinates to the unit cell
          (without cutting the molecules in half)

a) Start VMD
b) Load the trajectory in VMD (see above)
c) Menu  Extensions->Tk Console
d) Try entering these commands:

    pbc wrap -compound res -all
    pbc box

    ----- Optional ----
    Sometimes the solvent or membrane obscures the view of the solute.
    It can help to shift the location of the periodic boundary box 
    To shift the box in the y direction (for example) do this:

    pbc wrap -compound res -all -shiftcenterrel {-0.5 -0.5 -0.5}
    pbc box -shiftcenterrel {-0.5 -0.5 -0.5}

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

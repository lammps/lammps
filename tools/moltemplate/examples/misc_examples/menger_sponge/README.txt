NOTE: This example requires the "Al99.eam.alloy" file.
      (It was not included in this directory because if its large size.)
      As of 2012-11, I was able to obtain it here:
      http://www.ctcms.nist.gov/~cbecker/Download/Al-YM/Al99.eam.alloy
      Copy it to the directory containing this README file.
------------------------------------------------------------------------
3D fractal test

Moltemplate is useful for building larger molecular structures from smaller
pieces.  The purpose of this example is to show how to build large (many-level)
heirarchical objects (a Menger sponge) using moltemplate.

A Menger sponge is a fractal composed of subunits that resemble a Rubik's-cubes
with a hollow interior:

      ___
     ###|
     # #|
     ###'

      |
     \|/
      V

      _________
     /        /|
     ######### |
     # ## ## # |
     ######### |
     ###   ### |
     # #   # # |
     ###   ### |
     ######### |
     # ## ## #/
     #########

In this example, we will build a periodic lattice of Menger sponges.

The smallest cube subunits consist of 4 atoms of Aluminum
(arranged in a cubic FCC unit-cell for bulk Aluminum).
This was an arbitrary choice.  The resulting construct is not stable.

  --- To build the system ---

Carry out the instructions in README_setup.sh,
to generate the LAMMPS DATA file and input scripts you need:
system.data, system.in.init, system.in.settings.
(The run.in script contains references to these files.)

  --- To run LAMMPS, try a command like: ---

lmp_mpi -i run.in

    or (if you have mpi installed)

mpirun -np 4 lmp_mpi -i run.in

This will create an ordinary LAMMPS dump file you can visualize with VMD
traj.lammpstrj    (See README_visualize.txt)

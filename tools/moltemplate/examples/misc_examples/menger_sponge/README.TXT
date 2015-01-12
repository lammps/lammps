NOTE: This example requires the "Al99.eam.alloy" file.
      (It was not included in this directory because if its large size.)
      As of 2012-11, I was able to obtain it here:
      http://www.ctcms.nist.gov/~cbecker/Download/Al-YM/Al99.eam.alloy
      Copy it to the directory containing this README file.
------------------------------------------------------------------------
3D fractal test

Moltemplate is useful for building larger molecular structures from 
smaller pieces.  Although this simulation is of no scientific value, thiss 
example illustrates how to build large (many-level) heirarchical objects 
(Serpinski cubes) using moltemplate.  (This is also called a "Menger Sponge".)

The files in this directory demonstrate a way to build a periodic lattice of 
3-dimensional Serpinski-cubes (with 3 levels of recursive self-similarity).

In this example, the basic indivisible units are 4-atoms of Aluminum 
(arranged in a cubic FCC unit-cell for bulk Aluminum).
This was an arbitrary choice.  The resulting construct is not stable.
(But it makes pretty movies while collapsing.)

To understand what is going on with this example, look at this file:

./moltemplate_files/elegant_inefficient_version/serpinski_cubes.lt

(This approach uses too much memory to be practical for large simulaions.
The version I actually use is here: ./moltemplate_files/serpinski_cubes.lt)

  --- To build the system ---

Carry out the instructions in README_setup.sh,
to generate the LAMMPS DATA file and input scripts you need:
system.data, system.in.init, system.in.settings.
(The run.in script contains references to these files.)

  --- To run LAMMPS, try a command like: ---

lmp_linux -i run.in

    or (if you have mpi installed)

mpirun -np 4 lmp_linux -i run.in

This will create an ordinary LAMMPS dump file you can visualize with VMD
traj.lammpstrj    (See README_visualize.txt)

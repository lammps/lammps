This directory  used    to contain utility  scripts  for   using  VMD to
visualize and analyze LAMMPS  trajectories. As of  April 2010 all of the
scripts and many additional features have been merged into the topotools
plugin that is bundled with VMD. Updates between VMD releases are here:
http://sites.google.com/site/akohlmey/software/topotools
This page also contains detailed documentation and some tutorials.

These scripts within VMD and the plugin for native LAMMPS dump files are
are  maintained  by Axel Kohlmeyer <akohlmey@gmail.com>;  please contact
him through the LAMMPS mailing list in case of problems.

Below are a few comments on support for LAMMPS in VMD.

-------------------------

1. File formats and VMD limitations

   VMD currently supports reading several but not all file formats
   that LAMMPS can generate. Supported are: atom (text mode), custom
   (text mode, only some fields are directly supported, please see
   below for more details), dcd, xyz and xtc.  Cfg and binary native
   dump files are not supported (06/2012). The new molfile dump style
   in addition allows to use VMD molfile plugins to write dumps in
   any format that is supported by VMD.

   However VMD requires all frames of a file to have the same number of
   atoms.  If the number of atoms changes between two frames,  the file
   reader will stop. The topotools plugin has a special scripted file
   reader for .xyz files that can generate the necessary padding so that
   the file can still be read into VMD.  Whether an atom is real or 
   "invisible" is then flagged in the "user" field.  For efficiency
   reasons this script will not preserve atom identity between frames.

2. Topology files, a.k.a. as "data" files

   The topotools plugin also contains a read and write option for LAMMPS
   data files. This reader will try to preserve as much information as 
   possible and will also store useful information as comments upon
   writing. It does not store or read coefficient data.  In combination
   with other functionality in topotools complete topologies for rather
   complicated systems for LAMMPS can be build with VMD scripting.

3. Reading custom data fields into VMD

   At this moment VMD only supports reading coordinates and velocities
   (if present) as per timestep data. Everthing else is just taken
   from the first frame or whatever file was used to generate this
   structure information.  Through setting the environment variable
   LAMMPSREMAPFIELDS, custom properties can be mapped to the x, y, z,
   vx, vy, vz data fields and then accessed from within VMD. For
   example to store radius and charge of a particle in the vx and vy
   fields, respectively set this variable to "vx=radius,vy=q".  Future
   versions of VMD will allow more flexibility.

4. Recovering information about elements

   Colorization in VMD is by default based on atom names, yet LAMMPS
   requires identification of atoms by consecutive numbers starting at
   1. With the option of reading a LAMMPS data file, additional
   information is provided that can help to recover some of this
   data. 'topo guessatom element mass' will guess the atom's element
   name from it mass (with fuzz, where possible).

5. Reading files from the command line

   Converting a LAMMPS data file to a .psf file can be very convenient
   for loading trajectories from the command line. This conversion is
   done with

   topo readlammpsdata data.rhodo full
   animate write psf rhodo.psf

   In the future you can now load this PSF file first and then the
   LAMMPS dump file(s) (or a more compact and faster loading DCD or
   XTC file) with:

     vmd micelle.psf -lammpstrj dump.micelle

   Note how the -lammpstrj flag will tell VMD that dump.micelle is a
   lammps trajectory file.

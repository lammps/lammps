.. index:: dump

dump command
============

:doc:`dump vtk <dump_vtk>` command
==================================

:doc:`dump h5md <dump_h5md>` command
====================================

:doc:`dump molfile <dump_molfile>` command
==========================================

:doc:`dump netcdf <dump_netcdf>` command
========================================

:doc:`dump image <dump_image>` command
======================================

:doc:`dump movie <dump_image>` command
======================================

:doc:`dump atom/adios <dump_adios>` command
===========================================

:doc:`dump custom/adios <dump_adios>` command
=============================================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID style N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* style = *atom* or *atom/gz* or *atom/mpiio* or *cfg* or *cfg/gz* or
  *cfg/mpiio* or *custom* or *custom/gz* or *custom/mpiio* or *dcd* or *h5md* or *image* or *local* or *local/gz* or *molfile* or *movie* or *netcdf* or *netcdf/mpiio* or *vtk* or *xtc* or *xyz* or *xyz/gz* or *xyz/mpiio*
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = list of arguments for a particular style

  .. parsed-literal::

       *atom* args = none
       *atom/gz* args = none
       *atom/mpiio* args = none
       *atom/adios* args = none,  discussed on :doc:`dump atom/adios <dump_adios>` doc page
       *cfg* args = same as *custom* args, see below
       *cfg/gz* args = same as *custom* args, see below
       *cfg/mpiio* args = same as *custom* args, see below
       *custom*\ , *custom/gz*\ , *custom/mpiio* args = see below
       *custom/adios* args = same as *custom* args, discussed on :doc:`dump custom/adios <dump_adios>` doc page
       *dcd* args = none
       *h5md* args = discussed on :doc:`dump h5md <dump_h5md>` doc page
       *image* args = discussed on :doc:`dump image <dump_image>` doc page
       *local* args = see below
       *molfile* args = discussed on :doc:`dump molfile <dump_molfile>` doc page
       *movie* args = discussed on :doc:`dump image <dump_image>` doc page
       *netcdf* args = discussed on :doc:`dump netcdf <dump_netcdf>` doc page
       *netcdf/mpiio* args = discussed on :doc:`dump netcdf <dump_netcdf>` doc page
       *vtk* args = same as *custom* args, see below, also :doc:`dump vtk <dump_vtk>` doc page
       *xtc* args = none
       *xyz* args = none
       *xyz/gz* args = none
       *xyz/mpiio* args = none

* *custom* or *custom/gz* or *custom/mpiio* or *netcdf* or *netcdf/mpiio* args = list of atom attributes

  .. parsed-literal::

         possible attributes = id, mol, proc, procp1, type, element, mass,
                               x, y, z, xs, ys, zs, xu, yu, zu,
                               xsu, ysu, zsu, ix, iy, iz,
                               vx, vy, vz, fx, fy, fz,
                               q, mux, muy, muz, mu,
                               radius, diameter, omegax, omegay, omegaz,
                               angmomx, angmomy, angmomz, tqx, tqy, tqz,
                               c_ID, c_ID[N], f_ID, f_ID[N], v_name

  .. parsed-literal::

           id = atom ID
           mol = molecule ID
           proc = ID of processor that owns atom
           procp1 = ID+1 of processor that owns atom
           type = atom type
           element = name of atom element, as defined by :doc:`dump_modify <dump_modify>` command
           mass = atom mass
           x,y,z = unscaled atom coordinates
           xs,ys,zs = scaled atom coordinates
           xu,yu,zu = unwrapped atom coordinates
           xsu,ysu,zsu = scaled unwrapped atom coordinates
           ix,iy,iz = box image that the atom is in
           vx,vy,vz = atom velocities
           fx,fy,fz = forces on atoms
           q = atom charge
           mux,muy,muz = orientation of dipole moment of atom
           mu = magnitude of dipole moment of atom
           radius,diameter = radius,diameter of spherical particle
           omegax,omegay,omegaz = angular velocity of spherical particle
           angmomx,angmomy,angmomz = angular momentum of aspherical particle
           tqx,tqy,tqz = torque on finite-size particles
           c_ID = per-atom vector calculated by a compute with ID
           c_ID[I] = Ith column of per-atom array calculated by a compute with ID, I can include wildcard (see below)
           f_ID = per-atom vector calculated by a fix with ID
           f_ID[I] = Ith column of per-atom array calculated by a fix with ID, I can include wildcard (see below)
           v_name = per-atom vector calculated by an atom-style variable with name
           d_name = per-atom floating point vector with name, managed by fix property/atom
           i_name = per-atom integer vector with name, managed by fix property/atom

* *local* args = list of local attributes

  .. parsed-literal::

         possible attributes = index, c_ID, c_ID[I], f_ID, f_ID[I]
           index = enumeration of local values
           c_ID = local vector calculated by a compute with ID
           c_ID[I] = Ith column of local array calculated by a compute with ID, I can include wildcard (see below)
           f_ID = local vector calculated by a fix with ID
           f_ID[I] = Ith column of local array calculated by a fix with ID, I can include wildcard (see below)

Examples
""""""""

.. parsed-literal::

   dump myDump all atom 100 dump.atom
   dump myDump all atom/mpiio 100 dump.atom.mpiio
   dump myDump all atom/gz 100 dump.atom.gz
   dump 2 subgroup atom 50 dump.run.bin
   dump 2 subgroup atom 50 dump.run.mpiio.bin
   dump 4a all custom 100 dump.myforce.\* id type x y vx fx
   dump 4b flow custom 100 dump.%.myforce id type c_myF[3] v_ke
   dump 4b flow custom 100 dump.%.myforce id type c_myF[\*] v_ke
   dump 2 inner cfg 10 dump.snap.\*.cfg mass type xs ys zs vx vy vz
   dump snap all cfg 100 dump.config.\*.cfg mass type xs ys zs id type c_Stress[2]
   dump 1 all xtc 1000 file.xtc

Description
"""""""""""

Dump a snapshot of atom quantities to one or more files every N
timesteps in one of several styles.  The *image* and *movie* styles are
the exception: the *image* style renders a JPG, PNG, or PPM image file
of the atom configuration every N timesteps while the *movie* style
combines and compresses them into a movie file; both are discussed in
detail on the :doc:`dump image <dump_image>` doc page.  The timesteps on
which dump output is written can also be controlled by a variable.
See the :doc:`dump_modify every <dump_modify>` command.

Only information for atoms in the specified group is dumped.  The
:doc:`dump_modify thresh and region and refresh <dump_modify>` commands
can also alter what atoms are included.  Not all styles support
these options; see details on the :doc:`dump_modify <dump_modify>` doc
page.

As described below, the filename determines the kind of output (text
or binary or gzipped, one big file or one per timestep, one big file
or multiple smaller files).

.. note::

   Because periodic boundary conditions are enforced only on
   timesteps when neighbor lists are rebuilt, the coordinates of an atom
   written to a dump file may be slightly outside the simulation box.
   Re-neighbor timesteps will not typically coincide with the timesteps
   dump snapshots are written.  See the :doc:`dump_modify pbc <dump_modify>` command if you with to force coordinates to be
   strictly inside the simulation box.

.. note::

   Unless the :doc:`dump_modify sort <dump_modify>` option is
   invoked, the lines of atom information written to dump files
   (typically one line per atom) will be in an indeterminate order for
   each snapshot.  This is even true when running on a single processor,
   if the :doc:`atom_modify sort <atom_modify>` option is on, which it is
   by default.  In this case atoms are re-ordered periodically during a
   simulation, due to spatial sorting.  It is also true when running in
   parallel, because data for a single snapshot is collected from
   multiple processors, each of which owns a subset of the atoms.

For the *atom*\ , *custom*\ , *cfg*\ , and *local* styles, sorting is off by
default.  For the *dcd*\ , *xtc*\ , *xyz*\ , and *molfile* styles, sorting by
atom ID is on by default. See the :doc:`dump_modify <dump_modify>` doc
page for details.

The *atom/gz*\ , *cfg/gz*\ , *custom/gz*\ , and *xyz/gz* styles are identical
in command syntax to the corresponding styles without "gz", however,
they generate compressed files using the zlib library. Thus the filename
suffix ".gz" is mandatory. This is an alternative approach to writing
compressed files via a pipe, as done by the regular dump styles, which
may be required on clusters where the interface to the high-speed network
disallows using the fork() library call (which is needed for a pipe).
For the remainder of this doc page, you should thus consider the *atom*
and *atom/gz* styles (etc) to be inter-changeable, with the exception
of the required filename suffix.

As explained below, the *atom/mpiio*\ , *cfg/mpiio*\ , *custom/mpiio*\ , and
*xyz/mpiio* styles are identical in command syntax and in the format
of the dump files they create, to the corresponding styles without
"mpiio", except the single dump file they produce is written in
parallel via the MPI-IO library.  For the remainder of this doc page,
you should thus consider the *atom* and *atom/mpiio* styles (etc) to
be inter-changeable.  The one exception is how the filename is
specified for the MPI-IO styles, as explained below.

The precision of values output to text-based dump files can be
controlled by the :doc:`dump_modify format <dump_modify>` command and
its options.

----------

The *style* keyword determines what atom quantities are written to the
file and in what format.  Settings made via the
:doc:`dump_modify <dump_modify>` command can also alter the format of
individual values and the file itself.

The *atom*\ , *local*\ , and *custom* styles create files in a simple text
format that is self-explanatory when viewing a dump file.  Some of the
LAMMPS post-processing tools described on the :doc:`Tools <Tools>` doc
page, including `Pizza.py <http://www.sandia.gov/~sjplimp/pizza.html>`_,
work with this format, as does the :doc:`rerun <rerun>` command.

For post-processing purposes the *atom*\ , *local*\ , and *custom* text
files are self-describing in the following sense.

The dimensions of the simulation box are included in each snapshot.
For an orthogonal simulation box this information is formatted as:

.. parsed-literal::

   ITEM: BOX BOUNDS xx yy zz
   xlo xhi
   ylo yhi
   zlo zhi

where xlo,xhi are the maximum extents of the simulation box in the
x-dimension, and similarly for y and z.  The "xx yy zz" represent 6
characters that encode the style of boundary for each of the 6
simulation box boundaries (xlo,xhi and ylo,yhi and zlo,zhi).  Each of
the 6 characters is either p = periodic, f = fixed, s = shrink wrap,
or m = shrink wrapped with a minimum value.  See the
:doc:`boundary <boundary>` command for details.

For triclinic simulation boxes (non-orthogonal), an orthogonal
bounding box which encloses the triclinic simulation box is output,
along with the 3 tilt factors (xy, xz, yz) of the triclinic box,
formatted as follows:

.. parsed-literal::

   ITEM: BOX BOUNDS xy xz yz xx yy zz
   xlo_bound xhi_bound xy
   ylo_bound yhi_bound xz
   zlo_bound zhi_bound yz

The presence of the text "xy xz yz" in the ITEM line indicates that
the 3 tilt factors will be included on each of the 3 following lines.
This bounding box is convenient for many visualization programs.  The
meaning of the 6 character flags for "xx yy zz" is the same as above.

Note that the first two numbers on each line are now xlo\_bound instead
of xlo, etc, since they represent a bounding box.  See the :doc:`Howto triclinic <Howto_triclinic>` doc page for a geometric description
of triclinic boxes, as defined by LAMMPS, simple formulas for how the
6 bounding box extents (xlo\_bound,xhi\_bound,etc) are calculated from
the triclinic parameters, and how to transform those parameters to and
from other commonly used triclinic representations.

The "ITEM: ATOMS" line in each snapshot lists column descriptors for
the per-atom lines that follow.  For example, the descriptors would be
"id type xs ys zs" for the default *atom* style, and would be the atom
attributes you specify in the dump command for the *custom* style.

For style *atom*\ , atom coordinates are written to the file, along with
the atom ID and atom type.  By default, atom coords are written in a
scaled format (from 0 to 1).  I.e. an x value of 0.25 means the atom
is at a location 1/4 of the distance from xlo to xhi of the box
boundaries.  The format can be changed to unscaled coords via the
:doc:`dump_modify <dump_modify>` settings.  Image flags can also be
added for each atom via dump\_modify.

Style *custom* allows you to specify a list of atom attributes to be
written to the dump file for each atom.  Possible attributes are
listed above and will appear in the order specified.  You cannot
specify a quantity that is not defined for a particular simulation -
such as *q* for atom style *bond*\ , since that atom style doesn't
assign charges.  Dumps occur at the very end of a timestep, so atom
attributes will include effects due to fixes that are applied during
the timestep.  An explanation of the possible dump custom attributes
is given below.

For style *local*\ , local output generated by :doc:`computes <compute>`
and :doc:`fixes <fix>` is used to generate lines of output that is
written to the dump file.  This local data is typically calculated by
each processor based on the atoms it owns, but there may be zero or
more entities per atom, e.g. a list of bond distances.  An explanation
of the possible dump local attributes is given below.  Note that by
using input from the :doc:`compute property/local <compute_property_local>` command with dump local,
it is possible to generate information on bonds, angles, etc that can
be cut and pasted directly into a data file read by the
:doc:`read_data <read_data>` command.

Style *cfg* has the same command syntax as style *custom* and writes
extended CFG format files, as used by the
`AtomEye <http://mt.seas.upenn.edu/Archive/Graphics/A>`_ visualization
package.  Since the extended CFG format uses a single snapshot of the
system per file, a wildcard "\*" must be included in the filename, as
discussed below.  The list of atom attributes for style *cfg* must
begin with either "mass type xs ys zs" or "mass type xsu ysu zsu"
since these quantities are needed to write the CFG files in the
appropriate format (though the "mass" and "type" fields do not appear
explicitly in the file).  Any remaining attributes will be stored as
"auxiliary properties" in the CFG files.  Note that you will typically
want to use the :doc:`dump_modify element <dump_modify>` command with
CFG-formatted files, to associate element names with atom types, so
that AtomEye can render atoms appropriately. When unwrapped
coordinates *xsu*\ , *ysu*\ , and *zsu* are requested, the nominal AtomEye
periodic cell dimensions are expanded by a large factor UNWRAPEXPAND =
10.0, which ensures atoms that are displayed correctly for up to
UNWRAPEXPAND/2 periodic boundary crossings in any direction.  Beyond
this, AtomEye will rewrap the unwrapped coordinates.  The expansion
causes the atoms to be drawn farther away from the viewer, but it is
easy to zoom the atoms closer, and the interatomic distances are
unaffected.

The *dcd* style writes DCD files, a standard atomic trajectory format
used by the CHARMM, NAMD, and XPlor molecular dynamics packages.  DCD
files are binary and thus may not be portable to different machines.
The number of atoms per snapshot cannot change with the *dcd* style.
The *unwrap* option of the :doc:`dump_modify <dump_modify>` command
allows DCD coordinates to be written "unwrapped" by the image flags
for each atom.  Unwrapped means that if the atom has passed through
a periodic boundary one or more times, the value is printed for what
the coordinate would be if it had not been wrapped back into the
periodic box.  Note that these coordinates may thus be far outside
the box size stored with the snapshot.

The *xtc* style writes XTC files, a compressed trajectory format used
by the GROMACS molecular dynamics package, and described
`here <http://manual.gromacs.org/current/online/xtc.html>`_.
The precision used in XTC files can be adjusted via the
:doc:`dump_modify <dump_modify>` command.  The default value of 1000
means that coordinates are stored to 1/1000 nanometer accuracy.  XTC
files are portable binary files written in the NFS XDR data format,
so that any machine which supports XDR should be able to read them.
The number of atoms per snapshot cannot change with the *xtc* style.
The *unwrap* option of the :doc:`dump_modify <dump_modify>` command allows
XTC coordinates to be written "unwrapped" by the image flags for each
atom.  Unwrapped means that if the atom has passed through a periodic
boundary one or more times, the value is printed for what the
coordinate would be if it had not been wrapped back into the periodic
box.  Note that these coordinates may thus be far outside the box size
stored with the snapshot.

The *xyz* style writes XYZ files, which is a simple text-based
coordinate format that many codes can read. Specifically it has
a line with the number of atoms, then a comment line that is
usually ignored followed by one line per atom with the atom type
and the x-, y-, and z-coordinate of that atom. You can use the
:doc:`dump_modify element <dump_modify>` option to change the output
from using the (numerical) atom type to an element name (or some
other label). This will help many visualization programs to guess
bonds and colors.

Note that *atom*\ , *custom*\ , *dcd*\ , *xtc*\ , and *xyz* style dump files
can be read directly by `VMD <http://www.ks.uiuc.edu/Research/vmd>`_, a
popular molecular viewing program.

----------

Dumps are performed on timesteps that are a multiple of N (including
timestep 0) and on the last timestep of a minimization if the
minimization converges.  Note that this means a dump will not be
performed on the initial timestep after the dump command is invoked,
if the current timestep is not a multiple of N.  This behavior can be
changed via the :doc:`dump_modify first <dump_modify>` command, which
can also be useful if the dump command is invoked after a minimization
ended on an arbitrary timestep.  N can be changed between runs by
using the :doc:`dump_modify every <dump_modify>` command (not allowed
for *dcd* style).  The :doc:`dump_modify every <dump_modify>` command
also allows a variable to be used to determine the sequence of
timesteps on which dump files are written.  In this mode a dump on the
first timestep of a run will also not be written unless the
:doc:`dump_modify first <dump_modify>` command is used.

The specified filename determines how the dump file(s) is written.
The default is to write one large text file, which is opened when the
dump command is invoked and closed when an :doc:`undump <undump>`
command is used or when LAMMPS exits.  For the *dcd* and *xtc* styles,
this is a single large binary file.

Dump filenames can contain two wildcard characters.  If a "\*"
character appears in the filename, then one file per snapshot is
written and the "\*" character is replaced with the timestep value.
For example, tmp.dump.\* becomes tmp.dump.0, tmp.dump.10000,
tmp.dump.20000, etc.  This option is not available for the *dcd* and
*xtc* styles.  Note that the :doc:`dump_modify pad <dump_modify>`
command can be used to insure all timestep numbers are the same length
(e.g. 00010), which can make it easier to read a series of dump files
in order with some post-processing tools.

If a "%" character appears in the filename, then each of P processors
writes a portion of the dump file, and the "%" character is replaced
with the processor ID from 0 to P-1.  For example, tmp.dump.% becomes
tmp.dump.0, tmp.dump.1, ... tmp.dump.P-1, etc.  This creates smaller
files and can be a fast mode of output on parallel machines that
support parallel I/O for output. This option is not available for the
*dcd*\ , *xtc*\ , and *xyz* styles.

By default, P = the number of processors meaning one file per
processor, but P can be set to a smaller value via the *nfile* or
*fileper* keywords of the :doc:`dump_modify <dump_modify>` command.
These options can be the most efficient way of writing out dump files
when running on large numbers of processors.

Note that using the "\*" and "%" characters together can produce a
large number of small dump files!

For the *atom/mpiio*\ , *cfg/mpiio*\ , *custom/mpiio*\ , and *xyz/mpiio*
styles, a single dump file is written in parallel via the MPI-IO
library, which is part of the MPI standard for versions 2.0 and above.
Using MPI-IO requires two steps.  First, build LAMMPS with its MPIIO
package installed, e.g.

.. parsed-literal::

   make yes-mpiio    # installs the MPIIO package
   make mpi          # build LAMMPS for your platform

Second, use a dump filename which contains ".mpiio".  Note that it
does not have to end in ".mpiio", just contain those characters.
Unlike MPI-IO restart files, which must be both written and read using
MPI-IO, the dump files produced by these MPI-IO styles are identical
in format to the files produced by their non-MPI-IO style
counterparts.  This means you can write a dump file using MPI-IO and
use the :doc:`read_dump <read_dump>` command or perform other
post-processing, just as if the dump file was not written using
MPI-IO.

Note that MPI-IO dump files are one large file which all processors
write to.  You thus cannot use the "%" wildcard character described
above in the filename since that specifies generation of multiple
files.  You can use the ".bin" suffix described below in an MPI-IO
dump file; again this file will be written in parallel and have the
same binary format as if it were written without MPI-IO.

If the filename ends with ".bin", the dump file (or files, if "\*" or
"%" is also used) is written in binary format.  A binary dump file
will be about the same size as a text version, but will typically
write out much faster.  Of course, when post-processing, you will need
to convert it back to text format (see the :ref:`binary2txt tool <binary>`) or write your own code to read the binary
file.  The format of the binary file can be understood by looking at
the tools/binary2txt.cpp file.  This option is only available for the
*atom* and *custom* styles.

If the filename ends with ".gz", the dump file (or files, if "\*" or "%"
is also used) is written in gzipped format.  A gzipped dump file will
be about 3x smaller than the text version, but will also take longer
to write.  This option is not available for the *dcd* and *xtc*
styles.

----------

Note that in the discussion which follows, for styles which can
reference values from a compute or fix, like the *custom*\ , *cfg*\ , or
*local* styles, the bracketed index I can be specified using a
wildcard asterisk with the index to effectively specify multiple
values.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
size of the vector (for *mode* = scalar) or the number of columns in
the array (for *mode* = vector), then an asterisk with no numeric
values means all indices from 1 to N.  A leading asterisk means all
indices from 1 to n (inclusive).  A trailing asterisk means all
indices from n to N (inclusive).  A middle asterisk means all indices
from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 dump commands are
equivalent, since the :doc:`compute stress/atom <compute_stress_atom>`
command creates a per-atom array with 6 columns:

.. parsed-literal::

   compute myPress all stress/atom NULL
   dump 2 all custom 100 tmp.dump id myPress[\*]
   dump 2 all custom 100 tmp.dump id myPress[1] myPress[2] myPress[3] &
                                     myPress[4] myPress[5] myPress[6]

----------

This section explains the local attributes that can be specified as
part of the *local* style.

The *index* attribute can be used to generate an index number from 1
to N for each line written into the dump file, where N is the total
number of local datums from all processors, or lines of output that
will appear in the snapshot.  Note that because data from different
processors depend on what atoms they currently own, and atoms migrate
between processor, there is no guarantee that the same index will be
used for the same info (e.g. a particular bond) in successive
snapshots.

The *c\_ID* and *c\_ID[I]* attributes allow local vectors or arrays
calculated by a :doc:`compute <compute>` to be output.  The ID in the
attribute should be replaced by the actual ID of the compute that has
been defined previously in the input script.  See the
:doc:`compute <compute>` command for details.  There are computes for
calculating local information such as indices, types, and energies for
bonds and angles.

Note that computes which calculate global or per-atom quantities, as
opposed to local quantities, cannot be output in a dump local command.
Instead, global quantities can be output by the :doc:`thermo_style custom <thermo_style>` command, and per-atom quantities can be
output by the dump custom command.

If *c\_ID* is used as a attribute, then the local vector calculated by
the compute is printed.  If *c\_ID[I]* is used, then I must be in the
range from 1-M, which will print the Ith column of the local array
with M columns calculated by the compute.  See the discussion above
for how I can be specified with a wildcard asterisk to effectively
specify multiple values.

The *f\_ID* and *f\_ID[I]* attributes allow local vectors or arrays
calculated by a :doc:`fix <fix>` to be output.  The ID in the attribute
should be replaced by the actual ID of the fix that has been defined
previously in the input script.

If *f\_ID* is used as a attribute, then the local vector calculated by
the fix is printed.  If *f\_ID[I]* is used, then I must be in the
range from 1-M, which will print the Ith column of the local with M
columns calculated by the fix.  See the discussion above for how I can
be specified with a wildcard asterisk to effectively specify multiple
values.

Here is an example of how to dump bond info for a system, including
the distance and energy of each bond:

.. parsed-literal::

   compute 1 all property/local batom1 batom2 btype
   compute 2 all bond/local dist eng
   dump 1 all local 1000 tmp.dump index c_1[1] c_1[2] c_1[3] c_2[1] c_2[2]

----------

This section explains the atom attributes that can be specified as
part of the *custom* and *cfg* styles.

The *id*\ , *mol*\ , *proc*\ , *procp1*\ , *type*\ , *element*\ , *mass*\ , *vx*\ ,
*vy*\ , *vz*\ , *fx*\ , *fy*\ , *fz*\ , *q* attributes are self-explanatory.

*Id* is the atom ID.  *Mol* is the molecule ID, included in the data
file for molecular systems.  *Proc* is the ID of the processor (0 to
Nprocs-1) that currently owns the atom.  *Procp1* is the proc ID+1,
which can be convenient in place of a *type* attribute (1 to Ntypes)
for coloring atoms in a visualization program.  *Type* is the atom
type (1 to Ntypes).  *Element* is typically the chemical name of an
element, which you must assign to each type via the :doc:`dump_modify element <dump_modify>` command.  More generally, it can be any
string you wish to associated with an atom type.  *Mass* is the atom
mass.  *Vx*\ , *vy*\ , *vz*\ , *fx*\ , *fy*\ , *fz*\ , and *q* are components of
atom velocity and force and atomic charge.

There are several options for outputting atom coordinates.  The *x*\ ,
*y*\ , *z* attributes write atom coordinates "unscaled", in the
appropriate distance :doc:`units <units>` (Angstroms, sigma, etc).  Use
*xs*\ , *ys*\ , *zs* if you want the coordinates "scaled" to the box size,
so that each value is 0.0 to 1.0.  If the simulation box is triclinic
(tilted), then all atom coords will still be between 0.0 and 1.0.
I.e. actual unscaled (x,y,z) = xs\*A + ys\*B + zs\*C, where (A,B,C) are
the non-orthogonal vectors of the simulation box edges, as discussed
on the :doc:`Howto triclinic <Howto_triclinic>` doc page.

Use *xu*\ , *yu*\ , *zu* if you want the coordinates "unwrapped" by the
image flags for each atom.  Unwrapped means that if the atom has
passed through a periodic boundary one or more times, the value is
printed for what the coordinate would be if it had not been wrapped
back into the periodic box.  Note that using *xu*\ , *yu*\ , *zu* means
that the coordinate values may be far outside the box bounds printed
with the snapshot.  Using *xsu*\ , *ysu*\ , *zsu* is similar to using
*xu*\ , *yu*\ , *zu*\ , except that the unwrapped coordinates are scaled by
the box size. Atoms that have passed through a periodic boundary will
have the corresponding coordinate increased or decreased by 1.0.

The image flags can be printed directly using the *ix*\ , *iy*\ , *iz*
attributes.  For periodic dimensions, they specify which image of the
simulation box the atom is considered to be in.  An image of 0 means
it is inside the box as defined.  A value of 2 means add 2 box lengths
to get the true value.  A value of -1 means subtract 1 box length to
get the true value.  LAMMPS updates these flags as atoms cross
periodic boundaries during the simulation.

The *mux*\ , *muy*\ , *muz* attributes are specific to dipolar systems
defined with an atom style of *dipole*\ .  They give the orientation of
the atom's point dipole moment.  The *mu* attribute gives the
magnitude of the atom's dipole moment.

The *radius* and *diameter* attributes are specific to spherical
particles that have a finite size, such as those defined with an atom
style of *sphere*\ .

The *omegax*\ , *omegay*\ , and *omegaz* attributes are specific to
finite-size spherical particles that have an angular velocity.  Only
certain atom styles, such as *sphere* define this quantity.

The *angmomx*\ , *angmomy*\ , and *angmomz* attributes are specific to
finite-size aspherical particles that have an angular momentum.  Only
the *ellipsoid* atom style defines this quantity.

The *tqx*\ , *tqy*\ , *tqz* attributes are for finite-size particles that
can sustain a rotational torque due to interactions with other
particles.

The *c\_ID* and *c\_ID[I]* attributes allow per-atom vectors or arrays
calculated by a :doc:`compute <compute>` to be output.  The ID in the
attribute should be replaced by the actual ID of the compute that has
been defined previously in the input script.  See the
:doc:`compute <compute>` command for details.  There are computes for
calculating the per-atom energy, stress, centro-symmetry parameter,
and coordination number of individual atoms.

Note that computes which calculate global or local quantities, as
opposed to per-atom quantities, cannot be output in a dump custom
command.  Instead, global quantities can be output by the
:doc:`thermo_style custom <thermo_style>` command, and local quantities
can be output by the dump local command.

If *c\_ID* is used as a attribute, then the per-atom vector calculated
by the compute is printed.  If *c\_ID[I]* is used, then I must be in
the range from 1-M, which will print the Ith column of the per-atom
array with M columns calculated by the compute.  See the discussion
above for how I can be specified with a wildcard asterisk to
effectively specify multiple values.

The *f\_ID* and *f\_ID[I]* attributes allow vector or array per-atom
quantities calculated by a :doc:`fix <fix>` to be output.  The ID in the
attribute should be replaced by the actual ID of the fix that has been
defined previously in the input script.  The :doc:`fix ave/atom <fix_ave_atom>` command is one that calculates per-atom
quantities.  Since it can time-average per-atom quantities produced by
any :doc:`compute <compute>`, :doc:`fix <fix>`, or atom-style
:doc:`variable <variable>`, this allows those time-averaged results to
be written to a dump file.

If *f\_ID* is used as a attribute, then the per-atom vector calculated
by the fix is printed.  If *f\_ID[I]* is used, then I must be in the
range from 1-M, which will print the Ith column of the per-atom array
with M columns calculated by the fix.  See the discussion above for
how I can be specified with a wildcard asterisk to effectively specify
multiple values.

The *v\_name* attribute allows per-atom vectors calculated by a
:doc:`variable <variable>` to be output.  The name in the attribute
should be replaced by the actual name of the variable that has been
defined previously in the input script.  Only an atom-style variable
can be referenced, since it is the only style that generates per-atom
values.  Variables of style *atom* can reference individual atom
attributes, per-atom attributes, thermodynamic keywords, or
invoke other computes, fixes, or variables when they are evaluated, so
this is a very general means of creating quantities to output to a
dump file.

The *d\_name* and *i\_name* attributes allow to output custom per atom
floating point or integer properties that are managed by
:doc:`fix property/atom <fix_property_atom>`.

See the :doc:`Modify <Modify>` doc page for information on how to add
new compute and fix styles to LAMMPS to calculate per-atom quantities
which could then be output into dump files.

----------

Restrictions
""""""""""""

To write gzipped dump files, you must either compile LAMMPS with the
-DLAMMPS\_GZIP option or use the styles from the COMPRESS package.
See the :doc:`Build settings <Build_settings>` doc page for details.

The *atom/gz*\ , *cfg/gz*\ , *custom/gz*\ , and *xyz/gz* styles are part of
the COMPRESS package.  They are only enabled if LAMMPS was built with
that package.  See the :doc:`Build package <Build_package>` doc page for
more info.

The *atom/mpiio*\ , *cfg/mpiio*\ , *custom/mpiio*\ , and *xyz/mpiio* styles
are part of the MPIIO package.  They are only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

The *xtc* style is part of the MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`dump atom/adios <dump_adios>`, :doc:`dump custom/adios <dump_adios>`,
:doc:`dump h5md <dump_h5md>`, :doc:`dump image <dump_image>`,
:doc:`dump molfile <dump_molfile>`, :doc:`dump_modify <dump_modify>`,
:doc:`undump <undump>`

Default
"""""""

The defaults for the *image* and *movie* styles are listed on the
:doc:`dump image <dump_image>` doc page.

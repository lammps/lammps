.. index:: read_dump

read_dump command
=================

Syntax
""""""

.. code-block:: LAMMPS

   read_dump file Nstep field1 field2 ... keyword values ...

* file = name of dump file to read
* Nstep = snapshot timestep to read from file
* one or more fields may be appended

  .. parsed-literal::

     field = *x* or *y* or *z* or *vx* or *vy* or *vz* or *q* or *ix* or *iy* or *iz* or *fx* or *fy* or *fz*
       *x*,\ *y*,\ *z* = atom coordinates
       *vx*,\ *vy*,\ *vz* = velocity components
       *q* = charge
       *ix*,\ *iy*,\ *iz* = image flags in each dimension
       *fx*,\ *fy*,\ *fz* = force components

* zero or more keyword/value pairs may be appended
* keyword = *nfile* or *box* or *timestep* or *replace* or *purge* or *trim* or *add* or *label* or *scaled* or *wrapped* or *format*

  .. parsed-literal::

       *nfile* value = Nfiles = how many parallel dump files exist
       *box* value = *yes* or *no* = replace simulation box with dump box
       *timestep* value = *yes* or *no* = reset simulation timestep with dump timestep
       *replace* value = *yes* or *no* = overwrite atoms with dump atoms
       *purge* value = *yes* or *no* = delete all atoms before adding dump atoms
       *trim* value = *yes* or *no* = trim atoms not in dump snapshot
       *add* value = *yes* or *keep* or *no* = add new dump atoms to system
       *label* value = field column
         field = one of the listed fields or *id* or *type*
         column = label on corresponding column in dump file
       *scaled* value = *yes* or *no* = coords in dump file are scaled/unscaled
       *wrapped* value = *yes* or *no* = coords in dump file are wrapped/unwrapped
       *format* values = format of dump file, must be last keyword if used
         *native* = native LAMMPS dump file
         *xyz* = XYZ file
         *adios* [*timeout* value] = dump file written by the :doc:`dump adios <dump_adios>` command
           *timeout* = specify waiting time for the arrival of the timestep when running concurrently.
                     The value is a float number and is interpreted in seconds.
         *molfile* style path = VMD molfile plugin interface
           style = *dcd* or *xyz* or others supported by molfile plugins
           path = optional path for location of molfile plugins

Examples
""""""""

.. code-block:: LAMMPS

   read_dump dump.file 5000 x y z
   read_dump dump.xyz 5 x y z box no format xyz
   read_dump dump.xyz 10 x y z box no format molfile xyz "../plugins"
   read_dump dump.dcd 0 x y z box yes format molfile dcd
   read_dump dump.file 1000 x y z vx vy vz box yes format molfile lammpstrj /usr/local/lib/vmd/plugins/LINUXAMD64/plugins/molfile
   read_dump dump.file 5000 x y vx vy trim yes
   read_dump dump.file 5000 x y vx vy add yes box no timestep no
   read_dump ../run7/dump.file.gz 10000 x y z box yes
   read_dump dump.xyz 10 x y z box no format molfile xyz ../plugins
   read_dump dump.dcd 0 x y z format molfile dcd
   read_dump dump.file 1000 x y z vx vy vz format molfile lammpstrj /usr/local/lib/vmd/plugins/LINUXAMD64/plugins/molfile
   read_dump dump.bp 5000 x y z vx vy vz format adios
   read_dump dump.bp 5000 x y z vx vy vz format adios timeout 60.0

Description
"""""""""""

Read atom information from a dump file to overwrite the current atom
coordinates, and optionally the atom velocities and image flags, the
simulation timestep, and the simulation box dimensions.  This is useful
for restarting a run from a particular snapshot in a dump file.  See the
:doc:`read_restart <read_restart>` and :doc:`read_data <read_data>`
commands for alternative methods to do this.  Also see the
:doc:`rerun <rerun>` command for a means of reading multiple snapshots
from a dump file.

Note that a simulation box must already be defined before using the
read_dump command.  This can be done by the
:doc:`create_box <create_box>`, :doc:`read_data <read_data>`, or
:doc:`read_restart <read_restart>` commands.  The read_dump command can
reset the simulation box dimensions, as explained below.

Also note that reading per-atom information from a dump snapshot is
limited to the atom coordinates, velocities and image flags, as
explained below.  Other atom properties, which may be necessary to run
a valid simulation, such as atom charge, or bond topology information
for a molecular system, are not read from (or may not even be contained
in) dump files.  Thus this auxiliary information should be defined in
the usual way, e.g. in a data file read in by a :doc:`read_data <read_data>`
command, before using the read_dump command, or by the :doc:`set <set>`
command, after the dump snapshot is read.

----------

If the dump filename specified as *file* ends with ".gz", the dump
file is read in gzipped format.

You can read dump files that were written (in parallel) to multiple
files via the "%" wild-card character in the dump file name.  If any
specified dump file name contains a "%", they must all contain it.
See the :doc:`dump <dump>` command for details.
The "%" wild-card character is only supported by the *native* format
for dump files, described next.

If reading parallel dump files, you must also use the *nfile* keyword
to tell LAMMPS how many parallel files exist, via its specified
*Nfiles* value.

The format of the dump file is selected through the *format* keyword.
If specified, it must be the last keyword used, since all remaining
arguments are passed on to the dump reader.  The *native* format is
for native LAMMPS dump files, written with a :doc:`dump atom <dump>`
or :doc:`dump custom <dump>` command.  The *xyz* format is for generic XYZ
formatted dump files.  These formats take no additional values.

The *molfile* format supports reading data through using the `VMD <vmd_>`_
molfile plugin interface. This dump reader format is only available,
if the MOLFILE package has been installed when compiling
LAMMPS.

The *molfile* format takes one or two additional values.  The *style*
value determines the file format to be used and can be any format that
the molfile plugins support, such as DCD or XYZ.  Note that DCD dump
files can be written by LAMMPS via the :doc:`dump dcd <dump>` command.
The *path* value specifies a list of directories which LAMMPS will
search for the molfile plugins appropriate to the specified *style*\ .
The syntax of the *path* value is like other search paths: it can
contain multiple directories separated by a colon (or semi-colon on
windows).  The *path* keyword is optional and defaults to ".",
i.e. the current directory.

The *adios* format supports reading data that was written by the
:doc:`dump adios <dump_adios>` command. The
entire dump is read in parallel across all the processes, dividing
the atoms evenly among the processes. The number of writers that
has written the dump file does not matter. Using the adios style for
dump and read_dump is a convenient way to dump all atoms from *N*
writers and read it back by *M* readers. If one is running two
LAMMPS instances concurrently where one dumps data and the other is
reading it with the rerun command, the timeout option can be specified
to wait on the reader side for the arrival of the requested step.

Support for other dump format readers may be added in the future.

----------

Global information is first read from the dump file, namely timestep
and box information.

The dump file is scanned for a snapshot with a timestamp that matches
the specified *Nstep*\ .  This means the LAMMPS timestep the dump file
snapshot was written on for the *native* or *adios* formats.

The list of timestamps available in an adios .bp file is stored in the
variable *ntimestep*:

.. parsed-literal:: console

  $ bpls dump.bp -d ntimestep
    uint64_t  ntimestep  5*scalar
      (0)    0 50 100 150 200

Note that the *xyz* and *molfile* formats do not store the timestep.
For these formats, timesteps are numbered logically, in a sequential
manner, starting from 0.  Thus to access the 10th snapshot in an *xyz*
or *mofile* formatted dump file, use *Nstep* = 9.

The dimensions of the simulation box for the selected snapshot are
also read; see the *box* keyword discussion below.  For the *native*
format, an error is generated if the snapshot is for a triclinic box
and the current simulation box is orthogonal or vice versa.  A warning
will be generated if the snapshot box boundary conditions (periodic,
shrink-wrapped, etc) do not match the current simulation boundary
conditions, but the boundary condition information in the snapshot is
otherwise ignored.  See the "boundary" command for more details. The
*adios* reader does the same as the *native* format reader.

For the *xyz* format, no information about the box is available, so
you must set the *box* flag to *no*\ .  See details below.

For the *molfile* format, reading simulation box information is
typically supported, but the location of the simulation box origin is
lost and no explicit information about periodicity or
orthogonal/triclinic box shape is available.  The MOLFILE package
makes a best effort to guess based on heuristics, but this may not
always work perfectly.

----------

Per-atom information from the dump file snapshot is then read from the
dump file snapshot.  This corresponds to the specified *fields* listed
in the read_dump command.  It is an error to specify a z-dimension
field, namely *z*, *vz*, or *iz*, for a 2d simulation.

For dump files in *native* format, each column of per-atom data has a
text label listed in the file.  A matching label for each field must
appear, e.g. the label "vy" for the field *vy*\ .  For the *x*, *y*, *z*
fields any of the following labels are considered a match:

.. parsed-literal::

   x, xs, xu, xsu for field *x*
   y, ys, yu, ysu for field *y*
   z, zs, zu, zsu for field *z*

The meaning of xs (scaled), xu (unwrapped), and xsu (scaled and
unwrapped) is explained on the :doc:`dump <dump>` command doc page.
These labels are searched for in the list of column labels in the dump
file, in order, until a match is found.

The dump file must also contain atom IDs, with a column label of "id".

If the *add* keyword is specified with a value of *yes* or *keep*, as
discussed below, the dump file must contain atom types, with a column
label of "type".

If a column label you want to read from the dump file is not a match
to a specified field, the *label* keyword can be used to specify the
specific column label from the dump file to associate with that field.
An example is if a time-averaged coordinate is written to the dump
file via the :doc:`fix ave/atom <fix_ave_atom>` command.  The column
will then have a label corresponding to the fix-ID rather than "x" or
"xs".  The *label* keyword can also be used to specify new column
labels for fields *id* and *type*\ .

For dump files in *xyz* format, only the *x*, *y*, and *z* fields are
supported.  The dump file does not store atom IDs, so these are
assigned consecutively to the atoms as they appear in the dump file,
starting from 1.  Thus you should ensure that order of atoms is
consistent from snapshot to snapshot in the XYZ dump file.  See
the :doc:`dump_modify sort <dump_modify>` command if the XYZ dump file
was written by LAMMPS.

For dump files in *molfile* format, the *x*, *y*, *z*, *vx*, *vy*, and
*vz* fields can be specified.  However, not all molfile formats store
velocities, or their respective plugins may not support reading of
velocities.  The molfile dump files do not store atom IDs, so these
are assigned consecutively to the atoms as they appear in the dump
file, starting from 1.  Thus you should ensure that order of atoms are
consistent from snapshot to snapshot in the molfile dump file.
See the :doc:`dump_modify sort <dump_modify>` command if the dump file
was written by LAMMPS.

The *adios* format supports all fields that the *native* format supports
except for the *q* charge field.
The list of fields stored in an adios .bp file is recorded in the attributes
*columns* (array of short strings) and *columnstr* (space-separated single string).

.. parsed-literal:: console

  $ bpls -la dump.bp column*
    string    columns            attr   = {"id", "type", "x", "y", "z", "vx", "vy", "vz"}
    string    columnstr          attr   = "id type x y z vx vy vz "

----------

Information from the dump file snapshot is used to overwrite or
replace properties of the current system.  There are various options
for how this is done, determined by the specified fields and optional
keywords.

.. versionchanged:: 3Aug2022

The timestep of the snapshot becomes the current timestep for the
simulation unless the *timestep* keyword is specified with a *no* value
(default setting is *yes*).  See the :doc:`reset_timestep <reset_timestep>`
command if you wish to change this to a different value after the dump
snapshot is read.

If the *box* keyword is specified with a *yes* value, then the current
simulation box dimensions are replaced by the dump snapshot box
dimensions.  If the *box* keyword is specified with a *no* value, the
current simulation box is unchanged.

If the *purge* keyword is specified with a *yes* value, then all
current atoms in the system are deleted before any of the operations
invoked by the *replace*, *trim*, or *add* keywords take place.

If the *replace* keyword is specified with a *yes* value, then atoms
with IDs that are in both the current system and the dump snapshot
have their properties overwritten by field values.  If the *replace*
keyword is specified with a *no* value, atoms with IDs that are in
both the current system and the dump snapshot are not modified.

If the *trim* keyword is specified with a *yes* value, then atoms with
IDs that are in the current system but not in the dump snapshot are
deleted.  These atoms are unaffected if the *trim* keyword is
specified with a *no* value.

If the *add* keyword is specified with a *no* value (default), then
dump file atoms with IDs that are not in the current system are not
added to the system.  They are simply ignored.

If a *yes* value is specified, the atoms with new IDs are added to the
system but their atom IDs are not preserved.  Instead, after all the
atoms are added, new IDs are assigned to them in the same manner as is
described for the :doc:`create_atoms <create_atoms>` command.  Basically
the largest existing atom ID in the system is identified, and all the
added atoms are assigned IDs that consecutively follow the largest ID.

If a *keep* value is specified, the atoms with new IDs are added to
the system and their atom IDs are preserved.  This may lead to
non-contiguous IDs for the combined system.

Note that atoms added via the *add* keyword will only have the
attributes read from the dump file due to the *field* arguments.  For
example, if *x* or *y* or *z* or *q* is not specified as a field, a
value of 0.0 is used for added atoms.  Added atoms must have an atom
type, so this value must appear in the dump file.

Any other attributes (e.g. charge or particle diameter for spherical
particles) will be set to default values, the same as if the
:doc:`create_atoms <create_atoms>` command were used.

----------

Atom coordinates read from the dump file are first converted into
unscaled coordinates, relative to the box dimensions of the snapshot.
These coordinates are then be assigned to an existing or new atom in
the current simulation.  The coordinates will then be remapped to the
simulation box, whether it is the original box or the dump snapshot
box.  If periodic boundary conditions apply, this means the atom will
be remapped back into the simulation box if necessary.  If shrink-wrap
boundary conditions apply, the new coordinates may change the
simulation box dimensions.  If fixed boundary conditions apply, the
atom will be lost if it is outside the simulation box.

For *native* format dump files, the 3 xyz image flags for an atom in
the dump file are set to the corresponding values appearing in the
dump file if the *ix*, *iy*, *iz* fields are specified.  If not
specified, the image flags for replaced atoms are not changed and
image flags for new atoms are set to default values.  If coordinates
read from the dump file are in unwrapped format (e.g. *xu*\ ) then the
image flags for read-in atoms are also set to default values.  The
remapping procedure described in the previous paragraph will then
change images flags for all atoms (old and new) if periodic boundary
conditions are applied to remap an atom back into the simulation box.

.. note::

   If you get a warning about inconsistent image flags after
   reading in a dump snapshot, it means one or more pairs of bonded atoms
   now have inconsistent image flags.  As discussed on the :doc:`Errors common <Errors_common>` page this may or may not cause problems
   for subsequent simulations.  One way this can happen is if you read
   image flag fields from the dump file but do not also use the dump file
   box parameters.

LAMMPS knows how to compute unscaled and remapped coordinates for the
snapshot column labels discussed above, e.g. *x*, *xs*, *xu*, *xsu*\ .
If another column label is assigned to the *x* or *y* or *z* field via
the *label* keyword, e.g. for coordinates output by the :doc:`fix ave/atom <fix_ave_atom>` command, then LAMMPS needs to know whether
the coordinate information in the dump file is scaled and/or wrapped.
This can be set via the *scaled* and *wrapped* keywords.  Note that
the value of the *scaled* and *wrapped* keywords is ignored for fields
*x* or *y* or *z* if the *label* keyword is not used to assign a
column label to that field.

The scaled/unscaled and wrapped/unwrapped setting must be identical
for any of the *x*, *y*, *z* fields that are specified.  Thus you
cannot read *xs* and *yu* from the dump file.  Also, if the dump file
coordinates are scaled and the simulation box is triclinic, then all 3
of the *x*, *y*, *z* fields must be specified, since they are all
needed to generate absolute, unscaled coordinates.

----------

Restrictions
""""""""""""

To read gzipped dump files, you must compile LAMMPS with the
-DLAMMPS_GZIP option.  See the :doc:`Build settings <Build_settings>`
doc page for details.

The *molfile* dump file formats are part of the MOLFILE package.
They are only enabled if LAMMPS was built with that packages.  See the
:doc:`Build package <Build_package>` page for more info.

To write and read adios .bp files, you must compile LAMMPS with the
:ref:`ADIOS <PKG-ADIOS>` package.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump molfile <dump_molfile>`,
:doc:`dump adios <dump_adios>`,
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`,
:doc:`rerun <rerun>`

Default
"""""""

The option defaults are box = yes, timestep = yes, replace = yes, purge = no,
trim = no, add = no, scaled = no, wrapped = yes, and format = native.

.. _vmd: https://www.ks.uiuc.edu/Research/vmd

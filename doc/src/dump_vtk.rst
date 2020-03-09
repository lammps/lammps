.. index:: dump vtk

dump vtk command
================

Syntax
""""""


.. parsed-literal::

   dump ID group-ID vtk N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* vtk = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = same as arguments for :doc:`dump_style custom <dump>`

Examples
""""""""


.. parsed-literal::

   dump dmpvtk all vtk 100 dump\*.myforce.vtk id type vx fx
   dump dmpvtp flow vtk 100 dump\*.%.displace.vtp id type c_myD[1] c_myD[2] c_myD[3] v_ke

Description
"""""""""""

Dump a snapshot of atom quantities to one or more files every N
timesteps in a format readable by the `VTK visualization toolkit <http://www.vtk.org>`_ or other visualization tools that use it,
e.g. `ParaView <http://www.paraview.org>`_.  The timesteps on which dump
output is written can also be controlled by a variable; see the
:doc:`dump_modify every <dump_modify>` command for details.

This dump style is similar to :doc:`dump_style custom <dump>` but uses
the VTK library to write data to VTK simple legacy or XML format
depending on the filename extension specified for the dump file.  This
can be either *\*.vtk* for the legacy format or *\*.vtp* and *\*.vtu*,
respectively, for XML format; see the `VTK homepage <http://www.vtk.org/VTK/img/file-formats.pdf>`_ for a detailed
description of these formats.  Since this naming convention conflicts
with the way binary output is usually specified (see below), the
:doc:`dump_modify binary <dump_modify>` command allows setting of a
binary option for this dump style explicitly.

Only information for atoms in the specified group is dumped.  The
:doc:`dump_modify thresh and region <dump_modify>` commands can also
alter what atoms are included; see details below.

As described below, special characters ("\*", "%") in the filename
determine the kind of output.

.. warning::

   Because periodic boundary conditions are enforced only
   on timesteps when neighbor lists are rebuilt, the coordinates of an
   atom written to a dump file may be slightly outside the simulation
   box.

.. warning::

   Unless the :doc:`dump_modify sort <dump_modify>` option
   is invoked, the lines of atom information written to dump files will
   be in an indeterminate order for each snapshot.  This is even true
   when running on a single processor, if the :doc:`atom_modify sort <atom_modify>` option is on, which it is by default.  In this
   case atoms are re-ordered periodically during a simulation, due to
   spatial sorting.  It is also true when running in parallel, because
   data for a single snapshot is collected from multiple processors, each
   of which owns a subset of the atoms.

For the *vtk* style, sorting is off by default. See the
:doc:`dump_modify <dump_modify>` doc page for details.


----------


The dimensions of the simulation box are written to a separate file
for each snapshot (either in legacy VTK or XML format depending on the
format of the main dump file) with the suffix *\_boundingBox* appended
to the given dump filename.

For an orthogonal simulation box this information is saved as a
rectilinear grid (legacy .vtk or .vtr XML format).

Triclinic simulation boxes (non-orthogonal) are saved as
hexahedrons in either legacy .vtk or .vtu XML format.

Style *vtk* allows you to specify a list of atom attributes to be
written to the dump file for each atom.  The list of possible attributes
is the same as for the :doc:`dump_style custom <dump>` command; see
its doc page for a listing and an explanation of each attribute.

.. note::

   Since position data is required to write VTK files the atom
   attributes "x y z" do not have to be specified explicitly; they will
   be included in the dump file regardless.  Also, in contrast to the
   *custom* style, the specified *vtk* attributes are rearranged to
   ensure correct ordering of vector components (except for computes and
   fixes - these have to be given in the right order) and duplicate
   entries are removed.

The VTK format uses a single snapshot of the system per file, thus
a wildcard "\*" must be included in the filename, as discussed below.
Otherwise the dump files will get overwritten with the new snapshot
each time.


----------


Dumps are performed on timesteps that are a multiple of N (including
timestep 0) and on the last timestep of a minimization if the
minimization converges.  Note that this means a dump will not be
performed on the initial timestep after the dump command is invoked,
if the current timestep is not a multiple of N.  This behavior can be
changed via the :doc:`dump_modify first <dump_modify>` command, which
can also be useful if the dump command is invoked after a minimization
ended on an arbitrary timestep.  N can be changed between runs by
using the :doc:`dump_modify every <dump_modify>` command.
The :doc:`dump_modify every <dump_modify>` command
also allows a variable to be used to determine the sequence of
timesteps on which dump files are written.  In this mode a dump on the
first timestep of a run will also not be written unless the
:doc:`dump_modify first <dump_modify>` command is used.

Dump filenames can contain two wildcard characters.  If a "\*"
character appears in the filename, then one file per snapshot is
written and the "\*" character is replaced with the timestep value.
For example, tmp.dump\*.vtk becomes tmp.dump0.vtk, tmp.dump10000.vtk,
tmp.dump20000.vtk, etc.  Note that the :doc:`dump_modify pad <dump_modify>`
command can be used to insure all timestep numbers are the same length
(e.g. 00010), which can make it easier to read a series of dump files
in order with some post-processing tools.

If a "%" character appears in the filename, then each of P processors
writes a portion of the dump file, and the "%" character is replaced
with the processor ID from 0 to P-1 preceded by an underscore character.
For example, tmp.dump%.vtp becomes tmp.dump\_0.vtp, tmp.dump\_1.vtp, ...
tmp.dump\_P-1.vtp, etc.  This creates smaller files and can be a fast
mode of output on parallel machines that support parallel I/O for output.

By default, P = the number of processors meaning one file per
processor, but P can be set to a smaller value via the *nfile* or
*fileper* keywords of the :doc:`dump_modify <dump_modify>` command.
These options can be the most efficient way of writing out dump files
when running on large numbers of processors.

For the legacy VTK format "%" is ignored and P = 1, i.e., only
processor 0 does write files.

Note that using the "\*" and "%" characters together can produce a
large number of small dump files!

If *dump\_modify binary* is used, the dump file (or files, if "\*" or
"%" is also used) is written in binary format.  A binary dump file
will be about the same size as a text version, but will typically
write out much faster.


----------


Restrictions
""""""""""""


The *vtk* style does not support writing of gzipped dump files.

The *vtk* dump style is part of the USER-VTK package. It is only
enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

To use this dump style, you also must link to the VTK library.  See
the info in lib/vtk/README and insure the Makefile.lammps file in that
directory is appropriate for your machine.

The *vtk* dump style supports neither buffering or custom format
strings.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump image <dump_image>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`

Default
"""""""

By default, files are written in ASCII format. If the file extension
is not one of .vtk, .vtp or .vtu, the legacy VTK file format is used.

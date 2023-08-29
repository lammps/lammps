.. index:: dump molfile

dump molfile command
====================

Syntax
""""""

.. code-block:: LAMMPS

   dump ID group-ID molfile N file format path

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be imaged
* molfile = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write to
* format = file format to be used
* path = file path with plugins (optional)

Examples
""""""""

.. code-block:: LAMMPS

   dump mf1 all molfile 10 melt1.xml hoomd
   dump mf2 all molfile 10 melt2-*.pdb pdb .
   dump mf3 all molfile 50 melt3.xyz xyz .:/home/akohlmey/vmd/plugins/LINUX/molfile

Description
"""""""""""

Dump a snapshot of atom coordinates and selected additional quantities
to one or more files every N timesteps in one of several formats.
Only information for atoms in the specified group is dumped.  This
specific dump style uses molfile plugins that are bundled with the
`VMD <https://www.ks.uiuc.edu/Research/vmd>`_ molecular visualization and
analysis program.

Unless the filename contains a \* character, the output will be written
to one single file with the specified format. Otherwise there will be
one file per snapshot and the \* will be replaced by the time step number
when the snapshot is written.

.. note::

   Because periodic boundary conditions are enforced only on
   timesteps when neighbor lists are rebuilt, the coordinates of an atom
   written to a dump file may be slightly outside the simulation box.

The molfile plugin API has a few restrictions that have to be honored
by this dump style: the number of atoms must not change, the atoms
must be sorted, outside of the coordinates no change in atom properties
(like type, mass, charge) will be recorded.

----------

The *format* keyword determines what format is used to write out the
dump. For this to work, LAMMPS must be able to find and load a
compatible molfile plugin that supports this format.  Settings made via
the :doc:`dump_modify <dump_modify>` command can alter per atom properties
like element names.

The *path* keyword determines which in directories. This is a "path"
like other search paths, i.e. it can contain multiple directories
separated by a colon (or semicolon on Windows). This keyword is
optional and default to ".", the current directory.

The *unwrap* option of the :doc:`dump_modify <dump_modify>` command allows
coordinates to be written "unwrapped" by the image flags for each atom.
Unwrapped means that if the atom has passed through a periodic boundary
one or more times, the value is printed for what the coordinate would be
if it had not been wrapped back into the periodic box.  Note that these
coordinates may thus be far outside the box size stored with the
snapshot.

----------

Dumps are performed on timesteps that are a multiple of N (including
timestep 0) and on the last timestep of a minimization if the
minimization converges.  Note that this means a dump will not be
performed on the initial timestep after the dump command is invoked,
if the current timestep is not a multiple of N.  This behavior can be
changed via the :doc:`dump_modify first <dump_modify>` command, which can
be useful if the dump command is invoked after a minimization ended on
an arbitrary timestep.  N can be changed between runs by using the
:doc:`dump_modify every <dump_modify>` command. The :doc:`dump_modify every <dump_modify>` command also allows a variable to be used to
determine the sequence of timesteps on which dump files are written.

----------

Restrictions
""""""""""""

The *molfile* dump style is part of the MOLFILE package.  It is
only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Molfile plugins provide a consistent programming interface to read and
write file formats commonly used in molecular simulations. The
MOLFILE package only provides the interface code, not the plugins.
These can be obtained from a VMD installation which has to match the
platform that you are using to compile LAMMPS for. By adding plugins
to VMD, support for new file formats can be added to LAMMPS (or VMD
or other programs that use them) without having to re-compile the
application itself.  The plugins are installed in the directory:
<VMDHOME>/plugins/<VMDARCH>/molfile

.. note::

   while the programming interface (API) to the plugins is backward
   compatible, the binary interface (ABI) has been changing over time, so
   it is necessary to compile this package with the plugin header files
   from VMD that match the binary plugins.  These header files in the
   directory: <VMDHOME>/plugins/include For convenience, the package ships
   with a set of header files that are compatible with VMD 1.9 and 1.9.1
   (June 2012)

----------

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`

Default
"""""""

The default path is ".". All other properties have to be specified.

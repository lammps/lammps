.. index:: write\_dump

write\_dump command
===================

Syntax
""""""


.. parsed-literal::

   write_dump group-ID style file dump-args modify dump_modify-args

* group-ID = ID of the group of atoms to be dumped
* style = any of the supported :doc:`dump styles <dump>`
* file = name of file to write dump info to
* dump-args = any additional args needed for a particular :doc:`dump style <dump>`
* modify = all args after this keyword are passed to :doc:`dump_modify <dump_modify>` (optional)
* dump-modify-args = args for :doc:`dump_modify <dump_modify>` (optional)


Examples
""""""""


.. parsed-literal::

   write_dump all atom dump.atom
   write_dump subgroup atom dump.run.bin
   write_dump all custom dump.myforce.\* id type x y vx fx
   write_dump flow custom dump.%.myforce id type c_myF[3] v_ke modify sort id
   write_dump all xyz system.xyz modify sort id element O H
   write_dump all image snap\*.jpg type type size 960 960 modify backcolor white
   write_dump all image snap\*.jpg element element &
      bond atom 0.3 shiny 0.1 ssao yes 6345 0.2 size 1600 1600  &
      modify backcolor white element C C O H N C C C O H H S O H

Description
"""""""""""

Dump a single snapshot of atom quantities to one or more files for the
current state of the system.  This is a one-time immediate operation,
in contrast to the :doc:`dump <dump>` command which will will set up a
dump style to write out snapshots periodically during a running
simulation.

The syntax for this command is mostly identical to that of the
:doc:`dump <dump>` and :doc:`dump_modify <dump_modify>` commands as if
they were concatenated together, with the following exceptions: There
is no need for a dump ID or dump frequency and the keyword *modify* is
added.  The latter is so that the full range of
:doc:`dump_modify <dump_modify>` options can be specified for the single
snapshot, just as they can be for multiple snapshots.  The *modify*
keyword separates the arguments that would normally be passed to the
*dump* command from those that would be given the *dump\_modify*.  Both
support optional arguments and thus LAMMPS needs to be able to cleanly
separate the two sets of args.

Note that if the specified filename uses wildcard characters "\*" or
"%", as supported by the :doc:`dump <dump>` command, they will operate
in the same fashion to create the new filename(s).  Normally, :doc:`dump image <dump_image>` files require a filename with a "\*" character
for the timestep.  That is not the case for the write\_dump command; no
wildcard "\*" character is necessary.


----------


Restrictions
""""""""""""


All restrictions for the :doc:`dump <dump>` and
:doc:`dump_modify <dump_modify>` commands apply to this command as well,
with the exception of the :doc:`dump image <dump_image>` filename not
requiring a wildcard "\*" character, as noted above.

Since dumps are normally written during a :doc:`run <run>` or :doc:`energy minimization <minimize>`, the simulation has to be ready to run
before this command can be used.  Similarly, if the dump requires
information from a compute, fix, or variable, the information needs to
have been calculated for the current timestep (e.g. by a prior run),
else LAMMPS will generate an error message.

For example, it is not possible to dump per-atom energy with this
command before a run has been performed, since no energies and forces
have yet been calculated.  See the :doc:`variable <variable>` doc page
section on Variable Accuracy for more information on this topic.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump image <dump_image>`,
:doc:`dump_modify <dump_modify>`

Default
"""""""

The defaults are listed on the doc pages for the :doc:`dump <dump>` and
:doc:`dump image <dump_image>` and :doc:`dump_modify <dump_modify>`
commands.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

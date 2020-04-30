.. index:: fix nve/noforce

fix nve/noforce command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nve

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/noforce = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 wall nve/noforce

Description
"""""""""""

Perform updates of position, but not velocity for atoms in the group
each timestep.  In other words, the force on the atoms is ignored and
their velocity is not updated.  The atom velocities are used to update
their positions.

This can be useful for wall atoms, when you set their velocities, and
want the wall to move (or stay stationary) in a prescribed fashion.

This can also be accomplished via the :doc:`fix setforce <fix_setforce>`
command, but with fix nve/noforce, the forces on the wall atoms are
unchanged, and can thus be printed by the :doc:`dump <dump>` command or
queried with an equal-style :doc:`variable <variable>` that uses the
fcm() group function to compute the total force on the group of atoms.

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`

**Default:** none

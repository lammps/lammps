.. index:: fix enforce2d
.. index:: fix enforce2d/kk

fix enforce2d command
=====================

Accelerator Variants: *enforce2d/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID enforce2d

* ID, group-ID are documented in :doc:`fix <fix>` command
* enforce2d = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 5 all enforce2d

Description
"""""""""""

Zero out the z-dimension velocity and force on each atom in the group.
This is useful when running a 2d simulation to ensure that atoms do
not move from their initial z coordinate.

.. note::

   Fixes are applied in the order in which they are defined.  If another
   fix that adds forces or torques (e.g. a wall fix or :doc:`fix addforce
   <fix_addforce>`) is defined *after* fix enforce2d, it may add
   out-of-plane force components again after they were set to zero.
   Thus fix enforce2d should be defined after all such fixes.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

none


Default
"""""""

none

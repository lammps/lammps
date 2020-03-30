.. index:: fix store/force

fix store/force command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID store/force

* ID, group-ID are documented in :doc:`fix <fix>` command
* store/force = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all store/force

Description
"""""""""""

Store the forces on atoms in the group at the point during each
timestep when the fix is invoked, as described below.  This is useful
for storing forces before constraints or other boundary conditions are
computed which modify the forces, so that unmodified forces can be
:doc:`written to a dump file <dump>` or accessed by other :doc:`output commands <Howto_output>` that use per-atom quantities.

This fix is invoked at the point in the velocity-Verlet timestepping
immediately after :doc:`pair <pair_style>`, :doc:`bond <bond_style>`,
:doc:`angle <angle_style>`, :doc:`dihedral <dihedral_style>`,
:doc:`improper <improper_style>`, and :doc:`long-range <kspace_style>`
forces have been calculated.  It is the point in the timestep when
various fixes that compute constraint forces are calculated and
potentially modify the force on each atom.  Examples of such fixes are
:doc:`fix shake <fix_shake>`, :doc:`fix wall <fix_wall>`, and :doc:`fix indent <fix_indent>`.

.. note::

   The order in which various fixes are applied which operate at
   the same point during the timestep, is the same as the order they are
   specified in the input script.  Thus normally, if you want to store
   per-atom forces due to force field interactions, before constraints
   are applied, you should list this fix first within that set of fixes,
   i.e. before other fixes that apply constraints.  However, if you wish
   to include certain constraints (e.g. fix shake) in the stored force,
   then it could be specified after some fixes and before others.

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix produces a per-atom array which can be accessed by various
:doc:`output commands <Howto_output>`.  The number of columns for each
atom is 3, and the columns store the x,y,z forces on each atom.  The
per-atom values be accessed on any timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix store_state <fix_store_state>`

**Default:** none

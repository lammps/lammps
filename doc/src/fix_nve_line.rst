.. index:: fix nve/line

fix nve/line command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nve/line

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/line = style name of this fix command

Examples
""""""""


.. parsed-literal::

   fix 1 all nve/line

Description
"""""""""""

Perform constant NVE integration to update position, velocity,
orientation, and angular velocity for line segment particles in the
group each timestep.  V is volume; E is energy.  This creates a system
trajectory consistent with the microcanonical ensemble.  See :doc:`Howto spherical <Howto_spherical>` doc page for an overview of using line
segment particles.

This fix differs from the :doc:`fix nve <fix_nve>` command, which
assumes point particles and only updates their position and velocity.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the ASPHERE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix requires that particles be line segments as defined by the
:doc:`atom_style line <atom_style>` command.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nve/asphere <fix_nve_asphere>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

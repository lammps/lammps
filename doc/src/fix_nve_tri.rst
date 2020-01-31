.. index:: fix nve/tri

fix nve/tri command
===================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nve/tri

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/tri = style name of this fix command

Examples
""""""""


.. parsed-literal::

   fix 1 all nve/tri

Description
"""""""""""

Perform constant NVE integration to update position, velocity,
orientation, and angular momentum for triangular particles in the
group each timestep.  V is volume; E is energy.  This creates a system
trajectory consistent with the microcanonical ensemble.  See the
:doc:`Howto spherical <Howto_spherical>` doc page for an overview of
using triangular particles.

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

This fix requires that particles be triangles as defined by the
:doc:`atom_style tri <atom_style>` command.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nve/asphere <fix_nve_asphere>`

**Default:** none



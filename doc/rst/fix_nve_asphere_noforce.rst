.. index:: fix nve/asphere/noforce

fix nve/asphere/noforce command
===============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nve/asphere/noforce

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/asphere/noforce = style name of this fix command

Examples
""""""""

fix 1 all nve/asphere/noforce

Description
"""""""""""

Perform updates of position and orientation, but not velocity or
angular momentum for atoms in the group each timestep.  In other
words, the force and torque on the atoms is ignored and their velocity
and angular momentum are not updated.  The atom velocities and
angular momenta are used to update their positions and orientation.

This is useful as an implicit time integrator for Fast Lubrication
Dynamics, since the velocity and angular momentum are updated by the
:doc:`pair\_style lubricuteU <pair_lubricateU>` command.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the ASPHERE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the :doc:`atom\_style ellipsoid <atom_style>`
command.

All particles in the group must be finite-size.  They cannot be point
particles, but they can be aspherical or spherical as defined by their
shape attribute.

Related commands
""""""""""""""""

:doc:`fix nve/noforce <fix_nve_noforce>`, :doc:`fix nve/asphere <fix_nve_asphere>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

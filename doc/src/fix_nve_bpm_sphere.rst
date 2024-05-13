.. index:: fix nve/bpm/sphere

fix nve/bpm/sphere command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/bpm/sphere

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/bpm/sphere = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *disc*

  .. parsed-literal::

       *disc* value = none = treat particles as 2d discs, not spheres

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/bpm/sphere
   fix 1 all nve/bpm/sphere disc

Description
"""""""""""

.. versionadded:: 4May2022

Perform constant NVE integration to update position, velocity, angular
velocity, and quaternion orientation for finite-size spherical
particles in the group each timestep.  V is volume; E is energy.  This
creates a system trajectory consistent with the microcanonical
ensemble.

This fix differs from the :doc:`fix nve <fix_nve>` command, which
assumes point particles and only updates their position and velocity.
It also differs from the :doc:`fix nve/sphere <fix_nve_sphere>`
command which assumes finite-size spheroid particles which do not
store a quaternion.  It thus does not update a particle's orientation
or quaternion.

If the *disc* keyword is used, then each particle is treated as a 2d
disc (circle) instead of as a sphere.  This is only possible for 2d
simulations, as defined by the :doc:`dimension <dimension>` keyword.
The only difference between discs and spheres in this context is their
moment of inertia, as used in the time integration.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BPM package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

This fix requires that atoms store torque, angular velocity (omega), a
radius, and a quaternion as defined by the :doc:`atom_style bpm/sphere
<atom_style>` command.

All particles in the group must be finite-size spheres with
quaternions.  They cannot be point particles.

Use of the *disc* keyword is only allowed for 2d simulations, as
defined by the :doc:`dimension <dimension>` keyword.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nve/sphere <fix_nve_sphere>`

Default
"""""""

none


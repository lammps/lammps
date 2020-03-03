.. index:: fix smd/integrate_ulsph

fix smd/integrate_ulsph command
================================

Syntax
""""""


.. code-block:: LAMMPS

   fix ID group-ID smd/integrate_ulsph keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/integrate\_ulsph = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = adjust\_radius or limit\_velocity

adjust\_radius values = adjust\_radius\_factor min\_nn max\_nn
      adjust\_radius\_factor = factor which scale the smooth/kernel radius
      min\_nn = minimum number of neighbors
      max\_nn = maximum number of neighbors
limit\_velocity values = max\_velocity
      max\_velocity = maximum allowed velocity.

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all smd/integrate_ulsph adjust_radius 1.02 25 50
   fix 1 all smd/integrate_ulsph limit_velocity 1000

Description
"""""""""""

The fix performs explicit time integration for particles which
interact with the updated Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth Mach
Dynamics in LAMMPS.

The *adjust\_radius* keyword activates dynamic adjustment of the
per-particle SPH smoothing kernel radius such that the number of
neighbors per particles remains within the interval *min\_nn* to
*max\_nn*. The parameter *adjust\_radius\_factor* determines the amount
of adjustment per timestep. Typical values are *adjust\_radius\_factor*
=1.02, *min\_nn* =15, and *max\_nn* =20.

The *limit\_velocity* keyword will control the velocity, scaling the norm of
the velocity vector to max\_vel in case it exceeds this velocity limit.

**Restart, fix\_modify, output, run start/stop, minimize info:**

Currently, no part of USER-SMD supports restarting nor
minimization. This fix has no outputs.

Restrictions
""""""""""""


This fix is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

**Related commands:** none

**Default:** none

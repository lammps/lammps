.. index:: fix smd/integrate_ulsph

fix smd/integrate_ulsph command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID smd/integrate_ulsph keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/integrate_ulsph = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = adjust_radius or limit_velocity

adjust_radius values = adjust_radius_factor min_nn max_nn
      adjust_radius_factor = factor which scale the smooth/kernel radius
      min_nn = minimum number of neighbors
      max_nn = maximum number of neighbors
limit_velocity values = max_velocity
      max_velocity = maximum allowed velocity.

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all smd/integrate_ulsph adjust_radius 1.02 25 50
   fix 1 all smd/integrate_ulsph limit_velocity 1000

Description
"""""""""""

The fix performs explicit time integration for particles which
interact with the updated Lagrangian SPH pair style.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to using Smooth Mach
Dynamics in LAMMPS.

The *adjust_radius* keyword activates dynamic adjustment of the
per-particle SPH smoothing kernel radius such that the number of
neighbors per particles remains within the interval *min_nn* to
*max_nn*. The parameter *adjust_radius_factor* determines the amount
of adjustment per timestep. Typical values are *adjust_radius_factor*
=1.02, *min_nn* =15, and *max_nn* =20.

The *limit_velocity* keyword will control the velocity, scaling the norm of
the velocity vector to max_vel in case it exceeds this velocity limit.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Currently, no part of MACHDYN supports restarting nor
minimization. This fix has no outputs.

Restrictions
""""""""""""

This fix is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

none


Default
"""""""

none

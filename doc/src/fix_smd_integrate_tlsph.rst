.. index:: fix smd/integrate\_tlsph

fix smd/integrate\_tlsph command
================================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID smd/integrate_tlsph keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/integrate\_tlsph = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *limit\_velocity*


.. parsed-literal::

     *limit_velocity* value = max_vel
       max_vel = maximum allowed velocity

Examples
""""""""


.. parsed-literal::

   fix 1 all smd/integrate_tlsph
   fix 1 all smd/integrate_tlsph limit_velocity 1000

Description
"""""""""""

The fix performs explicit time integration for particles which
interact according with the Total-Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth Mach
Dynamics in LAMMPS.

The *limit\_velocity* keyword will control the velocity, scaling the
norm of the velocity vector to max\_vel in case it exceeds this
velocity limit.

**Restart, fix\_modify, output, run start/stop, minimize info:**

Currently, no part of USER-SMD supports restarting nor
minimization. This fix has no outputs.

Restrictions
""""""""""""


This fix is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`smd/integrate\_ulsph <fix_smd_integrate_ulsph>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

.. index:: fix smd/adjust_dt

fix smd/adjust_dt command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID smd/adjust_dt arg

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/adjust\_dt = style name of this fix command
* arg = *s\_fact*

  .. parsed-literal::

       *s_fact* = safety factor

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all smd/adjust_dt 0.1

Description
"""""""""""

The fix calculates a new stable time increment for use with the SMD
time integrators.

The stable time increment is based on multiple conditions. For the SPH
pair styles, a CFL criterion (Courant, Friedrichs & Lewy, 1928) is
evaluated, which determines the speed of sound cannot propagate
further than a typical spacing between particles within a single time
step to ensure no information is lost. For the contact pair styles, a
linear analysis of the pair potential determines a stable maximum time
step.

This fix inquires the minimum stable time increment across all
particles contained in the group for which this fix is defined. An
additional safety factor *s\_fact* is applied to the time increment.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to use Smooth Mach
Dynamics in LAMMPS.

**Restart, fix\_modify, output, run start/stop, minimize info:**

Currently, no part of USER-SMD supports restarting nor minimization.

Restrictions
""""""""""""

This fix is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`smd/tlsph\_dt <compute_smd_tlsph_dt>`

**Default:** none

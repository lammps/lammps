.. index:: compute smd/tlsph/dt

compute smd/tlsph/dt command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/tlsph/dt

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/tlsph/dt = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/tlsph/dt

Description
"""""""""""

Define a computation that outputs the CFL-stable time increment per
particle.  This time increment is essentially given by the speed of
sound, divided by the SPH smoothing length.  Because both the speed of
sound and the smoothing length typically change during the course of a
simulation, the stable time increment needs to be re-computed every
time step.  This calculation is performed automatically in the
relevant SPH pair styles and this compute only serves to make the
stable time increment accessible for output purposes.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-particle values will be given in :doc:`units <units>` of time.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This compute can only be used for particles interacting with the
Total-Lagrangian SPH pair style.

Related commands
""""""""""""""""

:doc:`smd/adjust/dt <fix_smd_adjust_dt>`

Default
"""""""

none

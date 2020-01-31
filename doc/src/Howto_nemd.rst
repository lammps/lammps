NEMD simulations
================

Non-equilibrium molecular dynamics or NEMD simulations are typically
used to measure a fluid's rheological properties such as viscosity.
In LAMMPS, such simulations can be performed by first setting up a
non-orthogonal simulation box (see the preceding Howto section).

A shear strain can be applied to the simulation box at a desired
strain rate by using the :doc:`fix deform <fix_deform>` command.  The
:doc:`fix nvt/sllod <fix_nvt_sllod>` command can be used to thermostat
the sheared fluid and integrate the SLLOD equations of motion for the
system.  Fix nvt/sllod uses :doc:`compute temp/deform <compute_temp_deform>` to compute a thermal temperature
by subtracting out the streaming velocity of the shearing atoms.  The
velocity profile or other properties of the fluid can be monitored via
the :doc:`fix ave/chunk <fix_ave_chunk>` command.

.. note::

   A recent (2017) book by :ref:`(Daivis and Todd) <Daivis-nemd>`
   discusses use of the SLLOD method and non-equilibrium MD (NEMD)
   thermostatting generally, for both simple and complex fluids,
   e.g. molecular systems.  The latter can be tricky to do correctly.

As discussed in the previous section on non-orthogonal simulation
boxes, the amount of tilt or skew that can be applied is limited by
LAMMPS for computational efficiency to be 1/2 of the parallel box
length.  However, :doc:`fix deform <fix_deform>` can continuously strain
a box by an arbitrary amount.  As discussed in the :doc:`fix deform <fix_deform>` command, when the tilt value reaches a limit,
the box is flipped to the opposite limit which is an equivalent tiling
of periodic space.  The strain rate can then continue to change as
before.  In a long NEMD simulation these box re-shaping events may
occur many times.

In a NEMD simulation, the "remap" option of :doc:`fix deform <fix_deform>` should be set to "remap v", since that is what
:doc:`fix nvt/sllod <fix_nvt_sllod>` assumes to generate a velocity
profile consistent with the applied shear strain rate.

An alternative method for calculating viscosities is provided via the
:doc:`fix viscosity <fix_viscosity>` command.

NEMD simulations can also be used to measure transport properties of a fluid
through a pore or channel. Simulations of steady-state flow can be performed
using the :doc:`fix flow/gauss <fix_flow_gauss>` command.


----------


.. _Daivis-nemd:



**(Daivis and Todd)** Daivis and Todd, Nonequilibrium Molecular Dynamics (book),
Cambridge University Press, https://doi.org/10.1017/9781139017848, (2017).



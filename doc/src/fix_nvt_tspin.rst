.. index:: fix nvt/tspin

fix nvt/tspin command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nvt/tspin temp Tstart Tstop Tdamp keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nvt/tspin = style name of this fix command
* *temp* values = Tstart Tstop Tdamp

  .. parsed-literal::

       Tstart, Tstop = desired spin temperature at start and end of the run
                       (temperature units, K in metal units)
       Tdamp = spin temperature damping parameter (time units)

* zero or more additional keyword/value pairs may be appended
* keyword = *tchain* or *tloop* or *drag* or *lattice* or *spin* or *spinmass*

  .. parsed-literal::

       *tchain* value = N
         N = length of the Nose-Hoover chain coupled to the spins
       *tloop* value = M
         M = number of sub-cycles of the chain integration per timestep
       *drag* value = Df
         Df = drag factor added to the chain (0 = no drag)
       *lattice*, *spin*, *spinmass* = as for :doc:`fix nve/tspin <fix_nve_tspin>`

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt/tspin temp 300.0 300.0 0.05 spinmass 0.0075
   fix 1 all nvt/tspin temp 10.0 600.0 0.1 lattice frozen tchain 5

Description
"""""""""""

.. versionadded:: TBD

Perform the same time integration of inertial spin dynamics as
:doc:`fix nve/tspin <fix_nve_tspin>`, and in addition couple the spin
velocities to a Nose-Hoover chain thermostat held at the requested spin
temperature.

The spin temperature is defined from the kinetic energy stored in the
spin degrees of freedom,

.. math::

   \frac{3}{2} N k_B T_s = \sum_i \frac{1}{2} m_s \left| \vec{v}^{s}_i \right|^2

where the sum runs over the magnetic atoms of the group and N is their
number.  The chain is integrated with the same scheme as :doc:`fix nvt
<fix_nh>`, so the *tchain*, *tloop* and *drag* keywords have the same
meaning as there.  Tstart and Tstop let the target spin temperature be
ramped linearly over the course of the run.

This fix thermostats the spin degrees of freedom only.  The atomic
positions and velocities are integrated without a thermostat when
*lattice* is *moving*; combine this fix with a separate :doc:`fix nvt
<fix_nh>` or :doc:`fix langevin <fix_langevin>` using *lattice frozen*
here if the lattice also has to be thermostatted.

The warning in :doc:`fix nve/tspin <fix_nve_tspin>` about the need for a
longitudinal potential on the spin modulus applies to this fix as well.

An alternative kinetic thermostat for the same degrees of freedom is
:doc:`fix langevin/tspin <fix_langevin_tspin>`.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the Nose-Hoover chain to :doc:`binary
restart files <restart>`, so that a simulation can continue correctly.
The *tchain* value must be the same in the restarted run.

This fix computes a global scalar, the cumulative energy change caused
by the thermostat, which can be accessed by various :doc:`output
commands <Howto_output>`.  Adding it to the potential energy and to the
spin kinetic energy from :doc:`compute ke/tspin <compute_ke_tspin>`
gives a quantity that should be conserved.  The scalar value is
"extensive".

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The *nvt/tspin* fix is part of the SPIN package.  This style is only
enabled if LAMMPS was built with this package.  See the :doc:`Build
package <Build_package>` page for more info.

This fix requires :doc:`atom_style tspin <atom_style>`, and cannot be
used with *spin frozen*.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`,
:doc:`compute ke/tspin <compute_ke_tspin>`,
:doc:`velocity/tspin <velocity_tspin>`,
:doc:`fix nvt <fix_nh>`

Default
""""""""

The default keyword values are tchain = 3, tloop = 1 and drag = 0.0.

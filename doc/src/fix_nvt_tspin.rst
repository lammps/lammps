.. index:: fix nvt/tspin
.. index:: fix npt/tspin
.. index:: fix nph/tspin

fix nvt/tspin command
=====================

fix npt/tspin command
=====================

fix nph/tspin command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style_name keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *nvt/tspin* or *npt/tspin* or *nph/tspin*
* these fixes accept the corresponding keywords of
  :doc:`fix nvt, npt and nph <fix_nh>`
* additional keyword = *spinmass* or *lattice* or *spin*

  .. parsed-literal::

       *spinmass* value = ms
         ms = spin mass multiplier relative to atomic mass (adim)
       *lattice* value = *moving* or *frozen*
         moving = integrate the atomic positions and velocities
         frozen = hold the atoms on a fixed lattice, spins only
       *spin* value = *moving* or *frozen*
         moving = integrate the spins
         frozen = hold the spins fixed, lattice only

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt/tspin temp 300.0 300.0 0.05 spinmass 0.0075
   fix 1 all nvt/tspin temp 10.0 600.0 0.1 lattice frozen tchain 5 spinmass 0.0075
   fix 1 all npt/tspin temp 300.0 300.0 0.1 iso 0.0 0.0 1.0 spinmass 0.0075
   fix 1 all nph/tspin iso 0.0 0.0 1.0 spinmass 0.0075

Description
"""""""""""

.. versionadded:: TBD

These fixes extend :doc:`fix nvt, npt and nph <fix_nh>` to inertial
spin dynamics.  As described for :doc:`fix nve/tspin <fix_nve_tspin>`,
the auxiliary spin dynamics are intended for equilibrium sampling and
do not represent physical spin precession.

The *nvt/tspin* and *npt/tspin* styles use a separate Nose-Hoover chain
to thermostat the auxiliary spin velocities.  The spin and lattice
chains share the target temperature, damping time, *tchain*, *tloop*,
and *drag* settings.  The *nph/tspin* style applies a barostat without
thermostatting either the lattice or the spins.

The spin thermostat uses the kinetic temperature defined by
:doc:`compute temp/tspin <compute_temp_tspin>`.  Spin masses can be
assigned with *spinmass* or :doc:`velocity/tspin <velocity_tspin>`.

The *lattice frozen* option allows spin-only sampling with *nvt/tspin*;
it cannot be combined with a barostat.  The longitudinal-potential and
magnetic-force requirements described for :doc:`fix nve/tspin
<fix_nve_tspin>` also apply here.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These fixes write their thermostat and barostat states to :doc:`binary
restart files <restart>`.  The *tchain* value must remain unchanged
when restarting.

These fixes compute a global scalar containing the thermostat and
barostat energy contribution, including the spin thermostat chain.
The scalar is "extensive".

The auxiliary spin kinetic energy is not included in *econserve*.
Add the energy from :doc:`compute ke/tspin <compute_ke_tspin>` to
monitor the extended conserved energy for constant temperature and
pressure targets, zero *drag*, and no other energy exchange.  All
relevant potential-energy contributions must be included in *pe*.

Standard :doc:`fix_modify <fix_modify>` options and other output
conventions are described in :doc:`fix nvt, npt and nph <fix_nh>`.
These fixes are not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

These fixes require :doc:`atom_style spin <atom_style>` and the SPIN
package.  :doc:`Dynamic groups <group>` are not supported, and at least
one of *lattice* and *spin* must be *moving*.  A barostat cannot be
used with *lattice frozen*.

There is no KOKKOS version.  Magnetic interactions must satisfy the
full-force requirement described for :doc:`fix nve/tspin <fix_nve_tspin>`.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`,
:doc:`compute ke/tspin <compute_ke_tspin>`,
:doc:`compute temp/tspin <compute_temp_tspin>`,
:doc:`velocity/tspin <velocity_tspin>`,
:doc:`fix nvt <fix_nh>`

Default
""""""""

The keyword defaults are the same as for :doc:`fix nvt, fix npt and fix nph
<fix_nh>`, plus lattice = moving and spin = moving.

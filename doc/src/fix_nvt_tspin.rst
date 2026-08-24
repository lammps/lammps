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
* these fixes take the same keywords as :doc:`fix nvt, fix npt and fix nph
  <fix_nh>`, in particular *temp*, *iso*, *aniso*, *tchain*, *pchain*,
  *tloop*, *drag* and *couple*
* additional keyword = *spinmass* or *lattice* or *spin*

  .. parsed-literal::

       *spinmass* value = ms
         ms = spin mass of each atom, as a multiple of its atomic mass (adim)
       *lattice* value = *moving* or *frozen*
         moving = integrate the atomic positions and velocities
         frozen = hold the atoms on a fixed lattice, spins only
       *spin* value = *moving*
         the spins are always integrated by these fixes

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt/tspin temp 300.0 300.0 0.05 spinmass 0.0075
   fix 1 all nvt/tspin temp 10.0 600.0 0.1 lattice frozen tchain 5
   fix 1 all npt/tspin temp 300.0 300.0 0.1 iso 0.0 0.0 1.0 spinmass 0.0075
   fix 1 all nph/tspin iso 0.0 0.0 1.0 spinmass 0.0075

Description
"""""""""""

.. versionadded:: TBD

These fixes are the inertial spin dynamics counterparts of :doc:`fix nvt,
fix npt and fix nph <fix_nh>`.  They integrate the lattice exactly as those
fixes do, integrate the spins as :doc:`fix nve/tspin <fix_nve_tspin>` does,
and in addition couple the spin velocities to a second Nose-Hoover chain.

That spin chain is independent of the lattice chain, but is driven by the
same target temperature and the same *tchain*, *tloop* and *drag* settings,
so the spin-lattice system is thermostatted as a whole.  *nph/tspin*
controls the pressure only and applies no thermostat at all, to either
subsystem.

The spin temperature is defined from the kinetic energy stored in the
spin degrees of freedom,

.. math::

   \frac{3}{2} N k_B T_s = \sum_i \frac{1}{2} m_s \left| \vec{v}^{s}_i \right|^2

where the sum runs over the magnetic atoms of the group and N is their
number.  The *spinmass* keyword sets the spin mass of every atom in the group as a
multiple of the atomic mass of its type; it can equivalently be set with
the :doc:`velocity/tspin <velocity_tspin>` command.

With *lattice frozen* the atoms are held fixed and only the spins are
integrated and thermostatted, which is the usual setup for pure spin
dynamics.  A barostat cannot be combined with *lattice frozen*.

The warning in :doc:`fix nve/tspin <fix_nve_tspin>` about the need for a
longitudinal potential on the spin modulus applies to this fix as well.

An alternative kinetic thermostat for the same degrees of freedom is
:doc:`fix langevin/tspin <fix_langevin_tspin>`.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These fixes write the state of both Nose-Hoover chains to :doc:`binary
restart files <restart>`, so that a simulation can continue correctly.
The *tchain* value must be the same in the restarted run.

These fixes compute a global scalar, the cumulative energy change caused
by the thermostat and barostat, which can be accessed by various
:doc:`output commands <Howto_output>`.  Adding it to the potential energy
and to the spin kinetic energy from :doc:`compute ke/tspin
<compute_ke_tspin>` gives a quantity that should be conserved.  The scalar
value is "extensive".

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

These fixes are part of the SPIN package.  They are only enabled if LAMMPS
was built with this package.  See the :doc:`Build package <Build_package>`
page for more info.

They require :doc:`atom_style tspin <atom_style>`.  A barostat cannot be
combined with *lattice frozen*.  The limitation of :doc:`fix nve/tspin
<fix_nve_tspin>` concerning the *zeeman* keyword of :doc:`fix
precession/spin <fix_precession_spin>` applies here as well.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`,
:doc:`compute ke/tspin <compute_ke_tspin>`,
:doc:`velocity/tspin <velocity_tspin>`,
:doc:`fix nvt <fix_nh>`

Default
""""""""

The keyword defaults are the same as for :doc:`fix nvt, fix npt and fix nph
<fix_nh>`, plus lattice = moving.

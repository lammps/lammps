.. index:: fix nve/tspin

fix nve/tspin command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/tspin keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/tspin = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *lattice* or *spin* or *spinmass*

  .. parsed-literal::

       *lattice* value = *moving* or *frozen*
         moving = integrate the atomic positions and velocities
         frozen = hold the atoms on a fixed lattice
       *spin* value = *moving* or *frozen*
         moving = integrate the spins
         frozen = hold the spins fixed
       *spinmass* value = ms
         ms = spin mass of each atom, as a multiple of its atomic mass (adim)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/tspin spinmass 0.0075
   fix 1 all nve/tspin lattice frozen spinmass 0.01
   fix 1 all nve/tspin lattice moving spin moving

Description
"""""""""""

.. versionadded:: TBD

Perform a constant NVE integration of inertial spin dynamics, as
described in :ref:`(Tspin) <Tspin1>`, in which the modulus of each
atomic spin is a dynamical degree of freedom rather than a constant of
the motion.

Each spin is written as the product of a direction and a modulus,

.. math::

   \vec{S}_i = \left| \vec{S}_i \right| \, \hat{s}_i

and the full vector :math:`\vec{S}_i` obeys the Newtonian equation of
motion

.. math::

   m_s \frac{d^2 \vec{S}_i}{dt^2} = \vec{F}^{m}_{i}

where :math:`m_s` is a spin mass and :math:`\vec{F}_{i}` is assembled from
the magnetic interactions as explained below.  A spin velocity
:math:`\vec{v}^{s}_i = d\vec{S}_i/dt` is stored for every atom, and the
equation of motion is integrated with the velocity-Verlet algorithm.
After each update the direction and the modulus of the spin are
recomputed, so that both evolve in time.

This is a different physical model from :doc:`fix nve/spin
<fix_nve_spin>`, which integrates the fixed-modulus Landau-Lifshitz
precession of a unit vector with a Suzuki-Trotter decomposition.  The
two fixes cannot be combined on the same group of atoms.  Because the
spin dynamics generated here is second order in time, its
characteristic frequencies scale as :math:`1/\sqrt{m_s}` and are not
the Landau-Lifshitz precession frequencies; the spin mass is a
numerical parameter that controls how fast the spin degrees of freedom
equilibrate, and static thermodynamic averages should be checked to be
independent of it.

How the magnetic force is assembled
"""""""""""""""""""""""""""""""""""""

These fixes run on the ordinary :doc:`atom_style spin <atom_style>`.  The spin
is the per-atom array *sp*, a direction *sp[0..2]* and a modulus *sp[3]* in
Bohr magnetons.  The magnetic interaction continues to use the per-atom *fm*
array in rad.THz, so no new atom style or magnetic-force array is required.

What inertial spin dynamics needs is the derivative of the energy with respect
to the *full* spin vector,

.. math::

   \vec{F}^{m}_{i} = -\frac{\partial U}{\partial \vec{S}_i}
   \qquad \left[ \mathrm{eV}/\mu_B \right]

To keep the existing storage and units, a TSPIN-compatible interaction encodes
this full force in *fm* as

.. math::

   \vec{fm}_i = \frac{\left| \vec{S}_i \right|}{\hbar}\vec{F}^{m}_i.

The fix then recovers the force with

.. math::

   \vec{F}^{m}_{i} = \frac{\hbar}{\left| \vec{S}_i \right|} \, \vec{fm}_i

keeping all three components.  This full-gradient encoding is a stronger
requirement than the torque convention used by fixed-modulus styles.  In the
Landau-Lifshitz equation a component of *fm* parallel to the spin drops out of
:math:`\vec{fm} \times \vec{S}`, so a fixed-modulus interaction need not define
that component as a physical radial force.  TSPIN does use it to drive the
modulus and therefore accepts only interactions that define all three
components through the equation above.  This is the convention used by
*pair_style deepspin* of the DeePMD-kit package, so that pair style needs no
change to its force output.

The spin velocity and the spin mass are kept in the custom per-atom
properties *d2_tspin_vs* and *d_tspin_smass*, created automatically through
:doc:`fix property/atom <fix_property_atom>` by whichever *tspin* style is
defined first, under the reserved fix ID *TSPIN_STATE*.  They migrate with the
atoms and are written to :doc:`binary restart files <restart>` but not to
:doc:`data files <write_data>`, and they can be output with :doc:`compute
property/atom <compute_property_atom>`::

   compute v all property/atom d2_tspin_vs[1] d2_tspin_vs[2] d2_tspin_vs[3]

.. warning::

   *TSPIN_STATE* holds the integrator state.  Do not :doc:`unfix <unfix>` it
   while a *tspin* style is defined; the styles bind to it again at the start
   of every run and stop with an error if it has gone.

The *spinmass* keyword sets the spin mass of every magnetic atom in the group
as a multiple of the atomic mass of its type, :math:`m_s = ms \times m_i`, and
does so again at the start of every run.  If the keyword is not used the fix
never writes the spin masses, and they must come from the
:doc:`velocity/tspin <velocity_tspin>` command or from a restart file; the fix
stops with an error if any magnetic atom in the group still has none.  Give the
keyword in one place only, otherwise the fix silently wins.

An atom with no spin mass and an initial spin modulus smaller than
:math:`10^{-8}` is treated as non-magnetic and is not integrated.  Once an atom
has a positive spin mass it remains a dynamical spin even if its modulus becomes
smaller than this value or passes through zero.  At exactly zero modulus the
*fm* encoding contains no recoverable force because it is proportional to
:math:`\left|\vec{S}_i\right|`; the integrator omits that isolated kick but
continues the Cartesian spin drift, so the zero-modulus state is not an
absorbing state.

The *lattice* keyword selects whether the atomic positions and
velocities are integrated as well (*moving*), or whether the atoms are
held on a fixed lattice (*frozen*) so that only the spins evolve.  The
*spin* keyword correspondingly selects whether the spins are
integrated.  For compatibility with :doc:`fix nve/spin <fix_nve_spin>`,
both keywords also accept the values *yes* and *no*.

.. warning::

   The spin modulus is a free coordinate and needs a potential of its
   own.  A magnetic potential that resolves the magnitude of the local
   moment supplies the longitudinal restoring force itself.  Otherwise
   add :doc:`fix spring/tspin <fix_spring_tspin>`, without which the
   modulus grows without bound.

The kinetic energy carried by the spin degrees of freedom is not part of the
thermodynamic keyword *ke*, and therefore not part of *etotal* or *econserve*
either.  Report it with :doc:`compute ke/tspin <compute_ke_tspin>` and add it
explicitly to get the conserved quantity:

.. code-block:: LAMMPS

   fix_modify    pin energy yes            # modulus spring energy into pe
   compute       ske all ke/tspin
   variable      H equal econserve+c_ske
   thermo_style  custom step temp pe c_ske v_H

.. note::

   It is tempting to fold the spin kinetic energy into a temperature compute
   so that *ke* and *econserve* pick it up on their own.  That does not work:
   :doc:`thermo_modify temp <thermo_modify>` also rewires the thermodynamic
   pressure compute, and :doc:`compute pressure <compute_pressure>` builds its
   scalar from ``dof * kB * T`` of that same compute.  A temperature that
   contains the spin kinetic energy therefore adds a spurious term to the
   scalar pressure while the pressure tensor stays correct.  Adding *c_ske*
   explicitly keeps the two separate.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix itself is written to :doc:`binary restart files
<restart>`.  The spin velocities and spin masses live in the custom per-atom
properties described above, owned by the internal :doc:`fix property/atom
<fix_property_atom>` with ID *TSPIN_STATE*, and that fix writes them to restart
files.  A restarted run continues from the same spin state as long as the
*tspin* style is defined again after :doc:`read_restart <read_restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix cannot be used with a :doc:`dynamic group <group>`.

There is no KOKKOS version.  Running the CPU *tspin* styles inside a KOKKOS
build has not been validated and is not supported.

:doc:`pair_style spin/dipole/cut <pair_spin_dipole>` can be used with this fix:
its energy depends on the complete spin vectors and its *fm* output follows the
full-gradient encoding above.  The other magnetic pair styles of the SPIN
package provide a Landau-Lifshitz torque but do not guarantee that the parallel
component of *fm* represents a radial force, so they stop with an error when a
*tspin* integrator is defined.  A direction-only interaction can be made
compatible by projecting its contribution onto the spin tangent plane before
it is added to *fm*.  The same caveat applies to :doc:`fix precession/spin
<fix_precession_spin>`, which is not checked automatically.

The *nve/tspin* fix is part of the SPIN package.  This style is only
enabled if LAMMPS was built with this package.  See the :doc:`Build
package <Build_package>` page for more info.

This fix requires :doc:`atom_style spin <atom_style>`.

Related commands
""""""""""""""""

:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`,
:doc:`fix spring/tspin <fix_spring_tspin>`,
:doc:`compute ke/tspin <compute_ke_tspin>`,
:doc:`velocity/tspin <velocity_tspin>`,
:doc:`fix nve/spin <fix_nve_spin>`

Default
""""""""

The default keyword values are lattice = moving and spin = moving.

----------

.. _Tspin1:

**(Tspin)** Huang, Bai, Wang, and Xu, Scalable Canonical and
Isothermal-Isobaric Sampling of Coupled Spin-Lattice Systems with
Machine-Learning Potentials, `arXiv:2506.12877
<https://arxiv.org/abs/2506.12877>`_ (2025),
`doi:10.48550/arXiv.2506.12877
<https://doi.org/10.48550/arXiv.2506.12877>`_.

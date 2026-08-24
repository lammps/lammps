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

The SPIN pair styles and :doc:`fix precession/spin <fix_precession_spin>`
do not compute a force.  They accumulate a precession field
:math:`\vec{\omega}_i` (the per-atom array *fm*, an angular frequency in
rad.THz) which drives the fixed-modulus equation
:math:`d\vec{S}_i/dt = \vec{\omega}_i \times \vec{S}_i`.  The energy
associated with such a term depends on the spin *direction* only, so the
corresponding force on the full spin vector is transverse and inversely
proportional to the modulus:

.. math::

   \vec{F}^{\perp}_{i} = \frac{\hbar}{\left| \vec{S}_i \right|}
   \left[ \vec{\omega}_i - \left( \vec{\omega}_i \cdot \hat{s}_i \right)
   \hat{s}_i \right]

Styles whose energy is a genuine function of the full spin vector, such as
:doc:`fix spring/tspin <fix_spring_tspin>` and :doc:`fix langevin/tspin
<fix_langevin_tspin>`, instead accumulate a true force in energy units into
the per-atom array *f_spin*, which is used as it is.  The total force is the
sum of the two contributions.  Treating the precession field directly as a
force, without the projection, would add a spurious radial component and the
total energy would not be conserved.

The *spinmass* keyword sets the spin mass of every atom in the group as
a multiple of the atomic mass of its type, :math:`m_s = ms \times m_i`.
If the keyword is not used, the spin masses must have been assigned
beforehand with the :doc:`velocity/tspin <velocity_tspin>` command; the
fix stops with an error if any magnetic atom in the group has no spin
mass.

The *lattice* keyword selects whether the atomic positions and
velocities are integrated as well (*moving*), or whether the atoms are
held on a fixed lattice (*frozen*) so that only the spins evolve.  The
*spin* keyword correspondingly selects whether the spins are
integrated.  For compatibility with :doc:`fix nve/spin <fix_nve_spin>`,
both keywords also accept the values *yes* and *no*.

.. warning::

   A simulation using inertial spin dynamics must include a
   longitudinal potential acting on the spin modulus, such as
   :doc:`fix spring/tspin <fix_spring_tspin>`.  The SPIN pair styles
   provide a Landau-Lifshitz effective field, whose component parallel
   to :math:`\vec{S}_i` does no work on a fixed-modulus spin but acts
   as a radial driving force once the modulus is free to evolve.
   Without a longitudinal potential the spin modulus is unbounded and
   the simulation diverges.

The kinetic energy carried by the spin degrees of freedom is not
included in the thermodynamic keyword *ke*.  Use :doc:`compute ke/tspin
<compute_ke_tspin>` and reference it with a *c_ID* keyword in
:doc:`thermo_style custom <thermo_style>` to monitor it, and to check
that the total energy is conserved.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  The spin velocities and spin masses are part of
:doc:`atom_style tspin <atom_style>` and are stored in restart files by
the atom style.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The *zeeman* keyword of :doc:`fix precession/spin
<fix_precession_spin>` describes an energy that is linear in the full spin
vector rather than a function of the spin direction alone, but it is
communicated through the precession field like every other term.  Its
contribution to inertial spin dynamics is therefore only approximate and
the total energy is not exactly conserved; the drift is small but does not
vanish as the timestep is reduced.  The SPIN pair styles and the
*anisotropy* keywords of fix precession/spin are not affected and conserve
the energy to the expected second order in the timestep.

The *nve/tspin* fix is part of the SPIN package.  This style is only
enabled if LAMMPS was built with this package.  See the :doc:`Build
package <Build_package>` page for more info.

This fix requires :doc:`atom_style tspin <atom_style>`.

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

**(Tspin)** AUTHORS_TBD, TITLE_TBD, PREPRINT_TBD.

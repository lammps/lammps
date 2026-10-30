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

Perform velocity-Verlet integration of auxiliary inertial spin dynamics
for equilibrium sampling, as described in :ref:`(Tspin) <Tspin1>`.
Both the direction and magnitude of each spin can evolve according to

.. math::

   m_{s,i}\frac{d^2\vec{S}_i}{dt^2}
   = \vec{F}^m_i
   = -\frac{\partial U}{\partial\vec{S}_i},

where :math:`\vec{S}_i` is the full spin vector and :math:`m_{s,i}`
is the auxiliary spin mass.  These dynamics do not describe physical
Landau-Lifshitz precession, as integrated by :doc:`fix nve/spin
<fix_nve_spin>`.

The *lattice* and *spin* keywords select whether atomic positions and
velocities, or spins, are integrated.  Both also accept *yes* and *no*.
The *spinmass* keyword sets each atom's auxiliary spin mass to the
specified multiple of its atomic mass.  Alternatively, spin masses can
be assigned with :doc:`velocity/tspin <velocity_tspin>` or restored from
a restart file.
If specified here, *spinmass* overrides previously assigned values at
the start of each run.

TSPIN requires the full Cartesian magnetic force, including its
longitudinal component.  A compatible interaction must encode

.. math::

   \vec{fm}_i = \frac{|\vec{S}_i|}{\hbar}
   \left(-\frac{\partial U}{\partial\vec{S}_i}\right)

in the SPIN *fm* array.  Compatible variable-moment potentials such as
*pair_style deepspin* of DeePMD-kit satisfy this requirement.
Fixed-modulus SPIN pair styles that do not supply this full gradient
are rejected.

Unless the magnetic potential provides a longitudinal restoring force,
add a potential such as :doc:`fix spring/tspin <fix_spring_tspin>` to
control the spin magnitude.

The auxiliary spin kinetic energy is not included in the
usual thermodynamic quantities *ke*, *etotal*, or *econserve*.
See :doc:`compute temp/tspin <compute_temp_tspin>` for how to
monitor the conserved energy in NVE simulations.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Spin velocities and masses are saved in :doc:`binary restart files
<restart>`, but not in :doc:`data files <write_data>`.  Redefine the
TSPIN integrator after :doc:`read_restart <read_restart>` to continue
from the saved spin state.

No :doc:`fix_modify <fix_modify>` options or :doc:`output quantities
<Howto_output>` are supported.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix requires :doc:`atom_style spin <atom_style>` and the SPIN
package; see :doc:`Build package <Build_package>`.
:doc:`Dynamic groups <group>` and :doc:`run_style respa <run_style>`
are not supported.  There is no KOKKOS version, and use inside a KOKKOS
build is not supported.

Do not combine this fix with :doc:`fix nve/spin <fix_nve_spin>` on the
same atoms.  With SPIN pair styles, this restriction also applies to
disjoint groups.

Magnetic interactions must supply the full Cartesian spin force described
above.  Compatibility of additional magnetic fixes,
such as :doc:`fix precession/spin <fix_precession_spin>`, is not checked
automatically.

Related commands
""""""""""""""""

:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`,
:doc:`fix spring/tspin <fix_spring_tspin>`,
:doc:`compute temp/tspin <compute_temp_tspin>`,
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

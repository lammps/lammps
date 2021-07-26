.. index:: fix dpd/energy
.. index:: fix dpd/energy/kk

fix dpd/energy command
======================

Accelerator Variants: *dpd/energy/kk*

Syntax
""""""

.. parsed-literal::

   fix ID group-ID dpd/energy

* ID, group-ID are documented in :doc:`fix <fix>` command
* dpd/energy = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all dpd/energy

Description
"""""""""""

Perform constant energy dissipative particle dynamics (DPD-E)
integration.  This fix updates the internal energies for particles in
the group at each timestep.  It must be used in conjunction with a
deterministic integrator (e.g. :doc:`fix nve <fix_nve>`) that updates
the particle positions and velocities.

For fix *dpd/energy*, the particle internal temperature is related to
the particle internal energy through a mesoparticle equation of state.
An additional fix must be specified that defines the equation of state
for each particle, e.g. :doc:`fix eos/cv <fix_eos_cv>`.

This fix must be used with the :doc:`pair_style dpd/fdt/energy <pair_style>` command.

Note that numerous variants of DPD can be specified by choosing an
appropriate combination of the integrator and :doc:`pair_style dpd/fdt/energy <pair_style>` command.  DPD under isoenergetic conditions
can be specified by using fix *dpd/energy*, fix *nve* and pair_style
*dpd/fdt/energy*\ .  DPD under isoenthalpic conditions can
be specified by using fix *dpd/energy*, fix *nph* and pair_style
*dpd/fdt/energy*\ .  Examples of each DPD variant are provided in the
examples/PACKAGES/dpd-react directory.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This command is part of the DPD-REACT package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This fix must be used with an additional fix that specifies time
integration, e.g. :doc:`fix nve <fix_nve>`.

The fix *dpd/energy* requires the *dpd* :doc:`atom_style <atom_style>`
to be used in order to properly account for the particle internal
energies and temperature.

The fix *dpd/energy* must be used with an additional fix that specifies the
mesoparticle equation of state for each particle.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>` :doc:`fix eos/cv <fix_eos_cv>`

Default
"""""""

none

----------

.. _Lisal1:

**(Lisal)** M. Lisal, J.K. Brennan, J. Bonet Avalos, "Dissipative
particle dynamics at isothermal, isobaric, isoenergetic, and
isoenthalpic conditions using Shardlow-like splitting algorithms.",
J. Chem. Phys., 135, 204105 (2011).

.. _Larentzos3:

**(Larentzos)** J.P. Larentzos, J.K. Brennan, J.D. Moore, and
W.D. Mattson, "LAMMPS Implementation of Constant Energy Dissipative
Particle Dynamics (DPD-E)", ARL-TR-6863, U.S. Army Research
Laboratory, Aberdeen Proving Ground, MD (2014).

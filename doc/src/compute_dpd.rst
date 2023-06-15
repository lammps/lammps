.. index:: compute dpd

compute dpd command
===================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID dpd

* ID, group-ID are documented in :doc:`compute <compute>` command
* dpd = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all dpd

Description
"""""""""""

Define a computation that accumulates the total internal conductive
energy (:math:`U^{\text{cond}}`), the total internal mechanical energy
(:math:`U^{\text{mech}}`), the total chemical energy (:math:`U^\text{chem}`)
and the *harmonic* average of the internal temperature (:math:`\theta_\text{avg}`)
for the entire system of particles.  See the
:doc:`compute dpd/atom <compute_dpd_atom>` command if you want
per-particle internal energies and internal temperatures.

The system internal properties are computed according to the following
relations:

.. math::

   U^\text{cond} = & \sum_{i=1}^{N} u_{i}^\text{cond} \\
   U^\text{mech} = & \sum_{i=1}^{N} u_{i}^\text{mech} \\
   U^\text{chem} = & \sum_{i=1}^{N} u_{i}^\text{chem} \\
               U = & \sum_{i=1}^{N} (u_{i}^\text{cond}
                     + u_{i}^\text{mech} + u_{i}^\text{chem}) \\
   \theta_{avg} = & \biggl(\frac{1}{N}\sum_{i=1}^{N}
                          \frac{1}{\theta_{i}}\biggr)^{-1} \\

where :math:`N` is the number of particles in the system.

----------

Output info
"""""""""""

This compute calculates a global vector of length 5 (:math:`U^\text{cond}`,
:math:`U^\text{mech}`, :math:`U^\text{chem}`, :math:`\theta_\text{avg}`,
:math:`N`), which can be accessed by indices 1 through 5.
See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The vector values will be in energy and temperature :doc:`units <units>`.

Restrictions
""""""""""""

This command is part of the DPD-REACT package.  It is only enabled if
LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

This command also requires use of the :doc:`atom_style dpd <atom_style>`
command.

Related commands
""""""""""""""""

:doc:`compute dpd/atom <compute_dpd_atom>`,
:doc:`thermo_style <thermo_style>`

Default
"""""""

none

----------

.. _Larentzos1:

**(Larentzos)** J.P. Larentzos, J.K. Brennan, J.D. Moore, and
W.D. Mattson, "LAMMPS Implementation of Constant Energy Dissipative
Particle Dynamics (DPD-E)", ARL-TR-6863, U.S. Army Research
Laboratory, Aberdeen Proving Ground, MD (2014).

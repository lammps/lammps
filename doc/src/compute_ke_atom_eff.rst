.. index:: compute ke/atom/eff

compute ke/atom/eff command
===========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID ke/atom/eff

* ID, group-ID are documented in :doc:`compute <compute>` command
* ke/atom/eff = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all ke/atom/eff

Description
"""""""""""

Define a computation that calculates the per-atom translational
(nuclei and electrons) and radial kinetic energy (electron only) in a
group.  The particles are assumed to be nuclei and electrons modeled
with the :doc:`electronic force field <pair_eff>`.

The kinetic energy for each nucleus is computed as 1/2 m v\^2, where m
corresponds to the corresponding nuclear mass, and the kinetic energy
for each electron is computed as 1/2 (me v\^2 + 3/4 me s\^2), where me
and v correspond to the mass and translational velocity of each
electron, and s to its radial velocity, respectively.

There is a subtle difference between the quantity calculated by this
compute and the kinetic energy calculated by the *ke* or *etotal*
keyword used in thermodynamic output, as specified by the
:doc:`thermo\_style <thermo_style>` command. For this compute, kinetic
energy is "translational" plus electronic "radial" kinetic energy,
calculated by the simple formula above. For thermodynamic output, the
*ke* keyword infers kinetic energy from the temperature of the system
with 1/2 Kb T of energy for each (nuclear-only) degree of freedom in
eFF.

.. note::

   The temperature in eFF should be monitored via the :doc:`compute temp/eff <compute_temp_eff>` command, which can be printed with
   thermodynamic output by using the :doc:`thermo\_modify <thermo_modify>`
   command, as shown in the following example:


.. parsed-literal::

   compute         effTemp all temp/eff
   thermo_style    custom step etotal pe ke temp press
   thermo_modify   temp effTemp

The value of the kinetic energy will be 0.0 for atoms (nuclei or
electrons) not in the specified compute group.

**Output info:**

This compute calculates a scalar quantity for each atom, which can be
accessed by any command that uses per-atom computes as input.  See the
:doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS
output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-EFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

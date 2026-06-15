.. index:: fix atom_weight/apip

fix atom_weight/apip command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID atom_weight/apip nevery fast_potential precise_potential lambda_input lambda_zone group_lambda_input [no_rescale]

* ID, group-ID are documented in :doc:`fix <fix>` command
* atom_weight/apip = style name of this fix command
* nevery = perform load calculation every this many steps
* fast_potential = *eam* or *ace* for time measurements of the corresponding pair_style or float for constant time
* precise_potential = *ace* for a time measurement of the pair_style pace/apip or float for constant time
* lambda_input = *lambda/input* for a time measurement of pair_style lambda/input/apip or float for constant time
* lambda_zone = *lambda/zone* for a time measurement of pair_style lambda/zone/apip or float for constant time
* group_lambda_input = group-ID of the group for which lambda_input is computed
* no_rescale = do not rescale the work per processor to the measured total force-computation time

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all atom_weight/apip 50 eam ace lambda/input lambda/zone all
   fix 2 all atom_weight/apip 50 1e-05 0.0004 4e-06 4e-06 all
   fix 2 all atom_weight/apip 50 ace ace 4e-06 4e-06 all no_rescale

Description
"""""""""""

This command approximates the load every atom causes when an
adaptive-precision interatomic potential (APIP) according to
:ref:`(Immel) <Immel2025_2>` is used.
This approximated load can be saved as atomic variable and
used as input for the dynamic load balancing via the
:doc:`fix balance <fix_balance>` command.

An adaptive-precision potential like :doc:`eam/apip <pair_eam_apip>`
and :doc:`pace/apip <pair_pace_apip>` is calculated only
for a subset of atoms.
The switching parameter that determines per atom, which potential energy is
used, can be also calculated by
:doc:`pair_style lambda/input/apip <pair_lambda_input_apip>`.
A spatial switching zone, that ensures a smooth transition between two
different interatomic potentials, can be calculated by
:doc:`pair_style lambda/zone/apip <pair_lambda_zone_apip>`.
Thus, there are up to four force-subroutines, that are computed only for a
subset of atoms and combined via the pair_style :doc:`hybrid/overlay <pair_hybrid>`.
For all four force-subroutines, the average work per atom is be measured
per processor by the corresponding pair_style.
This fix extracts these measurements of the pair styles every *nevery*
steps. The average compute times are used to calculates a per-atom vector with
the approximated atomic weight, whereas the average compute time of the four
subroutines contributes only to the load of atoms, for which the corresponding
subroutine was calculated.
If not disabled via *no_rescale*, the so calculated load is
rescaled per processor so that the total atomic compute time matches the
also measured total compute time of the whole pair_style.
This atomic weight is intended to be used
as input for :doc:`fix balance <fix_balance>`:

.. code-block:: LAMMPS

   variable nevery equal 10
   fix weight_atom all atom_weight/apip ${nevery} eam ace lambda/input lambda/zone all
   variable myweight atom f_weight_atom
   fix balance all balance ${nevery} 1.1 rcb weight var myweight

Furthermore, this fix provides the over the processors averaged compute time of the
four pair_styles, which are used to approximate the atomic weight, as vector.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to
:doc:`binary restart files <restart>`.  None of the
:doc:`fix_modify <fix_modify>` options are relevant to this fix.

This fix produces a per-atom vector that contains the atomic
weight of each atom.
The per-atom vector can only be accessed on timesteps that are multiples
of *nevery*.

Furthermore, this fix computes a global vector of length 4 with
statistical information about the four different (possibly)
measured compute times per force subroutine. The four
values in the vector are as follows:

  #. average compute time for one atom using the fast pair_style
  #. average compute time for one atom using the precise pair_style
  #. average compute time of lambda/input/apip for one atom
  #. average compute time of lambda/zone/apip for one atom

The compute times are computed as average of all processors that
measured at least one computation of the corresponding style.
The vector values calculated by this fix are "intensive" and
updated whenever the per-atom vector is computed, i.e., in
timesteps that are multiples of *nevery*.

The vector and the per-atom vector can be accessed by various
:doc:`output commands <Howto_output>`.


No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

This fix is part of the APIP package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix balance <fix_balance>`,
:doc:`fix lambda/apip <fix_lambda_apip>`,
:doc:`fix lambda_thermostat/apip <fix_lambda_thermostat_apip>`,
:doc:`pair_style lambda/zone/apip <pair_lambda_zone_apip>`,
:doc:`pair_style lambda/input/apip  <pair_lambda_input_apip>`,
:doc:`pair_style eam/apip <pair_eam_apip>`,
:doc:`pair_style pace/apip  <pair_pace_apip>`,

Default
"""""""

*no_rescale* is not used by default.

----------

.. _Immel2025_2:

**(Immel)** Immel, Drautz and Sutmann, J Chem Phys, 162, 114119 (2025)

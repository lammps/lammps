.. index:: pair_style lambda/zone/apip

pair_style lambda/zone/apip command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lambda/zone/apip cutoff

* lambda/zone/apip = style name of this pair style
* cutoff = global cutoff (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lambda/zone/apip 12.0

Description
"""""""""""

This pair_style calculates :math:`\lambda_{\text{min},i}`, which
is required for :doc:`fix lambda/apip <fix_lambda_apip>`.
The meaning of :math:`\lambda_{\text{min},i}` is documented in
:doc:`fix lambda/apip <fix_lambda_apip>`, as this pair_style is for use with
:doc:`fix lambda/apip <fix_lambda_apip>` only.

This pair_style requires only the global cutoff as argument.
The remaining quantities, that are required to calculate
:math:`\lambda_{\text{min},i}` are extracted from
:doc:`fix lambda/apip <fix_lambda_apip>` and, thus,
do not need to be passed to this pair_style as arguments.

.. warning::

   The cutoff given as argument to this pair style is only relevant for the
   neighbor list creation. The radii, which define :math:`r_{\lambda,\text{hi}}` and :math:`r_{\lambda,\text{lo}}` are defined by :doc:`fix lambda/apip <fix_lambda_apip>`.

The computation of :math:`\lambda_{\text{min},i}` is done by this
pair_style instead of by :doc:`fix lambda/apip <fix_lambda_apip>`, as this computation
takes time and this pair_style can be included in the load-balancing via
:doc:`fix atom_weight/apip <fix_atom_weight_apip>`.

A code example for the calculation of the switching parameter for an
adaptive-precision interatomic potential (APIP) is given in the following:
The adaptive-precision potential is created
by combining :doc:`pair_style eam/fs/apip <pair_eam_apip>`
and :doc:`pair_style pace/precise/apip <pair_pace_apip>`.
The input, from which the switching parameter is calculated, is provided
by :doc:`pair lambda/input/csp/apip <pair_lambda_input_apip>`.
The switching parameter is calculated by :doc:`fix lambda/apip <fix_lambda_apip>`,
whereas the spatial transition zone of the switching parameter is calculated
by this pair style.

.. code-block:: LAMMPS

   pair_style hybrid/overlay eam/fs/apip pace/precise/apip lambda/input/csp/apip fcc cutoff 5.0 lambda/zone/apip 12.0
   pair_coeff * * eam/fs/apip Cu.eam.fs Cu
   pair_coeff * * pace/precise/apip Cu_precise.yace Cu
   pair_coeff * * lambda/input/csp/apip
   pair_coeff * * lambda/zone/apip
   fix 2 all lambda/apip 3.0 3.5 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The cutoff distance for this pair style can be mixed.  The default mix
value is *geometric*\ .  See the "pair_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style writes no information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands need
to be specified in an input script that reads a restart file.

This pair style does not support the use of the *inner*, *middle*,
and *outer* keywords of the :doc:`run_style respa <run_style>` command.

----------

Restrictions
""""""""""""
This fix is part of the APIP package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix lambda/apip <fix_lambda_apip>`,
:doc:`fix atom_weight/apip <fix_atom_weight_apip>`
:doc:`pair_style lambda/input/apip  <pair_lambda_input_apip>`,
:doc:`pair_style eam/apip <pair_eam_apip>`,
:doc:`pair_style pace/apip  <pair_pace_apip>`,
:doc:`fix lambda_thermostat/apip <fix_lambda_thermostat_apip>`,

Default
"""""""

none

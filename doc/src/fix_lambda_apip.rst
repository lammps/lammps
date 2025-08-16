.. index:: fix lambda/apip

fix lambda/apip command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID lambda/apip thr_lo thr_hi keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lambda/apip = style name of this fix command
* thr_lo = value below which :math:`\lambda_i^\text{input}` results in a switching parameter of 1
* thr_hi = value above which :math:`\lambda_i^\text{input}` results in a switching parameter of 0
* zero or one keyword/args pairs may be appended
* keyword = *time_averaged_zone* or *min_delta_lambda* or *lambda_non_group* or *store_atomic_stats* or *dump_atomic_history* or *group_fast* or *group_precise* or *group_ignore_lambda_input*

  .. parsed-literal::

       *time_averaged_zone* args = cut_lo cut_hi history_len_lambda_input history_len_lambda
         cut_lo = distance at which the radial function decreases from 1
         cut_hi = distance from which on the radial function is 0
         history_len_lambda_input = number of time steps for which lambda_input is averaged
         history_len_lambda = number of time steps for which the switching parameter is averaged
       *min_delta_lambda* args = delta
         delta = value below which changes of the switching parameter are neglected (>= 0)
       *lambda_non_group* args = lambda_ng
         lambda_ng = *precise* or *fast* or float
           *precise* = assign a constant switching parameter of 0 to atoms, that are not in the group specified by group-ID
           *fast* = assign a constant switching parameter of 1 to atoms, that are not in the group specified by group-ID
           float = assign this constant switching parameter to atoms, that are not in the group specified by group-ID (0 <= float <= 1)
       *group_fast* args = group-ID-fast
         group-ID-fast = the switching parameter of 1 is used instead of the one computed by lambda_input for atoms in the group specified by group-ID-fast
       *group_precise* args = group-ID-precise
         group-ID-precise = the switching parameter of 0 is used instead of the one computed by lambda_input for atoms in the group specified by group-ID-precise
       *group_ignore_lambda_input* args = group-ID-ignore-lambda-input
         group-ID-ignore-lambda-input = the switching parameter of lambda_ng is used instead of the one computed by lambda_input for atoms in the group specified by group-ID-ignore-lambda-input
       *store_atomic_stats* args = none
       *dump_atomic_history* args = none

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all lambda/apip 3.0 3.5 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01
   fix 2 mobile lambda/apip 3.0 3.5 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01 group_ignore_lambda_input immobile lambda_non_group fast

Description
"""""""""""
The potential energy :math:`E_i` of an atom :math:`i` of an adaptive-precision
potential according to :ref:`(Immel) <Immel2025_3>` is given by

.. math::

   E_i = \lambda_i E_i^\text{(fast)} + (1-\lambda_i) E_i^\text{(precise)},

whereas :math:`E_i^\text{(fast)}` is the potential energy of atom :math:`i`
according to a fast interatomic potential like EAM,
:math:`E_i^\text{(precise)}` is the potential energy according to a precise
interatomic potential such as ACE and :math:`\lambda_i\in[0,1]` is the
switching parameter that decides which potential energy is used.
This fix calculates the switching parameter :math:`\lambda_i` based on the
input provided from :doc:`pair_style lambda/input/apip <pair_lambda_input_apip>`.

The calculation of the switching parameter is described in detail in
:ref:`(Immel) <Immel2025_3>`.
This fix calculates the switching parameter for all atoms in the
:doc:`group <group>`
described by group-ID, while the value of *lambda_non_group* is used
as switching parameter for all other atoms.

First, this fix calculates per atom :math:`i` the time averaged input
:math:`\lambda^\text{input}_{\text{avg},i}` from
:math:`\lambda^\text{input}_{i}`, whereas the number of averaged timesteps
can be set via *time_averaged_zone*.

.. note::

   :math:`\lambda^\text{input}_{i}` is calculated by
   :doc:`pair_style lambda/input/apip <pair_lambda_input_apip>`, which needs to be included
   in the input script as well.

The time averaged input :math:`\lambda^\text{input}_{\text{avg},i}` is then
used to calculate the switching parameter

.. math::

   \lambda_{0,i}(t) = f^\text{(cut)} \left(\frac{\lambda_{\text{avg},i}^\text{input}(t) - \lambda_\text{lo}^\text{input}}{\lambda_\text{hi}^\text{input} - \lambda_\text{lo}^\text{input}} \right)\,,

whereas the thresholds :math:`\lambda_\text{hi}^\text{input}`
and  :math:`\lambda_\text{lo}^\text{input}` are set by the
values provided as *thr_lo* and *thr_hi* and :math:`f^\text{(cut)}(x)` is a cutoff function
that is 1 for :math:`x\leq 0`, decays from 1 to 0 for :math:`x\in[0,1]`, and
is 0 for :math:`x\geq 1`.
If the *group_precise* argument is used, :math:`\lambda_{0,i}=0` is used for all
atoms :math:`i` assigned to the corresponding :doc:`group <group>`.
If the *group_fast* argument is used, :math:`\lambda_{0,i}=1` is used for all
atoms :math:`i` assigned to the corresponding :doc:`group <group>`.
If an atom is in the groups *group_fast* and *group_precise*,
:math:`\lambda_{0,i}=0` is used.
If the *group_ignore_lambda_input* argument is used,
:math:`\lambda_i^\text{input}` is not computed for all atoms :math:`i` assigned
to the corresponding :doc:`group <group>`; instead, if the value is not already
set by *group_fast* or *group_precise*, the value of *lambda_non_group* is
used.

.. note::

   The computation of :math:`\lambda_i^\text{input}` is not required for
   atoms that are in the groups *group_fast* and *group_precise*.
   Thus, one should use *group_ignore_lambda_input* and prevent the
   computation of :math:`\lambda_i^\text{input}` for all atoms, for
   which a constant input is used.

A spatial transition zone between the fast and the precise potential is
introduced via

.. math::

   \lambda_{\text{min},i}(t) = \text{min}\left(\left\{1 - (1 -\lambda_{0,j}(t)) f^\text{(cut)}\left(\frac{r_{ij}(t)-r_{\lambda,\text{lo}}}{r_{\lambda,\text{hi}} - r_{\lambda,\text{lo}}}\right) : j \in \Omega_{\lambda,i} \right\}\right)\,,

whereas the thresholds :math:`r_{\lambda,\text{lo}}` and
:math:`r_{\lambda,\text{hi}}`
of the cutoff function are set via *time_averaged_zone* and
:math:`\Omega_{\lambda,i}` is the set of
neighboring atoms of atom :math:`i`.

.. note::

   :math:`\lambda_{\text{min},i}` is calculated by
   :doc:`pair_style lambda/zone/apip <pair_lambda_zone_apip>`, which needs to be included
   in the input script as well.

The switching parameter is smoothed by the calculation of the time average

.. math::

   \lambda_{\text{avg},i}(t) = \frac{1}{N_{\lambda,\text{avg}}} \sum_{n=1}^{N_{\lambda,\text{avg}}} \lambda_{\text{min},i}(t - n \Delta t)\,,

whereas :math:`\Delta t` is the :doc:`timestep <timestep>` and
:math:`N_{\lambda,\text{avg}}` is the number of averaged timesteps, that
can be set via *time_averaged_zone*.

Finally, numerical fluctuations of the switching parameter are suppressed by the usage of

.. math::

   \lambda_{i}(t) = \left\{
   \begin{array}{ll}
   \lambda_{\text{avg},i}(t) & \text{ for } \left|\lambda_{\text{avg},i}(t) - \lambda_{i}(t-\Delta t)\right|\geq \Delta\lambda_\text{min} \text{ or } \lambda_{\text{avg},i}(t)\in\{0,1\}, \\
   \lambda_{i}(t-\Delta t) & \text{ otherwise}\,,
   \end{array}
   \right.

whereas the minimum change :math:`\Delta\lambda_\text{min}` is set by the
*min_delta_lambda* argument.

.. note::

   *group_fast* affects only :math:`\lambda_{0,i}(t)`. The switching parameter
   of atoms in this :doc:`group <group>` may change due to the calculation of the
   spatial switching zone.
   A switching parameter of 1 can be enforced by excluding the corresponding
   atoms from the :doc:`group <group>` described by group-ID and using *lambda_non_group* 1
   as argument.

----------

A code example for the calculation of the switching parameter for an
adaptive-precision potential is given in the following:
The adaptive-precision potential is created
by combining :doc:`pair_style eam/fs/apip <pair_eam_apip>`
and :doc:`pair_style pace/precise/apip <pair_pace_apip>`.
The input, from which the switching parameter is calculated, is provided
by :doc:`pair lambda/input/csp/apip <pair_lambda_input_apip>`.
The switching parameter is calculated by this fix, whereas the spatial
transition zone of the switching parameter is calculated by
:doc:`pair_style lambda/zone/apip <pair_lambda_zone_apip>`.

.. code-block:: LAMMPS

   pair_style hybrid/overlay eam/fs/apip pace/precise/apip lambda/input/csp/apip fcc cutoff 5.0 lambda/zone/apip 12.0
   pair_coeff * * eam/fs/apip Cu.eam.fs Cu
   pair_coeff * * pace/precise/apip Cu_precise.yace Cu
   pair_coeff * * lambda/input/csp/apip
   pair_coeff * * lambda/zone/apip
   fix 2 all lambda/apip 3.0 3.5 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01


----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The saved history of the switching parameter :math:`\lambda_i`
and the saved history of
:math:`\lambda_i^\text{input}` are written to
:doc:`binary restart files <restart>` allow a smooth restart of a simulation.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.

If the *store_atomic_stats* argument is used, basic statistics is provided as
per-atom array:

  #. :math:`\lambda_i^\text{input}(t)`
  #. :math:`\lambda_{\text{avg},i}^\text{input}(t)`
  #. :math:`\lambda_{0,i}(t)`
  #. :math:`\lambda_{\text{min},i}(t)`
  #. :math:`\lambda_{i}(t)`

If the *dump_atomic_history* argument is used, the whole saved history
of :math:`\lambda_i^\text{input}(t)` is appended to the previously
mentioned array per atom.

The per-atom vector can be accessed by various
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

:doc:`pair_style lambda/zone/apip <pair_lambda_zone_apip>`,
:doc:`pair_style lambda/input/apip  <pair_lambda_input_apip>`,
:doc:`pair_style eam/apip <pair_eam_apip>`,
:doc:`pair_style pace/apip  <pair_pace_apip>`,
:doc:`fix atom_weight/apip <fix_atom_weight_apip>`
:doc:`fix lambda_thermostat/apip <fix_lambda_thermostat_apip>`,

Default
"""""""

*min_delta_lambda* = 0,
*lambda_non_group* = 1,
*cut_lo* = 4.0,
*cut_hi* = 12.0,
*history_len_lambda_input* = 100,
*history_len_lambda* = 100,
*store_atomic_stats* is not used,
*dump_atomic_history* is not used,
*group_fast* is not used,
*group_precise* is not used,
*group_ignore_lambda_input* is not used

----------

.. _Immel2025_3:

**(Immel)** Immel, Drautz and Sutmann, J Chem Phys, 162, 114119 (2025)

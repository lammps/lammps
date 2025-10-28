Adaptive-precision interatomic potentials (APIP)
================================================

The :ref:`PKG-APIP <PKG-APIP>` enables use of adaptive-precision potentials
as described in :ref:`(Immel) <Immel2025_1>`.
In the context of this package, precision refers to the accuracy of an interatomic
potential.

Modern machine-learning (ML) potentials translate the accuracy of DFT
simulations into MD simulations, i.e., ML potentials are more accurate
compared to traditional empirical potentials.
However, this accuracy comes at a cost: there is a considerable performance
gap between the evaluation of classical and ML potentials, e.g., the force
calculation of a classical EAM potential is 100-1000 times faster compared
to the ML-based ACE method.
The evaluation time difference results in a conflict between large time and
length scales on the one hand and accuracy on the other.
This conflict is resolved by an APIP model for simulations, in which the highest precision
is required only locally but not globally.

An APIP model uses a precise but
expensive ML potential only for a subset of atoms, while a fast
potential is used for the remaining atoms.
Whether the precise or the fast potential is used is determined
by a continuous switching parameter :math:`\lambda_i` that can be defined for each
atom :math:`i`.
The switching parameter can be adjusted dynamically during a simulation or
kept constant as explained below.

The potential energy :math:`E_i` of an atom :math:`i` described by an
adaptive-precision
interatomic potential is given by :ref:`(Immel) <Immel2025_1>`

.. math::

   E_i = \lambda_i E_i^\text{(fast)} + (1-\lambda_i) E_i^\text{(precise)},

whereas :math:`E_i^\text{(fast)}` is the potential energy of atom :math:`i`
according to a fast interatomic potential,
:math:`E_i^\text{(precise)}` is the potential energy according to a precise
interatomic potential and :math:`\lambda_i\in[0,1]` is the
switching parameter that decides how the potential energies are weighted.

Adaptive-precision saves computation time when the computation of the
precise potential is not required for many atoms, i.e., when
:math:`\lambda_i=1` applies for many atoms.

The currently implemented potentials are:

.. list-table::
   :header-rows: 1

   * - Fast potential
     - Precise potential
   * - :doc:`ACE <pair_pace_apip>`
     - :doc:`ACE <pair_pace_apip>`
   * - :doc:`EAM <pair_eam_apip>`
     -

In theory, any short-range potential can be used for an adaptive-precision
interatomic potential. How to implement a new (fast or precise)
adaptive-precision
potential is explained in :ref:`here <implementing_new_apip_styles>`.

The switching parameter :math:`\lambda_i` that combines the two potentials
can be dynamically calculated during a
simulation.
Alternatively, one can set a constant switching parameter before the start
of a simulation.
To run a simulation with an adaptive-precision potential, one needs the
following components:

.. tabs::

   .. tab:: dynamic switching parameter

        #. :doc:`atom_style apip <atom_style>` so that the switching parameter :math:`\lambda_i` can be stored.
        #. A fast potential: :doc:`eam/apip <pair_eam_apip>` or :doc:`pace/fast/apip <pair_pace_apip>`.
        #. A precise potential: :doc:`pace/precise/apip <pair_pace_apip>`.
        #. :doc:`pair_style lambda/input/apip  <pair_lambda_input_apip>` to calculate :math:`\lambda_i^\text{input}`, from which :math:`\lambda_i` is calculated.
        #. :doc:`fix lambda/apip <fix_lambda_apip>` to calculate the switching parameter :math:`\lambda_i`.
        #. :doc:`pair_style lambda/zone/apip <pair_lambda_zone_apip>` to calculate the spatial transition zone of the switching parameter.
        #. :doc:`pair_style hybrid/overlay <pair_hybrid>` to combine the previously mentioned pair_styles.
        #. :doc:`fix lambda_thermostat/apip <fix_lambda_thermostat_apip>` to conserve the energy when switching parameters change.
        #. :doc:`fix atom_weight/apip <fix_atom_weight_apip>` to approximate the load caused by every atom, as the computations of the pair_styles are only required for a subset of atoms.
        #. :doc:`fix balance <fix_balance>` to perform dynamic load balancing with the calculated load.

   .. tab:: constant switching parameter

        #. :doc:`atom_style apip <atom_style>` so that the switching parameter :math:`\lambda_i` can be stored.
        #. A fast potential: :doc:`eam/apip <pair_eam_apip>` or :doc:`pace/fast/apip <pair_pace_apip>`.
        #. A precise potential: :doc:`pace/precise/apip <pair_pace_apip>`.
        #. :doc:`set <set>` command to set the switching parameter :math:`\lambda_i`.
        #. :doc:`pair_style hybrid/overlay <pair_hybrid>` to combine the previously mentioned pair_styles.
        #. :doc:`fix atom_weight/apip <fix_atom_weight_apip>` to approximate the load caused by every atom, as the computations of the pair_styles are only required for a subset of atoms.
        #. :doc:`fix balance <fix_balance>` to perform dynamic load balancing with the calculated load.

----------

Example
"""""""
.. note::

   How to select the values of the parameters of an adaptive-precision
   interatomic potential is discussed in detail in :ref:`(Immel) <Immel2025_1>`.


.. tabs::

   .. tab:: dynamic switching parameter

      Lines like these would appear in the input script:


      .. code-block:: LAMMPS

         atom_style apip
         comm_style tiled

         pair_style hybrid/overlay eam/fs/apip pace/precise/apip lambda/input/csp/apip fcc cutoff 5.0 lambda/zone/apip 12.0
         pair_coeff * * eam/fs/apip Cu.eam.fs Cu
         pair_coeff * * pace/precise/apip Cu.yace Cu
         pair_coeff * * lambda/input/csp/apip
         pair_coeff * * lambda/zone/apip

         fix 2 all lambda/apip 2.5 3.0 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01
         fix 3 all lambda_thermostat/apip N_rescaling 200
         fix 4 all atom_weight/apip 100 eam ace lambda/input lambda/zone all

         variable myweight atom f_4

         fix 5 all balance 100 1.1 rcb weight var myweight

      First, the :doc:`atom_style apip <atom_style>` and the communication style are set.

      .. note::
         Note, that :doc:`comm_style <comm_style>` *tiled* is required for the style *rcb* of
         :doc:`fix balance <fix_balance>`, but not for APIP.
         However, the flexibility offered by the balancing style *rcb*, compared to the
         balancing style *shift*, is advantageous for APIP.

      An adaptive-precision EAM-ACE potential, for which the switching parameter
      :math:`\lambda` is calculated from the CSP, is defined via
      :doc:`pair_style hybrid/overlay <pair_hybrid>`.
      The fixes ensure that the switching parameter is calculated, the energy conserved,
      the weight for the load balancing calculated and the load-balancing itself is done.

   .. tab:: constant switching parameter

      Lines like these would appear in the input script:

      .. code-block:: LAMMPS

         atom_style apip
         comm_style tiled

         pair_style hybrid/overlay eam/fs/apip pace/precise/apip
         pair_coeff * * eam/fs/apip Cu.eam.fs Cu
         pair_coeff * * pace/precise/apip Cu.yace Cu

         # calculate lambda somehow
         variable lambda atom ...
         set group all apip/lambda v_lambda

         fix 4 all atom_weight/apip 100 eam ace lambda/input lambda/zone all

         variable myweight atom f_4

         fix 5 all balance 100 1.1 rcb weight var myweight

      First, the :doc:`atom_style apip <atom_style>` and the communication style are set.

      .. note::
         Note, that :doc:`comm_style <comm_style>` *tiled* is required for the style *rcb* of
         :doc:`fix balance <fix_balance>`, but not for APIP.
         However, the flexibility offered by the balancing style *rcb*, compared to the
         balancing style *shift*, is advantageous for APIP.

      An adaptive-precision EAM-ACE potential is defined via
      :doc:`pair_style hybrid/overlay <pair_hybrid>`.
      The switching parameter :math:`\lambda_i` of the adaptive-precision
      EAM-ACE potential is set via the :doc:`set command <set>`.
      The parameter is not updated during the simulation.
      Therefore, the potential is conservative.
      The fixes ensure that the weight for the load balancing is calculated
      and the load-balancing itself is done.

----------

.. _implementing_new_apip_styles:

Implementing new APIP pair styles
"""""""""""""""""""""""""""""""""

One can introduce adaptive-precision to an existing pair style by modifying
the original pair style.
One should calculate the force
:math:`F_i = - \nabla_i \sum_j E_j^\text{original}` for a fast potential or
:math:`F_i = - (1-\nabla_i) \sum_j E_j^\text{original}` for a precise
potential from the original potential
energy :math:`E_j^\text{original}` to see where the switching parameter
:math:`\lambda_i` needs to be introduced in the force calculation.
The switching parameter :math:`\lambda_i` is known for all atoms :math:`i`
in force calculation routine.
One needs to introduce an abortion criterion based on :math:`\lambda_i` to
ensure that all not required calculations are skipped and compute time can
be saved.
Furthermore, one needs to provide the number of calculations and measure the
computation time.
Communication within the force calculation needs to be prevented to allow
effective load-balancing.
With communication, the load balancer cannot balance few calculations of the
precise potential on one processor with many computations of the fast
potential on another processor.

All changes in the pair_style pace/apip compared to the pair_style pace
are annotated and commented.
Thus, the pair_style pace/apip can serve as an example for the implementation
of new adaptive-precision potentials.

----------

.. _Immel2025_1:

**(Immel)** Immel, Drautz and Sutmann, J Chem Phys, 162, 114119 (2025)

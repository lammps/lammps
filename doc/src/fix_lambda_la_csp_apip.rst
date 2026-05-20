.. index:: fix lambda/la/csp/apip

fix lambda/la/csp/apip command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID lambda/la/csp/apip thr_lo thr_hi cut_lo cut_hi lattice keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lambda/la/csp/apip = style name of this fix command
* thr_lo = value below which the differentiable CSP results in a switching parameter of 1 (squared distance units)
* thr_hi = value above which the differentiable CSP results in a switching parameter of 0 (squared distance units)
* cut_lo = radius below which the weighting function for the local averaging is 1 (distance units)
* cut_hi = radius above which the weighting function for the local averaging is 0 (distance units)
* lattice = *fcc* or *bcc* or N = # of neighbors per atom to include in the CSP calculation
* zero or more keyword/args pairs may be appended
* keyword = *csp_cut* or *csp_mode* or *forces* or *lambda_non_group* or *store_peratom*

  .. parsed-literal::

       *csp_cut* args = float
           float = neighboring atoms outside of this cutoff radius may not be considered for the CSP calculation
       *csp_mode* args = *dynamic* or *static*
           *dynamic* = use the differentiable CSP to calculate the switching parameter
           *static* = use the non-differentiable CSP instead of the differentiable one
       *forces* args = *no* or *yes*
           *yes* = compute the forces caused by the differentiation of the switching parameter
           *no* = do not compute the forces caused by the differentiation of the switching parameter
       *lambda_non_group* args = *precise* or *fast* or float
           *precise* = assign a constant switching parameter of 0 to atoms, that are not in the group specified by group-ID
           *fast* = assign a constant switching parameter of 1 to atoms, that are not in the group specified by group-ID
           float = assign this constant switching parameter to atoms, that are not in the group specified by group-ID (0 <= float <= 1)
       *store_peratom* args = integer
           integer = provide per-atom output every this many timesteps

Examples
""""""""

.. code-block:: LAMMPS

   fix lambda_la all lambda/la/csp/apip 0.25 1.5 15.0 16.0 bcc
   fix lambda_la mobile lambda/la/csp/apip 0.24 1.5 11.0 12.0 bcc lambda_non_group fast

Description
"""""""""""

.. versionadded:: 30Mar2026

The potential energy :math:`E_i` of an atom :math:`i` according to an
adaptive-precision potential is given by :ref:`(Immel2025) <Immel2025_8>`

.. math::

   E_i = \lambda_i E_i^\text{(fast)} + (1-\lambda_i) E_i^\text{(precise)}\,,

where :math:`E_i^\text{(fast)}` is the potential energy of atom
:math:`i` according to a fast computable interatomic potential,
:math:`E_i^\text{(precise)}` is the potential energy according to a
precise interatomic potential and :math:`\lambda_i\in[0,1]` is the
switching parameter that decides which potential energy is used.

This fix calculates the switching parameter :math:`\lambda_i` based on
local averaging of a descriptor according to :ref:`(Immel2026)
<Immel2026_1>` and parts of the conservatively calculated force as will
be discussed later.  The descriptor is averaged within the cutoff radius
provided as *cut_hi*.

Per default, a differentiable version of the centro-symmetry parameter
(CSP) is used as descriptor. This differentiable version is described in
detail in :ref:`(Immel2026) <Immel2026_1>`.  The usage of a
differentiable CSP results in a conservative potential, that conserves
(in the absence of external forces) energy and momentum by design.  The
force :math:`\pmb{F}_i=-\nabla_i\sum_kE_k` following from the
adaptive-precision potential energy is given by

.. math::
   \pmb{F}_i  = \sum_k \big(- \lambda_k \nabla_i
   E_k^\text{(fast)} - (1 - \lambda_k) \nabla_i E_k^\text{(precise)}
   + (\nabla_i\lambda_k) (E_k^\text{(precise)} - E_k^\text{(fast)})\big)\,.

This fix calculates the terms :math:`(\nabla_i\lambda_k)
(E_k^\text{(precise)} - E_k^\text{(fast)})` based on potential energies
that are provided by pair styles.  This force-calculation is enabled by
default, but prevented by *forces = no*. Thus, one can use this fix to
only calculate switching parameters.

.. note::

   The fast potential and the precise potential are combined via
   :doc:`pair_style hybrid/overlay <pair_hybrid>` as shown in the code
   example below.

The original CSP :ref:`(Kelchner) <Kelchner_3>` is used instead of the
differentiable CSP with the option *csp_mode dynamic*.

.. warning::

   Note that the original CSP is not differentiable and does not result
   in a conservative potential.

This fix calculates the switching parameter for all atoms in the
:doc:`group <group>` described by group-ID, while the value of
*lambda_non_group* is used as switching parameter for all other atoms.


----------

A code example for the usage of a conservative adaptive-precision
interatomic potential is given in the following: This fix calculates the
switching parameter.  The :doc:`pair_style hybrid/overlay <pair_hybrid>`
is used to combine the two pair styles :doc:`pair_style eam/fs/apip
<pair_eam_apip>` and :doc:`pair_style pace/precise/apip
<pair_pace_apip>`, which calculate the potential energies
:math:`E_k^\text{(precise)}` and :math:`E_k^\text{(fast)}`.  The
conservative force is calculated by both the fix and the pair styles.

.. code-block:: LAMMPS

   fix lambda_la all lambda/la/csp/apip 0.24 1.5 11.0 12.0 bcc
   pair_style hybrid/overlay eam/fs/apip pace/apip
   pair_coeff * * eam/fs/apip W.eam.fs W
   pair_coeff * * pace/apip W.yace W

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The CSP neighbors of the differentiable CSP are saved and written to
:doc:`binary restart files <restart>` to allow a smooth restart of a
simulation.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

This fix can compute the force caused by the differentiation of the
switching parameters.  Therefore, the :doc:`fix_modify <fix_modify>`
*virial* option is supported by this fix to add the contribution of
these forces to both the global pressure and per-atom stress of the
system via the :doc:`compute pressure <compute_pressure>` and
:doc:`compute stress/atom <compute_stress_atom>` commands.  The former
can be accessed by :doc:`thermodynamic output <thermo_style>`.  The
default setting for this fix is :doc:`fix_modify virial yes
<fix_modify>`.

If the (non-conservative) dynamic CSP-pairs are used, this fix returns
the accumulated number of changed CSP-pairs as global scalar.  The
potential is conservative as long as this number of changed CSP-pairs is
constant (every changed CSP pair may change the total energy of the
system).  The global scalar can be accessed by various :doc:`output
commands <Howto_output>`.

If *store_peratom* is used, 5 quantities are provided as per-atom
vector: 1.-3. the force that is caused by the differentiation of the
switching parameter, 4. the differentiable CSP and 5. the locally
averaged differentiable CSP that was used to calculate the switching
parameter.  The per-atom vector can be accessed by various :doc:`output
commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

This fix is part of the APIP package. It is only enabled if LAMMPS was
built with that package. See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`pair_style eam/apip <pair_eam_apip>`,
:doc:`pair_style pace/apip  <pair_pace_apip>`

Default
"""""""

*forces* = *yes*,
*csp_cut* = 5.0,
*csp_mode* = *static*,
*lambda_non_group* = 1,
*store_peratom* = :math:`\infty`

----------

.. _Immel2025_8:

**(Immel2025)** Immel, Drautz and Sutmann, J Chem Phys, 162, 114119 (2025)

.. _Immel2026_1:

**(Immel2026)** Immel, Drautz and Sutmann, arXiv:2512.07693

.. _Kelchner_3:

**(Kelchner)** Kelchner, Plimpton, Hamilton, Phys Rev B, 58, 11085 (1998)

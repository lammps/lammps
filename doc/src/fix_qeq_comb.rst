.. index:: fix qeq/comb
.. index:: fix qeq/comb/omp

fix qeq/comb command
====================

Accelerator Variants: *qeq/comb/omp*

Syntax
""""""

.. parsed-literal::

   fix ID group-ID qeq/comb Nevery precision keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* qeq/comb = style name of this fix command
* Nevery = perform charge equilibration every this many steps
* precision = convergence criterion for charge equilibration
* zero or more keyword/value pairs may be appended
* keyword = *file*

  .. parsed-literal::

       *file* value = filename
         filename = name of file to write QEQ equilibration info to

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 surface qeq/comb 10 0.0001

Description
"""""""""""

Perform charge equilibration (QeQ) in conjunction with the COMB
(Charge-Optimized Many-Body) potential as described in
:ref:`(COMB_1) <COMB_1>` and :ref:`(COMB_2) <COMB_2>`.  It performs the charge
equilibration portion of the calculation using the so-called QEq
method, whereby the charge on each atom is adjusted to minimize the
energy of the system.  This fix can only be used with the COMB
potential; see the :doc:`fix qeq/reaxff <fix_qeq_reaxff>` command for a QeQ
calculation that can be used with any potential.

Only charges on the atoms in the specified group are equilibrated.
The fix relies on the pair style (COMB in this case) to calculate the
per-atom electronegativity (effective force on the charges).  An
electronegativity equalization calculation (or QEq) is performed in an
iterative fashion, which in parallel requires communication at each
iteration for processors to exchange charge information about nearby
atoms with each other.  See :ref:`Rappe_and_Goddard <Rappe_and_Goddard>` and
:ref:`Rick_and_Stuart <Rick_and_Stuart>` for details.

During a run, charge equilibration is performed every *Nevery* time
steps.  Charge equilibration is also always enforced on the first step
of each run.  The *precision* argument controls the tolerance for the
difference in electronegativity for all atoms during charge
equilibration.  *Precision* is a trade-off between the cost of
performing charge equilibration (more iterations) and accuracy.

If the *file* keyword is used, then information about each
equilibration calculation is written to the specified file.

.. note::

   In order to solve the self-consistent equations for electronegativity
   equalization, LAMMPS imposes the additional constraint that all the
   charges in the fix group must add up to zero.  The initial charge
   assignments should also satisfy this constraint.  LAMMPS will print a
   warning if that is not the case.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is performing charge equilibration. Default is
the outermost level.

This fix produces a per-atom vector which can be accessed by various
:doc:`output commands <Howto_output>`.  The vector stores the gradient
of the charge on each atom.  The per-atom values be accessed on any
timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix can be invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix command currently only supports :doc:`pair style *comb*\ <pair_comb>`.

Related commands
""""""""""""""""

:doc:`pair_style comb <pair_comb>`

Default
"""""""

No file output is performed.

----------

.. _COMB_1:

**(COMB_1)** J. Yu, S. B. Sinnott, S. R. Phillpot, Phys Rev B, 75, 085311 (2007),

.. _COMB_2:

**(COMB_2)** T.-R. Shan, B. D. Devine, T. W. Kemper, S. B. Sinnott, S. R.
Phillpot, Phys Rev B, 81, 125328 (2010).

.. _Rappe_and_Goddard:

**(Rappe_and_Goddard)** A. K. Rappe, W. A. Goddard, J Phys Chem 95, 3358
(1991).

.. _Rick_and_Stuart:

**(Rick_and_Stuart)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chem Phys
101, 16141 (1994).

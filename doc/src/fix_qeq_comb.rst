.. index:: fix qeq/comb

fix qeq/comb command
====================

fix qeq/comb/omp command
========================

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
:ref:`(COMB\_1) <COMB_1>` and :ref:`(COMB\_2) <COMB_2>`.  It performs the charge
equilibration portion of the calculation using the so-called QEq
method, whereby the charge on each atom is adjusted to minimize the
energy of the system.  This fix can only be used with the COMB
potential; see the :doc:`fix qeq/reax <fix_qeq_reax>` command for a QeQ
calculation that can be used with any potential.

Only charges on the atoms in the specified group are equilibrated.
The fix relies on the pair style (COMB in this case) to calculate the
per-atom electronegativity (effective force on the charges).  An
electronegativity equalization calculation (or QEq) is performed in an
iterative fashion, which in parallel requires communication at each
iteration for processors to exchange charge information about nearby
atoms with each other.  See :ref:`Rappe\_and\_Goddard <Rappe_and_Goddard>` and
:ref:`Rick\_and\_Stuart <Rick_and_Stuart>` for details.

During a run, charge equilibration is performed every *Nevery* time
steps.  Charge equilibration is also always enforced on the first step
of each run.  The *precision* argument controls the tolerance for the
difference in electronegativity for all atoms during charge
equilibration.  *Precision* is a trade-off between the cost of
performing charge equilibration (more iterations) and accuracy.

If the *file* keyword is used, then information about each
equilibration calculation is written to the specified file.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

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

.. _COMB\_1:

**(COMB\_1)** J. Yu, S. B. Sinnott, S. R. Phillpot, Phys Rev B, 75, 085311 (2007),

.. _COMB\_2:

**(COMB\_2)** T.-R. Shan, B. D. Devine, T. W. Kemper, S. B. Sinnott, S. R.
Phillpot, Phys Rev B, 81, 125328 (2010).

.. _Rappe\_and\_Goddard:

**(Rappe\_and\_Goddard)** A. K. Rappe, W. A. Goddard, J Phys Chem 95, 3358
(1991).

.. _Rick\_and\_Stuart:

**(Rick\_and\_Stuart)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chem Phys
101, 16141 (1994).

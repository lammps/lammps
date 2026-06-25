.. index:: fix ilves

fix ilves command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ilves tol iter N selectors ... [keyword args ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* tol = convergence tolerance on the relative bond-length error
* iter = maximum number of Newton iterations per time step
* N = print constraint statistics every this many time steps (0 = never)
* selectors = one or more of *b*, *t*, *m*, each followed by one or more numeric values

  .. parsed-literal::

       *b* values = one or more bond types
       *t* values = one or more atom types
       *m* values = one or more atom masses

* zero or more keyword/value pairs may be appended

  .. parsed-literal::

       *variant* value = *fast* or *full*
         *fast* = symmetric quasi-Newton with banded LDLT (default)
         *full* = exact-Newton (asymmetric) with LU decomposition
       *mode* value = *converge* or *fixed*
         *converge* = iterate until the tolerance is met (default)
         *fixed* = always perform *iter* Newton iterations
       *store* value = *yes* or *no*
         *yes* exposes the per-atom constraint forces via a per-atom array

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves 1.0e-6 25 0 b 4 6 8 10
   fix 2 wat ilves 1.0e-8 25 1000 t 1 m 1.008 store yes
   fix 3 all ilves 1.0e-6 25 0 b 1 variant full

Description
"""""""""""

.. versionadded:: TBD

Apply bond-length constraints using the ILVES algorithm of
:ref:`(Lopez-Villellas) <Lopez-Villellas2025>`.  ILVES enforces holonomic
distance constraints with Newton's method on a sparse system of nonlinear
equations.  Unlike :doc:`fix shake <fix_shake>`, ILVES handles arbitrarily
large connected constraint clusters --- for example all the C-C backbone bonds
of a long polymer or protein chain --- in a single solve.

This command is a LAMMPS port of the reference ILVES implementation that the
algorithm authors integrated into GROMACS; the constraint solver itself
(the parallel Schur-complement sparse direct solver and the constraint
topology handling) is reused largely unchanged, while the interface to the
LAMMPS data structures and the time integration follow :doc:`fix shake
<fix_shake>` and :doc:`fix rattle <fix_shake>`.

.. note::

   In this version *fix ilves* constrains bond lengths only.  Angle
   constraints (imposed as a "virtual" distance constraint between the two
   outer atoms of an angle) are not yet supported.

User interface
^^^^^^^^^^^^^^

The user interface of *fix ilves* follows that of :doc:`fix shake
<fix_shake>`.  The three arguments after the *ilves* keyword are the
tolerance :math:`\frac{|g_k|}{d_k^2}` (where :math:`g_k = \frac{1}{2}(|s_k|^2
- d_k^2)` and :math:`d_k` is the target length of constraint *k*), the maximum
number of Newton iterations per step, and the frequency of statistics output
(0 turns it off).

Then one or more groups of selectors (``b``, ``t``, ``m`` lists) pick which
bonds get constrained.  A bond is constrained when **both** of its atoms are
in the fix group AND at least one of the selectors matches:

* the bond type is in the *b* list, or
* either atom type is in the *t* list, or
* either atom mass is within 0.1 mass units of any value in the *m* list.

For each bond that is selected, *fix ilves* sets the corresponding
``bond_type`` to its negative value so that the configured :doc:`bond_style
<bond_style>` skips the (now rigid) interaction and thus avoids
double-counting of bonded forces.  This mirrors how :doc:`fix shake
<fix_shake>` handles the same problem and is reversed automatically when the
fix is deleted.

Algorithm
^^^^^^^^^

At each time step, after the force computation, *fix ilves*:

1. predicts the unconstrained position
   :math:`x'_i = x_i + \Delta t\, v_i + \frac{\Delta t^2}{m_i} f_i` for every
   atom (identical to :doc:`fix shake <fix_shake>`),
2. solves the constraint equations :math:`g_k(x') = 0` with Newton's method,
   using a sparse direct (Schur-complement) solver, and
3. converts the resulting Lagrange multipliers into constraint forces that are
   added to ``atom->f``, so that the following integration step produces the
   constrained positions.

After the velocity update, a RATTLE-style velocity projection removes the
component of relative velocity along each constrained bond, so that both the
bond lengths and their time derivatives satisfy the constraints.

The constraint topology is distributed: each MPI rank builds its constraint
list from its own local bond storage plus the ghost atoms in the standard
communication shell.  The Newton iteration is global --- its convergence test
is the maximum relative bond-length violation reduced over all ranks --- and
the cross-rank coupling is resolved by communicating the predicted positions
to the ghosts and reverse-summing the per-atom corrections to their owners
each iteration.  The results are therefore independent of the number of MPI
ranks (to the solver tolerance).

Solver variants
^^^^^^^^^^^^^^^

Two solver variants are available via the ``variant`` keyword:

* ``fast`` (default): a structurally symmetric, quasi-Newton step whose
  Jacobian is constant within a time step and is therefore factored only once
  per step with a banded LDLT factorization.  This is the faster choice for
  the common case.

* ``full``: the exact Newton Jacobian, which changes every iteration and is
  re-factored with an LU decomposition each iteration.  This converges in
  fewer iterations for strongly coupled or ill-conditioned clusters at a
  higher per-iteration cost.

Both variants converge to the same solution of the constraint equations.

By default the Newton iteration runs until the global maximum relative
bond-length violation falls below *tol* (or *iter* iterations are reached).
With ``mode fixed`` the solver instead performs exactly *iter* iterations every
step.  This removes the per-iteration global reduction used by the convergence
test and can improve parallel efficiency, at the cost of not guaranteeing the
tolerance.  Because ILVES converges quadratically, a small fixed count is
usually sufficient; run first with the default ``mode converge`` and a nonzero
statistics interval *N* to see how many iterations are actually needed (the
maximum is reported with the statistics) before choosing a fixed count.

Output info
^^^^^^^^^^^

When ``N > 0``, every ``N`` time steps the fix prints a summary line per
constrained bond type giving the count, the average constrained length, and
the spread (max - min) of the lengths across all constraints of that type.

With the ``store yes`` keyword, this fix exposes a *per-atom array* with 3
columns containing the constraint-force components ``(fx, fy, fz)`` added by
the fix on the current time step, accessible as ``f_<ID>[1]``, ``f_<ID>[2]``,
and ``f_<ID>[3]``.

By default the constraint forces are not added to the pressure virial; use
:doc:`fix_modify <fix_modify>` *virial yes* to include their contribution to the
global pressure.

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`; the constraint topology is rebuilt from ``atom->bond_atom`` at the
first reneighbor.  Redeclare the same ``fix ilves`` command after
:doc:`read_restart <read_restart>` to restore the constraints (the negated
bond types are preserved in the restart file).

The :doc:`fix_modify <fix_modify>` *virial* option is supported.  This fix
computes no global or per-atom energy and is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the RIGID package.  It is only enabled if LAMMPS was built
with that package.  See the :doc:`Build package <Build_package>` page for more
info.

The molecular topology (bonds) must be defined; an :doc:`atom_style
<atom_style>` such as *full*, *molecular*, or *bond* (or a hybrid atom style
including them) is required, and an :doc:`atom map <atom_modify>` must be
enabled.

Only one ``fix ilves`` instance may be defined at a time.  ``fix ilves`` and
:doc:`fix shake <fix_shake>` must not be used together for overlapping sets of
constrained atoms.

``fix ilves`` does not support :doc:`run_style respa <run_style>`.

All atoms of a constraint cluster must lie within the communication cutoff of
each other on every rank.  For small clusters (water, methyl, hydrogen-only
constraints) this is satisfied automatically; for clusters that span large
distances increase the cutoff with :doc:`comm_modify cutoff <comm_modify>`.
The fix stops with an error if a constraint partner is not available locally.

``fix ilves`` does not support dynamic topology.  Fixes or commands that add,
remove, or change constrained bonds during a run (for example :doc:`fix
bond/create <fix_bond_create>`, :doc:`fix bond/break <fix_bond_break>`, or
:doc:`fix bond/react <fix_bond_react>`) must not be applied to the constrained
atoms.

Related commands
""""""""""""""""

:doc:`fix shake <fix_shake>`, :doc:`fix rattle <fix_shake>`,
:doc:`fix restrain <fix_restrain>`, :doc:`fix rigid <fix_rigid>`

Default
"""""""

The keyword defaults are *variant* = *fast*, *mode* = *converge*, and
*store* = *no*.

----------

.. _Lopez-Villellas2025:

**(Lopez-Villellas)** L. Lopez-Villellas, C. C. K. Mikkelsen,
J. J. Galano-Frutos, S. Marco-Sola, J. Alastruey-Bende, P. Ibanez,
P. Echenique, M. Moreto, M. C. De Rosa, and P. Garcia-Risueno, "ILVES:
Accurate and Efficient Bond Length and Angle Constraints in Molecular
Dynamics", J. Chem. Theory Comput. 21, 8711-8719 (2025),
https://doi.org/10.1021/acs.jctc.5c01376

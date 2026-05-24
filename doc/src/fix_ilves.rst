.. index:: fix ilves

fix ilves command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ilves tol iter N selectors ... [keyword args ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* ilves = style name of this fix command
* tol = convergence tolerance on relative bond-length error
  (``|g_k| / d_k^2`` where ``g_k = 0.5*(|s_k|^2 - d_k^2)``)
* iter = maximum number of Newton iterations per time step
* N = print constraint statistics every this many time steps (0 = never)
* selectors = one or more of *b*, *a*, *t*, *m*, each followed by one or
  more numeric values

  .. parsed-literal::

       *b* values = one or more bond types
       *a* values = one or more angle types
       *t* values = one or more atom types
       *m* values = one or more atom masses

* zero or more keyword/value pairs may be appended

  .. parsed-literal::

       *kbond* value = K
         K = spring constant (energy/distance^2 units) used when substituting
             constraints with harmonic restraints during minimization
       *store* value = *yes* or *no*
         *yes* exposes per-atom constraint forces via array_atom (3 columns)
       *variant* value = *ilves* or *ilvesf*
         *ilves*  = full ILVES (true Newton iteration on the exact Jacobian)
         *ilvesf* = ILVES-F (structurally-symmetric quasi-Newton, faster)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves 1.0e-4 25 100 b 4 6 8 10 12 14 18 a 31
   fix 2 wat ilves 1.0e-6 30 1000 t 1 m 1.008 store yes
   fix 3 all ilves 1.0e-5 20 0 b 1 a 1 variant ilves

Description
"""""""""""

.. versionadded:: TBD

Apply bond-length (and optionally end-to-end angle "virtual bond")
constraints using the ILVES algorithm of
:ref:`(Garcia-Risueno) <Garcia-Risueno2026>`.  ILVES enforces holonomic
distance constraints using Newton's method on a sparse system of nonlinear
equations.  Unlike :doc:`fix shake <fix_shake>`, ILVES handles arbitrarily
large connected constraint clusters --- including the full set of C-C
backbone bonds of a long polymer or protein chain --- in a single solve.

User interface
^^^^^^^^^^^^^^

The user interface follows :doc:`fix shake <fix_shake>` closely:

* The first three arguments are the convergence tolerance, the maximum
  Newton iteration count, and the statistics-print interval.
* One or more groups of selectors (``b``, ``a``, ``t``, ``m`` lists)
  pick which bonds/angles get constrained.

A bond is constrained when **both** of its atoms are in the fix group
AND at least one of the selectors matches:

* the bond type is in the *b* list, or
* either atom type is in the *t* list, or
* either atom mass (within 0.1 mass units of any value in the *m* list).

An angle of type *aa* listed via *a* contributes a "virtual" A-C distance
constraint when **all three** atoms are in the fix group AND both flanking
bonds (A-B, B-C) are themselves constrained per the criteria above.  The
A-C target distance is computed from the law of cosines using the two
bond equilibrium distances and the angle equilibrium value, identical to
the construction used by :doc:`fix shake <fix_shake>`.

For each bond/angle that is selected, ``fix ilves`` sets the corresponding
``bond_type`` / ``angle_type`` to its negative value so that the configured
:doc:`bond_style <bond_style>` and :doc:`angle_style <angle_style>` skip
the (now-rigid) interaction --- avoiding double-counting of bonded forces.
This mirrors how :doc:`fix shake <fix_shake>` handles the same problem and
is reversed automatically when the fix is destroyed.

Algorithm
^^^^^^^^^

At each time step, after the integrator has updated positions and forces
but before the velocity half-step finalizes, ``fix ilves`` runs the
following sequence:

1. Predict the unconstrained position ``xshake_i = x_i + dt*v_i + dt^2/m_i * f_i``
   for every atom.
2. For each connected constraint cluster, assemble the structurally
   symmetric Jacobian

   .. parsed-literal::

       A_kk = (1/m_a + 1/m_b) * (s_k . r_k)
       A_kl = sigma_k^p * sigma_l^p * (s_k . r_l) / m_p   (k != l, sharing atom p)

   where ``r_k = x[a]-x[b]`` is the bond vector at the start of the step,
   ``s_k = xshake[a]-xshake[b]`` is the current Newton iterate, and
   ``sigma_k^p`` is +1 if atom *p* is the first atom of constraint *k*,
   -1 otherwise.  Both vectors are minimum-image-corrected for PBC.
3. Solve ``A * dlambda = -g`` (where ``g_k = 0.5*(|s_k|^2 - d_k^2)``) by
   dense LU decomposition with partial pivoting.
4. Apply ``dlambda``: ``xshake[a] += dlambda/m_a * r_k``,
   ``xshake[b] -= dlambda/m_b * r_k``, accumulate the Lagrange multipliers.
5. Forward-comm ``xshake`` to neighboring ranks; iterate steps 2-5
   until the global max relative residual falls below ``tol`` (or
   ``iter`` is reached).
6. Convert the accumulated multipliers into per-atom forces
   ``f[a] += lambda/dt^2 * r_k``, ``f[b] -= lambda/dt^2 * r_k`` so that
   the next initial_integrate produces the constrained positions.

After the integrator's final_integrate step, ``end_of_step`` runs the
RATTLE-style velocity projection: solve a similar (symmetric) linear
system to remove the component of relative atom velocity along each
constrained bond direction.

During minimization, the standard integration-based machinery cannot
apply --- there are no velocities or time steps.  Instead, every
constraint *k* is replaced by a strong harmonic potential
``V_k = 0.5 * kbond * (|r_k| - d_k)^2``, contributed to the minimizer
just like a bond.

Output info
^^^^^^^^^^^

This fix computes a *global scalar* equal to the harmonic-restraint
energy accumulated during minimization (zero outside of minimization).
The scalar is accessible via ``f_<ID>`` from compute, variable, dump,
and thermo commands.

With the *store yes* keyword, this fix also exposes a *per-atom array*
with 3 columns containing the constraint-force components ``(fx, fy, fz)``
that were added by the fix on this time step.  The per-atom array is
accessible as ``f_<ID>[1]``, ``f_<ID>[2]``, ``f_<ID>[3]``.

When ``N > 0``, every ``N`` time steps the fix prints a summary line
per constrained bond type and angle type giving the average distance,
its spread (max - min) across constraints, and the count.

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`fix_modify <fix_modify>` *energy*, *virial*, and *peratom*
options are supported by this fix to add the contribution of the
constraint forces to the global potential energy, the global pressure
virial, and the per-atom energy / per-atom virial of the system.

No information about this fix is written to :doc:`binary restart files
<restart>` because the constraint topology is rebuilt from
``atom->bond_atom`` / ``atom->angle_atom`` at the first reneighbor after
restart.  If you redeclare the same ``fix ilves`` command after
``read_restart``, the constraints are restored automatically (the
negated bond/angle types are preserved in the restart file).

This fix is not invoked during :doc:`energy minimization <minimize>` via
the constraint solver; instead, strong harmonic-restraint forces
substitute for the constraints, with spring constant ``kbond`` (see
keyword options).

.. note::

   The implementation gathers the complete bond/angle topology onto
   every MPI rank at init (one ``MPI_Allgatherv`` for bonds and one
   for angles), then on each rank includes every constraint that
   belongs to a cluster intersecting at least one local atom -- even
   constraints between two ghost atoms.  This ensures all participating
   ranks of a cross-rank cluster compute the same Lagrange multipliers
   from synchronized positions, and each rank applies the constraint
   force only to its local atoms.  Both ``newton on`` and ``newton off``
   are supported without loss of accuracy.

   Two solver variants are available via the ``variant`` keyword:

   * ``ilvesf`` (default): symmetric Jacobian using
     :math:`(r_k\cdot r_l)/m_p` in the off-diagonals -- a quasi-Newton
     step that uses reference (start-of-step) bond vectors in the
     coupling terms.  The matrix is genuinely symmetric.
   * ``ilves``: the exact Newton Jacobian using
     :math:`(s_k\cdot r_l)/m_p` in the off-diagonals -- structurally
     but not numerically symmetric.

   Both variants converge to the exact solution of the constraint
   equations; the symmetric variant has slightly cheaper per-iteration
   matrix assembly and would permit a Cholesky solve (currently both
   variants use the same LU with partial pivoting).

Restrictions
""""""""""""

This fix is part of the RIGID package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

The molecular topology (bonds, optionally angles) must be defined; an
:doc:`atom_style <atom_style>` such as *full*, *molecular*, or *bond* is
required.

Only one ``fix ilves`` instance may be defined.  ``fix ilves`` and
:doc:`fix shake <fix_shake>` must not be used together for overlapping
sets of constrained bonds.

Related commands
""""""""""""""""

:doc:`fix shake <fix_shake>`, :doc:`fix restrain <fix_restrain>`,
:doc:`fix rigid <fix_rigid>`, :doc:`fix ehex <fix_ehex>`

Default
"""""""

The keyword defaults are:

* *kbond* = ``1.0e9 * boltz`` (same default as :doc:`fix shake <fix_shake>`)
* *store* = *no*
* *variant* = *ilvesf*

----------

.. _Garcia-Risueno2026:

**(Garcia-Risueno)** P. Garcia-Risueno, A. de la Vega de Leon,
T. Asikainen, L. A. Aguado-Gomez, V. Beltran-Sanchez, and L. Vazquez,
"ILVES: a high-performance constraint solver for molecular dynamics
simulations", J. Chem. Theory Comput. (2026),
DOI: 10.1021/acs.jctc.5c01376

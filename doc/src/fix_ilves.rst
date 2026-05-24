.. index:: fix ilves/local
.. index:: fix ilves/global

fix ilves/local command
=======================

fix ilves/global command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style tol iter N selectors ... [keyword args ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *ilves/local* or *ilves/global*
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
         *ilvesf* = ILVES-F (symmetric quasi-Newton, default; solved with
                    Cholesky for ~15-20% lower solver cost than LU)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves/global 1.0e-4 25 100 b 4 6 8 10 12 14 18 a 31
   fix 2 wat ilves/local  1.0e-6 30 1000 t 1 m 1.008 store yes
   fix 3 all ilves/global 1.0e-5 20 0 b 1 a 1 variant ilves

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

Two variants are provided that differ in how each MPI rank discovers
the constraint topology:

* ``ilves/local`` builds the constraint list from each rank's own
  bond/angle storage only.  Per-rank memory scales with the local atom
  count, not with the total system size.  Every constraint cluster
  must fit within a single subdomain plus the communication cutoff
  (the atoms of the cluster must be either local or available as
  ghosts on the rank that owns the cluster's center).  Suitable for
  water HOH constraints, methyl/methane groups, isolated C-H bonds ---
  any star or 1+1 topology where each cluster is short-ranged in
  space.  Requires ``newton off bond`` (see the
  :doc:`newton <newton>` command).
* ``ilves/global`` gathers the complete bond / angle topology onto
  every MPI rank at init via ``MPI_Allgatherv``.  Per-rank memory
  therefore grows with the total system size, but the variant
  supports constraint clusters of any topology and any spatial extent
  --- including a polymer backbone constrained across many subdomains.
  Use this variant when ``ilves/local`` reports a cluster larger than
  the communication cutoff or a non-star (multi-center) cluster.

User interface
^^^^^^^^^^^^^^

The user interface follows :doc:`fix shake <fix_shake>` closely and is
identical for both variants:

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
but before the velocity half-step finalizes, the fix runs the following
sequence:

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
   dense LU decomposition with partial pivoting (or Cholesky for the
   *ilvesf* variant -- see below).
4. Apply ``dlambda``: ``xshake[a] += dlambda/m_a * r_k``,
   ``xshake[b] -= dlambda/m_b * r_k``, accumulate the Lagrange multipliers.
5. ``ilves/global``: forward-comm ``xshake`` to neighboring ranks after each
   iteration so the redundant cluster solves on multiple participating
   ranks stay in sync; iterate steps 2-5 until the global max relative
   residual falls below ``tol``.

   ``ilves/local``: each cluster is owned by a single rank (the rank
   holding the cluster's star center, or the lower-tag endpoint of a
   1+1 pair).  The owner runs Newton iteration using the local plus
   ghost-cell view; no per-iteration ``forward_comm`` is needed.
6. Convert the accumulated multipliers into per-atom forces
   ``f[a] += lambda/dt^2 * r_k``, ``f[b] -= lambda/dt^2 * r_k`` so that
   the next initial_integrate produces the constrained positions.
   ``ilves/local`` runs ``reverse_comm`` to deliver the constraint
   contributions on ghost atoms back to their owners' local atom.

After the integrator's ``final_integrate`` step, ``end_of_step`` runs the
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
restart.  If you redeclare the same ``fix ilves/local`` or
``fix ilves/global`` command after ``read_restart``, the constraints are
restored automatically (the negated bond/angle types are preserved in the
restart file).

This fix is not invoked during :doc:`energy minimization <minimize>` via
the constraint solver; instead, strong harmonic-restraint forces
substitute for the constraints, with spring constant ``kbond`` (see
keyword options).

Solver variants
^^^^^^^^^^^^^^^

Two solver variants are available via the ``variant`` keyword (both
work with either ``ilves/local`` or ``ilves/global``):

* ``ilvesf`` (default): symmetric Jacobian using
  :math:`(r_k\cdot r_l)/m_p` in the off-diagonals -- a quasi-Newton
  step that uses reference (start-of-step) bond vectors in the
  coupling terms.  The matrix is genuinely symmetric.
* ``ilves``: the exact Newton Jacobian using
  :math:`(s_k\cdot r_l)/m_p` in the off-diagonals -- structurally
  but not numerically symmetric.

Both variants converge to the exact solution of the constraint
equations.  The symmetric ``ilvesf`` variant is solved with an
in-place Cholesky factorization (``A = L L^T``); on a 2004-atom
peptide-in-water test with all bonds + the water HOH angle
constrained, the Cholesky path is ~14% faster than LU on the same
matrix in serial and ~21% faster on 4 MPI ranks.  If Cholesky
encounters a non-positive pivot (degenerate cluster, linearly
dependent constraints) the solver re-assembles the matrix and falls
back to LU with partial pivoting -- the per-step ``output_every``
stats line reports the fallback count when ``N > 0``.  The
asymmetric ``ilves`` variant uses LU directly.  The velocity
projection in ``end_of_step`` is also Cholesky-based (the velocity
matrix is always symmetric regardless of variant).

Choosing between ``ilves/local`` and ``ilves/global``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Prefer ``ilves/local`` whenever the constraint topology consists of
small clusters that each fit within a single subdomain plus the
communication cutoff:

* Water HOH (2 O-H bonds + 1 H-O-H angle = 3-atom star).
* Methyl / methane groups (C with 3 or 4 H bonds = 4/5-atom star).
* Isolated C-H / N-H / O-H bond constraints (1+1 pairs).

For these cases, ``ilves/local``:

* Avoids the ``MPI_Allgatherv`` of the global topology at init.
* Stores no ``natoms``-sized cluster-tag table.
* Uses single-rank cluster ownership, eliminating redundant per-cluster
  Newton solves on neighboring ranks.

Choose ``ilves/global`` when:

* A cluster is genuinely large -- e.g., constraining every backbone
  C-C bond of a polymer or protein, so the cluster's atoms span
  multiple subdomains.
* You want to constrain non-star cluster topologies (rare; most
  hydrogen-only constraints form stars).
* Convenience: the global variant does not require
  ``newton off bond`` and handles arbitrary cluster topologies
  without preconditions on the simulation setup.

.. note::

   Each rank in ``ilves/global`` gathers the complete bond / angle
   topology and stores it in replicated tables.  In addition, a
   cluster-tag lookup table of size ``atom->natoms+1`` is allocated.
   For very large simulations (tens of millions of atoms or more)
   this replication can exceed available memory at high rank counts.
   Prefer ``ilves/local`` or :doc:`fix shake <fix_shake>` in that
   regime unless the cluster topology really requires the global
   variant.

Restrictions
""""""""""""

These fixes are part of the RIGID package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

The molecular topology (bonds, optionally angles) must be defined; an
:doc:`atom_style <atom_style>` such as *full*, *molecular*, or *bond* is
required.

Only one ``fix ilves/*`` instance may be defined at a time
(``ilves/local`` and ``ilves/global`` may not be combined).
``fix ilves/*`` and :doc:`fix shake <fix_shake>` must not be used
together for overlapping sets of constrained bonds.

``fix ilves/local`` additionally requires:

* ``newton off bond`` (use the :doc:`newton <newton>` command).
* Every constraint cluster must fit within this rank's local subdomain
  plus the communication cutoff.  Cluster atoms unreachable as either
  a local atom or a ghost trigger a hard error at the first neighbor
  build; switch to ``ilves/global`` or enlarge the cutoff via
  ``comm_modify cutoff <value>``.
* Cluster topology must be a star (one central atom with k bonds to
  leaves, k >= 1) or a 1+1 pair (single bond).  Multi-center clusters
  (e.g. constraining an interior C-C bond of a polymer chain with
  multiple branching points) are rejected with an error message
  recommending ``ilves/global``.

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

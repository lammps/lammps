.. index:: fix ilves

fix ilves command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ilves tol iter N selectors ... [keyword args ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
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
       *variant* value = *full* or *fast*
         *fast* = symmetric quasi-Newton with banded Cholesky (default)
         *full* = exact-Newton (asymmetric) with LU decomposition
       *linearangle* value = theta_deg [Lmin]
         theta_deg = threshold in degrees
         Lmin      = minimum \|B-M\| target length

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves 1.0e-4 25 100 b 4 6 8 10 12 14 18 a 31
   fix 2 wat ilves 1.0e-6 30 1000 t 1 m 1.008 store yes
   fix 3 all ilves 1.0e-5 20 0 b 1 a 1 variant full
   fix 4 all ilves 1.0e-6 30 100 b 1 a 1 variant fast linearangle 170 0.05

Description
"""""""""""

.. versionadded:: TBD

Apply bond-length (and optionally end-to-end angle "virtual bond")
constraints using the ILVES algorithm of :ref:`(Garcia-Risueno)
<Garcia-Risueno2026>`.  ILVES enforces holonomic distance constraints
using Newton's method on a sparse system of nonlinear equations.  Unlike
:doc:`fix shake <fix_shake>`, ILVES handles arbitrarily large connected
constraint clusters --- including the full set of C-C backbone bonds of
a long polymer or protein chain --- in a single solve.

The fix gathers the complete bond / angle topology onto every MPI rank
at init via ``MPI_Allgatherv`` and builds a per-rank constraint list
that includes every cluster intersecting the rank's local atoms.  Each
participating rank solves intersecting clusters redundantly so that the
same Lagrange multipliers are computed from synchronized positions; the
constraint force is applied locally to each rank's owned atoms only.
Per-rank topology storage is sparse --- only tags involved in at least
one selected bond or angle are kept, so the memory cost scales with the
constrained subset rather than with ``natoms``.

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
but before the velocity half-step finalizes, the fix runs the following
sequence:

1. Predict the unconstrained position ``xshake_i = x_i + dt*v_i + dt^2/m_i * f_i``
   for every atom.
2. For each connected constraint cluster, assemble the Jacobian and
   compute its banded Cholesky factor (``A = L L^T``).  The clusters
   are pre-permuted by reverse Cuthill-McKee at constraint-list build
   time so that the band is as narrow as possible.  Storage cost is
   ``O(n_c * bw_c)`` per cluster rather than ``O(n_c^2)`` for the dense
   form.
3. Solve ``A * dlambda = -g`` (where ``g_k = 0.5*(|s_k|^2 - d_k^2)``) by
   forward and back substitution against the cached Cholesky factor.
4. Apply ``dlambda`` to ``xshake`` and accumulate the Lagrange
   multipliers.
5. ``forward_comm`` ``xshake`` so the redundantly-solving ranks stay in
   sync, and iterate steps 2-5 until the global max relative residual
   falls below ``tol``.  Within a step, only the right-hand side ``g``
   and the Newton iterate ``s_k = xshake[a] - xshake[b]`` change; the
   Jacobian is step-constant for ``variant fast`` and is factored only
   once per step per cluster.
6. Convert the accumulated multipliers into per-atom forces
   ``f[a] += lambda/dt^2 * r_k``, ``f[b] -= lambda/dt^2 * r_k`` so that
   the next ``initial_integrate`` produces the constrained positions.

After ``final_integrate``, ``end_of_step`` runs the RATTLE-style
velocity projection: solve a similar (symmetric) linear system to remove
the component of relative atom velocity along each constrained bond
direction.

During minimization, the standard integration-based machinery cannot
apply --- there are no velocities or time steps.  Instead, every
constraint *k* is replaced by a strong harmonic potential
``V_k = 0.5 * kbond * (|r_k| - d_k)^2``, contributed to the minimizer
just like a bond.

Solver variants
^^^^^^^^^^^^^^^

Two solver variants are available via the ``variant`` keyword:

* ``fast`` (default): the structurally-symmetric Jacobian using
  :math:`(r_k\cdot r_l)/m_p` in the off-diagonals -- a quasi-Newton
  step that uses reference (start-of-step) bond vectors in the
  coupling terms.  The matrix is genuinely symmetric and is solved
  with an in-place **banded Cholesky** factorization that is cached
  and reused across all Newton iterations within a step.  Combined
  with the RCM band reduction, this is up to ~5x faster than the
  ``full`` variant on the all-backbone polymer / protein benchmark.
  If the Cholesky path encounters a non-positive pivot (a degenerate
  cluster --- typically a co-linear angle), the solver falls back to
  LU with partial pivoting for that cluster; the ``output_every``
  stats line reports the fallback count.

* ``full``: the exact Newton Jacobian using
  :math:`(s_k\cdot r_l)/m_p` in the off-diagonals -- structurally
  but not numerically symmetric (changes each iteration), solved
  with LU on the dense matrix.  Use when ``fast`` repeatedly falls
  back to LU on the same clusters (degenerate constraint topology
  not handled by Cholesky) and you want to skip the factor-attempt
  overhead.

Both variants converge to the same solution of the constraint
equations.

The velocity projection in ``end_of_step`` always uses dense Cholesky
(the velocity matrix is genuinely symmetric regardless of variant).

Near-linear angles
^^^^^^^^^^^^^^^^^^

For an angle A-B-C with equilibrium :math:`\theta_0` close to
180 degrees (common in coarse-grain polymer backbones), constraining
the three legs :math:`\{|AB|, |BC|, |AC|\}` makes the constraint
Jacobian rank-deficient: the triangle inequality saturates
(:math:`|AC| = |AB| + |BC|`), so the three constraints become linearly
dependent at exactly 180 degrees and very ill-conditioned nearby.

To handle this regime, the ``linearangle`` keyword switches the
constraint topology for affected angle types.  An angle type is
classified as near-linear when its equilibrium angle (recovered from
the bond and A-C virtual distances via the law of cosines) is at or
above the user-supplied threshold.  For those angles:

* The real bond between B and the **higher-tag endpoint** of {A, C}
  is dropped.
* The A-C virtual bond is kept unchanged.
* A 3-atom **B-M virtual constraint** is added with target

  .. math::

     |B - M|^2 = \tfrac{1}{2}(|AB|^2 + |BC|^2) - \tfrac{1}{4}|AC|^2

  where :math:`M = (A + C) / 2` is the midpoint of the A-C virtual
  bond.  The retained set keeps :math:`|B-C|` within solver tolerance
  via the triangle geometry while remaining well-conditioned at
  :math:`\theta` near 180 degrees.

For an exactly-180-degree symmetric angle (:math:`|AB| = |BC|`,
:math:`\theta_0 = 180^\circ`) the natural :math:`|B-M|` target is
zero and the constraint Jacobian row vanishes.  The ``Lmin`` argument
to ``linearangle`` clamps the target up from zero in that limit; the
default value of 0.01 length units bends the angle by less than
1 degree off 180 degrees for typical bond lengths, well below thermal
fluctuation amplitudes.

Set ``linearangle 180`` to disable near-linear handling entirely (the
default threshold is 165 degrees).

.. note::

   The ``linearangle`` feature is a **workaround**, not a full
   solution.  At :math:`\theta` near 180 degrees the constraint
   manifold is geometrically nearly singular and any iterative
   distance-constraint solver (SHAKE, RATTLE, P-SHAKE, ILVES) is
   ill-conditioned there.  ``linearangle`` widens the regime in which
   the iterative solver still converges, but it does so by
   substituting an approximate constraint -- the :math:`|B-M|` clamp
   biases the equilibrium angle slightly off 180 degrees (by an angle
   that increases with ``Lmin``) and only constrains the median
   length, not the bond :math:`|B-C|` directly (which is determined
   only indirectly through the other three legs).

   This feature is most useful for large systems with a **small
   fraction** of near-linear angles -- e.g. a coarse-grain polymer or
   protein system where most angles are well-behaved (tetrahedral,
   sp2 trigonal, water-like) and only a handful of backbone or
   special-purpose angles approach 180 degrees.  In that regime the
   small error from the ``Lmin`` clamp on a few constraints is
   negligible compared to the thermal fluctuations that would
   otherwise unconstrain the rest of the system.

   For the simulation of intrinsically rigid linear molecules
   (e.g. CO\ :sub:`2`, CS\ :sub:`2`, HCN, C\ :sub:`2`\ H\ :sub:`2`),
   :doc:`fix rigid/small <fix_rigid>` is the recommended approach.
   It treats each molecule as a rigid body and has no iterative
   convergence problems at the 180 degree singularity.  The ILVES
   ``linearangle`` workaround should not be used as a substitute for
   ``fix rigid/small`` in those cases -- the angle bias from the
   ``Lmin`` clamp and the larger iteration counts both make ILVES the
   wrong tool for that workload.

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
its spread (max - min) across constraints, and the count.  3-atom B-M
virtual constraints are not included in the per-type angle stats
because their length (the median to the A-C side) is geometrically
distinct from the A-C virtual length reported for the angle.

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

Restrictions
""""""""""""

This fix is part of the RIGID package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

The molecular topology (bonds, optionally angles) must be defined; an
:doc:`atom_style <atom_style>` such as *full*, *molecular*, or *bond* is
required.

Only one ``fix ilves`` instance may be defined at a time.  ``fix ilves``
and :doc:`fix shake <fix_shake>` must not be used together for
overlapping sets of constrained bonds.

For exactly-180-degree symmetric angle constraints (e.g. rigid linear
triatomics like CO2 with bond_length(O-C) = bond_length(C-O)) the
Jacobian is irreducibly rank-deficient regardless of the constraint
topology.  ``fix ilves`` runs in that regime with the ``Lmin`` clamp
adding a small constraint bias, but :doc:`fix rigid <fix_rigid>` is
the recommended approach.

Memory requirements
^^^^^^^^^^^^^^^^^^^

The replicated bond / angle topology is gathered onto every MPI rank
at init.  The full table is stored on **every** rank, so the
per-rank memory cost does **not** decrease as the rank count grows;
on the contrary, the aggregate (Nranks * this) is what the host
machine must accommodate.

At init time the fix prints a line such as::

   Fix ilves: gathered global topology with NB bonds, NA selected angles, NT involved atoms
   Fix ilves: replicated topology storage ~ X bytes/rank (X.XX MB)

where the printed estimate is the sum of:

* bond table (``gb_a``, ``gb_b``, ``gb_type``): ``2*NB*sizeof(tagint) + NB*sizeof(int)``
* angle table (``ga1/2/3``, ``ga_type``): ``3*NA*sizeof(tagint) + NA*sizeof(int)``
* the ``tag_cluster`` map: roughly ``48*NT`` bytes (per-entry overhead
  for ``std::unordered_map`` -- node + bucket + allocator overhead,
  measured against libstdc++; libc++ and MSVC implementations are
  similar to within ~30%)

For the default ``-DLAMMPS_SMALLBIG`` build (32-bit ``tagint``,
32-bit ``int``) this works out to:

* ~12 bytes per selected bond
* ~16 bytes per selected angle
* ~48 bytes per atom involved in some selected bond / angle

So a system with one constraint per bonded hydrogen (~3 per heavy
atom typical for proteins) and 100M atoms would store roughly:

* bonds:   3 * 100M * 12 B   = 3.6 GB / rank
* angles:  1 * 100M * 16 B   = 1.6 GB / rank
* tags:    4 * 100M * 48 B   = 19.2 GB / rank
* total: ~24 GB / rank

Far too large for typical compute nodes.  An all-bond-constrained
system at that scale would also blow up.  Rules of thumb:

* **< 10M selected bonds total**: fits comfortably in a few hundred
  MB / rank; no concerns.
* **10M -- 100M selected bonds**: the storage is hundreds of MB to a
  few GB per rank.  Check ``Fix ilves: replicated topology storage``
  against your host's per-rank memory budget before committing to a
  long run.
* **> 100M selected bonds**: ``fix ilves`` in its current
  implementation is likely not viable.  If the constraint clusters
  are small and strictly local (water HOH, methyl groups, isolated
  C-H bonds), use :doc:`fix shake <fix_shake>` -- its per-rank
  memory scales with the local atom count, not the global atom count.
  If the clusters span the whole simulation cell (e.g. a polymer
  backbone connected across all subdomains), the current
  implementation has no escape hatch and is the wrong tool at this
  scale.

The per-rank ``memory_usage`` reported by ``info memory`` and
``thermo_style custom ... memory`` includes the replicated topology
plus the per-rank constraint-list arrays (rebuilt each reneighbor)
and the per-cluster Cholesky factor cache (``variant fast`` only).
For dynamic-load-balanced or non-uniform systems, the cluster
factor cache size on each rank depends on how many cluster atoms
that rank touches; consult ``info fixes`` for the breakdown.

Related commands
""""""""""""""""

:doc:`fix shake <fix_shake>`, :doc:`fix rattle <fix_shake>`,
:doc:`fix restrain <fix_restrain>`, :doc:`fix rigid <fix_rigid>`

Default
"""""""

The keyword defaults are:

* *kbond* = :math:`1.0^9 k_B` (same default as :doc:`fix shake <fix_shake>`)
* *store* = *no*
* *variant* = *fast*
* *linearangle* = 165 degrees with Lmin = 1/30 length units

----------

.. _Garcia-Risueno2026:

**(Garcia-Risueno)** P. Garcia-Risueno, A. de la Vega de Leon,
T. Asikainen, L. A. Aguado-Gomez, V. Beltran-Sanchez, and L. Vazquez,
"ILVES: a high-performance constraint solver for molecular dynamics
simulations", J. Chem. Theory Comput. (2026),
DOI: 10.1021/acs.jctc.5c01376

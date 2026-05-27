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
       *linearangle* values = theta_deg [K]
         theta_deg = threshold in degrees above which an angle constraint
                     is skipped
         K         = optional stiff angle force constant (energy units).
                     When K > 0, near-linear angles are still tagged
                     as "no AC constraint", but additionally the
                     ``angle_type`` slot is negated (so the user's
                     :doc:`angle_style <angle_style>` skips them) and
                     ``fix ilves`` applies its own stiff angle force
                     E = K * (1 + cos(theta)) -- the canonical
                     "cosine" form with no 1/sin singularity at
                     theta = 180.  K = 0 (default) leaves the
                     angle_style force in place.
       *warmstart* value = *yes* or *no*
         *yes* initializes the Lagrange multipliers of each Newton solve
         from the converged values of the previous timestep (with the
         corresponding ``xshake`` position correction).  Reduces Newton
         iteration counts on multi-rank runs whose convergence is slowed
         by additive-Schwarz iteration between subdomains.  Discarded
         automatically on reneighbor (constraint indexing is invalidated).
         Default is *no*, which preserves bit-for-bit equivalence with
         the original cold-start path used by the unit-test references;
         enabling produces ULP-level (~5e-11) trajectory differences
         that compound chaotically over many steps.
       *relax* value = omega
         omega = Newton-update damping factor in (0, 1], default 1.0.
         Each Newton iteration's ``dlambda`` is scaled by omega before
         being applied.  Values below 1 damp Schwarz-iteration
         oscillations between ranks at the cost of slower per-iteration
         progress; in practice on well-conditioned clusters
         ``relax = 1.0`` converges fastest.

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves 1.0e-4 25 100 b 4 6 8 10 12 14 18 a 31
   fix 2 wat ilves 1.0e-6 30 1000 t 1 m 1.008 store yes
   fix 3 all ilves 1.0e-5 20 0 b 1 a 1 variant full
   fix 4 all ilves 1.0e-6 30 100 b 1 a 1 variant fast linearangle 170

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
Newton iterations on such an ill-conditioned system produce large
Lagrange multipliers that translate into excessive forces, often
showing up as ``Atoms moved too far for minimum image`` errors a few
steps later.

To avoid this regime, ``fix ilves`` classifies every selected angle
type whose equilibrium :math:`\theta_0` (recovered from the bond and
A-C virtual distances via the law of cosines) is at or above the
``linearangle`` threshold as **near-linear**.  For those angle types:

* The angle's A-C virtual bond constraint is **not added** to the
  constraint list.
* The angle's ``angle_type`` slot is **not negated**, so the standard
  :doc:`angle_style <angle_style>` force-field term continues to act
  on the angle and keeps it close to :math:`\theta_0`.
* The two flanking bonds (A-B, B-C) are constrained as normal,
  provided their bond types are themselves selected via the ``b``
  selector.

The default threshold is 165 degrees; the smallest allowed value is
150 degrees.  Set *linearangle* to 180 to disable the bailout
entirely, so that the A-C constraint is added for every selected
angle including the ill-conditioned ones (likely to produce
``max_iter`` warnings or atom-ejection failures).

For near-linear angles you can also opt into a stiff angle restraint
applied by ``fix ilves`` itself, in place of the user's
:doc:`angle_style <angle_style>`.  Pass an optional ``K`` after the
threshold:

.. code-block:: LAMMPS

   fix cstr all ilves 1.0e-6 30 0 b 1 2 a 1 variant fast linearangle 165 5.0

With ``K > 0``, near-linear angle types have their ``angle_type`` slot
negated (so the user's angle_style skips them, the same way it does
for ordinary constrained angles), and ``fix ilves`` applies
:math:`E = K \,(1 + \cos\theta)` -- the standard "cosine" angle form,
identical to LAMMPS :doc:`angle_style cosine <angle_cosine>`.  This
form has its minimum at :math:`\theta = 180^\circ` and a non-vanishing
restoring curvature there, with no :math:`1/\sin(\theta)` singularity
in the force.  Useful when the rest of the system uses harmonic angles
(which would blow up at :math:`\theta = 180^\circ`).

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

Restrictions
""""""""""""

This fix is part of the RIGID package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

The molecular topology (bonds, optionally angles) must be defined; an
:doc:`atom_style <atom_style>` such as *full*, *molecular*, *bond*, or
angle (or a hybrid atom style including them) is required.

Only one ``fix ilves`` instance may be defined at a time.  ``fix ilves``
and :doc:`fix shake <fix_shake>` must not be used together for
overlapping sets of atoms participating in a constraint.

``fix ilves`` works under both :doc:`newton on bond <newton>` (the
LAMMPS default) and ``newton off bond``.  Under newton off the
constraint solver looks up bond types directly from local bond
storage; under newton on it uses an MPI-reduced per-angle-type cache
(``angle_btype1`` / ``angle_btype2``) built at init for the
cross-rank lookup.  The optional stiff angle force
(``linearangle <theta> <K>`` with ``K > 0``) works under both newton
modes -- under newton on its per-angle force trio is accumulated into
a fix-owned per-atom buffer at the middle-atom owner and reverse_comm
sends ghost-side contributions back to owner ranks before being added
to ``atom->f``.

``fix ilves`` does not support dynamic topology.  The bond and angle
tables are read from local atom storage at init and assumed to remain
valid for the duration of the run.  Fixes or commands that mutate the molecular
topology *during* a run -- for example :doc:`fix bond/create
<fix_bond_create>`, :doc:`fix bond/break <fix_bond_break>`, :doc:`fix
bond/react <fix_bond_react>`, :doc:`fix deposit <fix_deposit>` must not
add, remove, or change any constrained atoms.  Any detected change
aborts the run with an error message naming the most likely cause.
Note, that this test does not detect rearrangements that keep the counts
identical.

Angle types whose equilibrium :math:`\theta_0` is at or above the
``linearangle`` threshold (default 165 degrees) are silently skipped
when building the constraint list.

Memory requirements
^^^^^^^^^^^^^^^^^^^

``fix ilves`` stores topology in **distributed** form: each MPI rank
builds its constraint list from its own local atom storage
(``atom->bond_*``, ``atom->angle_*``) plus ghost atoms in the standard
communication shell.  No replicated global table is held, so per-rank
memory shrinks roughly linearly with rank count for partial-constraint
workloads (hydrogens-only, methyl, water) and approximately linearly
for cluster-spanning workloads as long as clusters fit within the
ghost shell.

A reference benchmark on the rhodopsin all-bonds case (~31k selected
constraints total, single dominant cluster) shows per-rank memory:

============  =====================  ======================
MPI ranks     ``info memory`` MB     constraints/proc
============  =====================  ======================
1             193                    31955
2             137                    16200
4              69                     8197
8              61                     4172
============  =====================  ======================

The per-rank ``memory_usage`` reported by ``info memory`` and
``thermo_style custom ... memory`` includes the per-rank constraint-list
arrays (rebuilt each reneighbor) and the per-cluster Cholesky factor
cache (``variant fast`` only).  Both scale with this rank's own atoms.
Consult ``info fixes`` for the per-fix breakdown.

For clusters whose extent exceeds the standard LAMMPS ghost cutoff
(e.g. a fully-bond-constrained polymer chain spanning a whole
subdomain), the constraint solver still works -- each rank handles
its local piece of the cluster, with Newton iteration converging via
ghost-atom updates between iterations.  Newton iteration count may
increase modestly for large spanning clusters compared to the prior
redundant-solve model.

Related commands
""""""""""""""""

:doc:`fix shake <fix_shake>`, :doc:`fix rattle <fix_shake>`,
:doc:`fix restrain <fix_restrain>`, :doc:`fix rigid <fix_rigid>`

Default
"""""""

The keyword defaults are:

* *kbond* = 1.0e9*k_B (same default as :doc:`fix shake <fix_shake>`)
* *store* = *no*
* *variant* = *fast*
* *linearangle* = 165 degrees

----------

.. _Garcia-Risueno2026:

**(Garcia-Risueno)** P. Garcia-Risueno, A. de la Vega de Leon,
T. Asikainen, L. A. Aguado-Gomez, V. Beltran-Sanchez, and L. Vazquez,
"ILVES: a high-performance constraint solver for molecular dynamics
simulations", J. Chem. Theory Comput. (2026),
DOI: 10.1021/acs.jctc.5c01376

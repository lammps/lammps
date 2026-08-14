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
* selectors = one or more of *b*, *a*, *t*, *m*, each followed by one or more numeric values

  .. parsed-literal::

       *b* values = one or more bond types (may use type labels)
       *a* values = one or more angle types (may use type labels)
       *t* values = one or more atom types (may use type labels)
       *m* values = one or more atom masses

* zero or more keyword/value pairs may be appended
* keyword = *mode* or *linearangle* or *kbond* or *store*

  .. parsed-literal::

       *mode* value = *converge* or *fixed*
         *converge* = iterate until the tolerance is met (default)
         *fixed* = always perform *iter* Newton iterations
       *linearangle* values = style threshold
         style = *error* or *skip* or *restrain*
           *error* = stop if a selected angle type is near-linear (default)
           *skip* = do not constrain near-linear angle types
           *restrain* = replace the near-linear A-C constraint with a stiff harmonic restraint
         threshold = equilibrium angle (in degrees) at or above which an angle type is near-linear
       *kbond* value = dynamics force constant of the harmonic-bond substitute used
         for the *linearangle restrain* mode (the energy-minimization substitute is a
         fixed factor stiffer); in the active unit system
       *store* value = *yes* or *no*
         *yes* exposes the per-atom constraint forces via a per-atom array

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves 1.0e-6 25 0 b 4 6 8 10
   fix 2 wat ilves 1.0e-8 25 1000 b 18 a 31
   fix 3 all ilves 1.0e-6 25 0 b 1 mode fixed
   fix 4 sol ilves 1.0e-8 25 0 t 1 m 1.008 store yes
   fix 5 co2 ilves 1.0e-8 25 0 b 1 a 1 linearangle restrain 170 kbond 2000

Description
"""""""""""

.. versionadded:: 4Jul2026

Apply bond-length and angle constraints using the ILVES algorithm of
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
<fix_shake>`.

User interface
^^^^^^^^^^^^^^

The user interface of *fix ilves* follows that of :doc:`fix shake
<fix_shake>`.  The three arguments after the *ilves* keyword are the
tolerance :math:`\frac{|g_k|}{d_k^2}`
(where :math:`g_k = \frac{1}{2}(|s_k|^2 - d_k^2)` and :math:`d_k` is
the target length of constraint *k*), the maximum number of Newton
iterations per step, and the frequency of statistics output
(0 turns it off).

Then one or more groups of selectors (``b``, ``a``, ``t``, ``m`` lists)
pick which bonds and angles get constrained.  This selection is
*identical* to :doc:`fix shake <fix_shake>`.  A bond is constrained when
*both* of its atoms are in the fix group AND at least one of the
selectors matches:

* the bond type is in the *b* list, or
* either atom type is in the *t* list (i.e. all bonds connected to an
  atom of a listed type are constrained), or
* either atom mass is within a fudge factor of MASSDELTA (0.1 mass
  units, defined in ``src/RIGID/fix_ilves.cpp``) of a value in the *m*
  list (i.e. all bonds connected to an atom of a listed mass are
  constrained).

The types may be given as type labels *only* if there is no atom, bond,
or angle type label named *b*, *a*, *t*, or *m* defined in the
simulation.  If that is the case, type labels cannot be used as
constraint type index with these two fixes, because the type labels
would be incorrectly treated as a new type of constraint instead.  Thus,
LAMMPS will print a warning and type label handling is disabled and
numeric types must be used.

An angle whose type is in the *a* list contributes a "virtual bond"
constraint on the distance between its two outer atoms (A and C of an
A-B-C angle), which together with the two constrained legs makes the
angle rigid.  The angle is constrained when all three of its atoms are
in the fix group, its type is selected, and *both* flanking bonds (A-B
and B-C) are themselves constrained.  The A-C target distance is
computed from the two bond equilibrium lengths and the angle equilibrium
value via the law of cosines, identical to :doc:`fix shake <fix_shake>`.

For each bond or angle that is selected, *fix ilves* sets the
corresponding ``bond_type`` or ``angle_type`` to its negative value so
that the configured :doc:`bond_style <bond_style>` / :doc:`angle_style
<angle_style>` skips the (now rigid) interaction and thus avoids
double-counting of bonded forces.  This mirrors how :doc:`fix shake
<fix_shake>` handles the same problem and is reversed automatically when
the fix is deleted.

Unlike :doc:`fix shake <fix_shake>`, which only supports small isolated
clusters, the legs and the A-C virtual bond of an angle become part of
the same connected constraint cluster that the ILVES solver handles
directly, so angles sharing atoms (or angles within larger constrained
networks) need no special casing.

Near-linear angles
^^^^^^^^^^^^^^^^^^

As the equilibrium angle approaches 180 degrees the A-C virtual bond
becomes rank-deficient: the A-C distance is nearly stationary with
respect to the angle, so the constraint is ill-conditioned and the solve
degrades and eventually fails.  Angle types whose equilibrium angle is
at or above the *linearangle* threshold (default 175 degrees) are
therefore treated specially according to the *linearangle* style:

* *error* (default): stop with an error identifying the offending angle type(s).
* *skip*: do not constrain those angle types at all.  Their two flanking bonds
  are still held rigid, but the angle itself is left to the configured
  :doc:`angle_style <angle_style>`, which is the recommended choice when the
  force field provides a (well-behaved) angle term for the near-linear angle.
* *restrain*: replace the rank-deficient A-C constraint with a stiff harmonic
  restraint on the A-C distance, :math:`E = k\,(r_{AC} - d_{AC})^2`, with force
  constant *k* set by the *kbond* keyword (default 1000 in the active unit
  system).  Because the A-C distance is insensitive to the angle near 180
  degrees, this is a soft control on the angle itself; *skip* gives tighter
  angle behavior when a force-field angle term is available.

Unlike a constrained angle, a *skip* or *restrain* angle keeps its
bending degree of freedom: *skip* leaves it free, and the *restrain*
substitute is a soft potential whose A-C coordinate still vibrates and
carries kinetic energy.  Only the two rigid legs remove degrees of
freedom in these modes, and the reported temperature accounts for this
automatically.  For *restrain*, choose *kbond* soft enough that this
vibration is resolved by the timestep; an excessively stiff value is
neither stable nor properly thermostatted.

The restraint energy of the *restrain* substitute is available as the
global scalar of this fix and, with :doc:`fix_modify <fix_modify>`
*energy yes*, is added to the potential energy and the thermodynamic
output.

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

As with :doc:`fix shake <fix_shake>`, the constraints are enforced on
the positions only; the velocities are projected onto the constraint
manifold once at the start of a run (the analogue of SHAKE's initial
velocity correction) but not every step.  ``fix ilves`` therefore has no
special requirement on its placement relative to other time-integration
fixes, unlike a full RATTLE integrator.

The constraint topology is distributed: each MPI rank builds its
constraint list from its own local bond storage plus the ghost atoms in
the standard communication shell.  The Newton iteration is global ---
its convergence test is the maximum relative bond-length violation
reduced over all ranks --- and the cross-rank coupling is resolved by
communicating the predicted positions to the ghosts and reverse-summing
the per-atom corrections to their owners each iteration.  The results
are therefore independent of the number of MPI ranks (to the solver
tolerance).

By default the Newton iteration runs until the global maximum relative
bond-length violation falls below *tol* (or *iter* iterations are
reached).  With ``mode fixed`` the solver instead performs exactly
*iter* iterations every step.  This removes the per-iteration global
reduction used by the convergence test and can improve parallel
efficiency, at the cost of not guaranteeing the tolerance.  Because
ILVES converges quadratically, a small fixed count is usually
sufficient; run first with the default ``mode converge`` and a nonzero
statistics interval *N* to see how many iterations are actually needed
(the maximum is reported with the statistics) before choosing a fixed
count.

Output info
^^^^^^^^^^^

When ``N > 0``, every ``N`` time steps the fix prints constraint
statistics in the same layout as :doc:`fix shake <fix_shake>`, so the
two can be parsed interchangeably.  Each line gives the bond or angle
type, the average value (bond length for ``Bond:`` lines, bend angle in
degrees for ``Angle:`` lines), the spread (max - min) of that value
across all constraints of the type, and the count.  The header line
additionally reports the largest number of Newton iterations used since
the previous output, which is specific to ILVES and is useful for
choosing a fixed iteration count (see ``mode fixed`` above).

With the ``store yes`` keyword, this fix exposes a *per-atom array* with
3 columns containing the constraint-force components ``(fx, fy, fz)``
added by the fix on the current time step, accessible as ``f_<ID>[1]``,
``f_<ID>[2]``, and ``f_<ID>[3]``.

This fix computes a global scalar, the potential energy of the
*linearangle restrain* substitute (zero unless that mode is active),
accessible as ``f_<ID>``.

The constraint forces contribute to the global pressure virial by
default, as for :doc:`fix shake <fix_shake>`; use :doc:`fix_modify
<fix_modify>` *virial no* to exclude their contribution.  The
contribution is computed from the converged constraint multipliers and
the start-of-step bond vectors, so it is reproduced exactly after a
:doc:`restart <read_restart>` and does not depend on whether the masses
are defined per atom type or per atom.  The *linearangle restrain* and
minimization restraint energy is likewise included in the potential
energy by default.

Choosing between *fix shake* and *fix ilves*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a rule of thumb, use :doc:`fix shake <fix_shake>` wherever you can
without running into its limitation of constraining only small isolated
clusters, and use *fix ilves* elsewhere.  Both fixes constrain the same
bonds and angles to the same rigid geometry with the same tolerance.

:doc:`fix shake <fix_shake>` is limited to small, isolated clusters (a
central atom bonded to at most three others, or two bonds and an angle)
and the code is optimized for this use case by minimizing communication
and solving the constraints for each cluster independently.

Using *fix ilves* instead of *fix shake* avoids the small cluster
limitation and can, for example, handle chains of constrained bonds, but
that requires a different kind of constraint solver -- one that needs
communication between MPI ranks and, in parallel, more iterations and
thus is not as efficient.

If you have a system that has both cases, e.g. a peptide with a
constrained backbone and a rigid water solvent, you can consider
using both fixes at the same time. That is, use *fix ilves* for
the peptide and *fix shake* for the water.

------------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`; the constraint topology is rebuilt from ``atom->bond_atom``
at the first reneighbor.  Redeclare the same ``fix ilves`` command after
:doc:`read_restart <read_restart>` to restore the constraints (the
negated bond types are preserved in the restart file).

The :doc:`fix_modify <fix_modify>` *virial* and *energy* options are
supported and both are enabled by default: the constraint contribution
to the global pressure is included in the virial (use *fix_modify virial
no* to exclude it), and the restraint and minimization energy is
included in the potential energy (use *fix_modify energy no* to exclude
it).

During :doc:`energy minimization <minimize>` there is no time
integration, so the holonomic constraints cannot be enforced as during
dynamics.  Instead each constrained bond (and angle A-C virtual bond) is
replaced by a stiff harmonic bond :math:`E = k_\text{bond} (r - d)^2`
whose energy and forces are added to the minimization, exactly as
:doc:`fix shake <fix_shake>` does.

The *kbond* keyword sets the dynamics force constant (used by the
*linearangle restrain* substitute); the minimization substitute is
automatically a fixed factor of 2000 stiffer, because minimization has
no timestep stability limit while a minimization-stiff restraint would
be unstable in dynamics.  Both defaults are expressed as multiples of
the Boltzmann constant so they scale with the unit system rather than
being tied to one unit choice: if *kbond* is not given it defaults to
:math:`5\times10^5 k_B` for dynamics (about 1000 in real units), so the
minimization default works out to :math:`10^9 k_B`, the value used by
:doc:`fix shake <fix_shake>`.  A *kbond* value set explicitly is
likewise the dynamics constant and is scaled by the same factor for
minimization.  Minimize first, then run, so the holonomic solver
tightens the bonds exactly once dynamics begins.

Restrictions
""""""""""""

This fix is part of the RIGID package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

The molecular topology (bonds) must be defined; an :doc:`atom_style
<atom_style>` such as *full*, *molecular*, or *bond* (or a hybrid atom
style including them) is required, and an :doc:`atom map <atom_modify>`
must be enabled.

Only one ``fix ilves`` instance may be defined at a time.  ``fix ilves``
and :doc:`fix shake <fix_shake>` may be used at the same time, but the
sets of constrained atoms *must not* overlap.

``fix ilves`` supports :doc:`run_style respa <run_style>`.  As for
:doc:`fix shake <fix_shake>`, the constraints are enforced at every
r-RESPA level using the level-dependent effective timestep.

All atoms of a constraint cluster must lie within the communication
cutoff of each other on every rank.  For small clusters (water, methyl,
hydrogen-only constraints) this is satisfied automatically; for clusters
that span large distances increase the cutoff with :doc:`comm_modify
cutoff <comm_modify>`.  The fix stops with an error if a constraint
partner is not available locally.

``fix ilves`` does not support dynamic topologies.  Fixes or commands
that add, remove, or change constrained bonds during a run (for example
:doc:`fix bond/create <fix_bond_create>`, :doc:`fix bond/break
<fix_bond_break>`, or :doc:`fix bond/react <fix_bond_react>`) must not
be applied to the constrained atoms.

Related commands
""""""""""""""""

:doc:`fix shake <fix_shake>`, :doc:`fix rattle <fix_shake>`,
:doc:`fix restrain <fix_restrain>`, :doc:`fix rigid <fix_rigid>`

Default
"""""""

The keyword defaults are *mode* = *converge*, *linearangle* = *error*
175, and *store* = *no*.  If *kbond* is not set, the dynamics restrain
substitute uses :math:`5\times10^5 k_B` (about 1000 in real units) and
the minimization substitute is a factor of 2000 stiffer (:math:`10^9
k_B`, as :doc:`fix shake <fix_shake>`); a *kbond* value sets the
dynamics constant and is scaled by the same factor for minimization.

----------

.. _Lopez-Villellas2025:

**(Lopez-Villellas)** L. Lopez-Villellas, C. C. K. Mikkelsen,
J. J. Galano-Frutos, S. Marco-Sola, J. Alastruey-Benede, P. Ibanez,
P. Echenique, M. Moreto, M. C. De Rosa, and P. Garcia-Risueno, "ILVES:
Accurate and Efficient Bond Length and Angle Constraints in Molecular
Dynamics", J. Chem. Theory Comput. 21, 8711-8719 (2025),
https://doi.org/10.1021/acs.jctc.5c01376

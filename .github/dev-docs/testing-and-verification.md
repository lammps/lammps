# Testing and Verification Methodology

How to validate LAMMPS changes beyond the standard unit/regression suites: choosing
the right gate for refactors vs re-implementations, force verification, benchmark
construction, and memory-checking conventions.

## Choose the right correctness gate

**Behavior-preserving refactor -> bit-identical gate.** When restructuring code that
must not change behavior (style consolidation, template refactors), gate on
byte-identical thermo output: small decks per style x {serial, 4 ranks} x {newton
on/off} x accelerator suffixes, thermo every step for ~50 steps; extract the thermo
table (`awk '/^   Step/{f=1} /^Loop time/{f=0} f'`) and `cmp` against a baseline built
from the pre-change commit.  Gotchas:

- Keep floating-point summation ORDER identical: when caching in-place accumulation in
  a local, seed the local with the current array value (`fxtmp = f[i][0]; ...;
  f[i][0] = fxtmp;`) -- zero-seeding plus `+=` at the end regroups the summation and
  changes ulps whenever the array is pre-accumulated (hybrid sub-styles).
  Template-int parameters (`if (FLAG)` on `template <int FLAG>`) fold at compile time
  and stay bit-identical.
- WARNING lines can land inside the thermo table and embed FLERR line numbers that
  shift with the source; strip them before comparing.
- `-sf X` silently falls back to the base style when no `/X` variant exists; once a
  port adds the variant, those runs legitimately change -- exclude them from the
  byte-compare and validate numerically instead.
- GPU mixed precision: expect ~1e-5 relative deviation vs CPU; compare against an
  established style's gpu-vs-cpu deviation as the reference scale.
- A test program that segfaults AFTER `[ PASSED ]` is destroying globals in undefined
  order at exit (for Kokkos: `lammps_kokkos_finalize()` never called).

**Re-implementation with an analytic spec -> gate on the spec.** When new code has an
independent closed-form oracle (analytic model results), the specification is the
pass/fail gate; a bit-comparison against the legacy implementation is only a
diagnostic to localize behavior changes.  A new-vs-legacy discrepancy where the new
code matches the closed form is a legacy bug to report, not a regression -- do not
tune to hide it.  Add a permanent regression test for every bug found.

## Force verification with fix numdiff

To check that forces equal -grad(E), use `fix ID all numdiff Nevery delta`
(EXTRA-FIX): it finite-differences the per-atom PE into `f_ID[1..3]` while saving and
restoring the analytic forces.  Needs `atom_modify map array` (granular/atomic styles
have no map by default).  Idiom: `variable ferr atom sqrt((f_nd[1]-fx)^2+...)` +
`compute maxerr all reduce max v_ferr`, `run 0`.  A conservative force gives err
~1e-7 (delta^2); a non-conservative one is orders of magnitude larger.  Always include
the stock style as a positive control in the same geometry and an equilibrium
configuration as a negative control.

**Equilibrium configurations hide cutoff-coordination bugs.** In cutoff-coordination
many-body potentials (AIREBO/REBO family), the coordination force is proportional to
the cutoff-function derivative, which is zero unless a neighbor distance falls inside
the transition region (rcmin..rcmax).  At equilibrium, bonds sit inside rcmin and
non-bonds beyond rcmax, so such bugs are invisible.  Expose them by straining the
system (`change_box all z scale 1.2 remap`) or using fractional-coordination clusters.

## Benchmark construction

Grow a benchmark from a small example by inserting a `replicate Nx Ny Nz` command
rather than hand-building geometry.  Watch for vacuum regions: replication tiles the
vacuum too, giving some MPI ranks little work.  Use the `processors` keyword to align
the decomposition with the material and/or a `balance` command -- or better, prefer a
solid periodic block with no vacuum for clean performance comparisons.

## Valgrind suppressions

Suppressions live in `tools/valgrind/*.supp` (globbed and concatenated into the build
directory at CMake configure time).  Author them GENERICALLY: anchor each on one
stable high-level frame (e.g. `fun:PMPI_Init`) and wildcard the rest with `...` --
narrow stacks break on the next MPICH/libfabric/Python update.  A suppression block
needs at least one concrete `fun:`/`obj:` frame (an all-`...` stack is a fatal syntax
error).  Validate against an existing binary with `valgrind --suppressions=...` and
confirm the finding moves into the "suppressed" count; the small masking risk is an
accepted trade-off.

# OpenMP (OPENMP package) Porting and Audit Guide

How to add `/omp` accelerator variants to LAMMPS styles and how to keep the
existing ones honest.  The OPENMP package (~330 files) threads force styles,
neighbor builds, and a growing set of computes; a full-package audit (PR #5163,
2026-08) produced the bug-class catalog at the end of this guide.  There are TWO
distinct porting models: force styles use the ThrOMP machinery, computes
deliberately do not.

## Architecture: FixOMP, ThrData, ThrOMP

- `-sf omp` (or `package omp N`) instantiates the internal fix `OMP` (`FixOMP`,
  `src/OPENMP/fix_omp.cpp`).  It owns one `ThrData` object per thread (per-thread
  force/torque/energy/virial accumulators) and calls `omp_set_num_threads()`, so
  the OpenMP runtime thread count always equals `comm->nthreads` -- bare
  `#pragma omp parallel` regions anywhere in LAMMPS get the right thread count
  without any per-style setup.
- Force styles derive DOUBLY: `class PairXOMP : public PairX, public ThrOMP`
  (`ThrOMP(lmp, THR_PAIR)`; bond/angle/dihedral/improper/kspace use their `THR_*`
  constant).  The constructor sets `suffix_flag |= Suffix::OMP;`.
- Inside `compute()`, the canonical shape (see `pair_lj_cut_omp.cpp`) is:
  a parallel region wrapped in `LMP_DEFAULT_NONE LMP_SHARED(...)` (compat macros
  from `omp_compat.h` -- always use them, plain clauses break some compilers),
  then per thread: `loop_setup_thr(ifrom, ito, tid, inum, nthreads)`,
  `ThrData *thr = fix->get_thr(tid)`, `ev_setup_thr(...)`, a templated
  `eval<EVFLAG,EFLAG,NEWTON_PAIR>(ifrom, ito, thr)` writing forces via
  `thr->get_f()` and tallying via the `*_tally_thr()` helpers, and finally
  `reduce_thr(this, eflag, vflag, thr)` to fold the per-thread arrays back.
- The per-thread force reduction changes summation order: `-sf omp` TRAJECTORIES
  are reproducible only at a FIXED thread count, never across thread counts.
- Set `respa_enable = 0` in the ctor unless the respa path is actually ported.

## Registration mechanics (all automatic -- do NOT edit build files)

- CMake `RegisterStylesExt` (`StyleHeaderUtils.cmake`) scans `*_omp.h` style
  headers; `src/OPENMP/Install.sh` step 1 auto-installs `*_omp.{h,cpp}` and
  skips children of uninstalled parent packages; `src/.gitignore` already has
  `/*_omp.h` `/*_omp.cpp` wildcards.  Just drop the files into `src/OPENMP/`
  and RECONFIGURE (`cmake -S cmake -B <build>`) so the glob re-registers.
- Docs for a new `/omp` variant: add the `o` letter in the matching
  `Commands_*.rst` cell (no commas, e.g. `(ko)`), an `.. index:: <name>/omp`
  entry, and append `, *<name>/omp*` to the style page's `Accelerator Variants:`
  line.  Do NOT add `.. versionadded::` for an accelerator variant.

## Porting a compute style: NO ThrOMP

Computes write per-atom OUTPUT, never `atom->f`, so the ThrOMP force-reduction
machinery is irrelevant.  Rules established with the 2026 compute-porting batch
(a dozen `/omp` computes in-tree to copy from, e.g. `compute_coord_atom_omp`):

- Subclass the base compute in `src/OPENMP/compute_<name>_omp.{h,cpp}`; thread
  with bare `#pragma omp parallel for` -- no ThrOMP, no ThrData, no `fix omp`
  interaction.  `Compute` has NO `suffix_flag` member (unlike `Pair`) -- do not
  set one.  Suffix resolution happens in `Modify::add_compute`.
- Most base headers need `private:` -> `protected:` for the subclass to reach
  members and helpers.  File-static constants/enums in the base `.cpp` are
  simply re-declared in the `/omp` `.cpp`.
- LINKAGE GOTCHA: a base helper marked `inline` but DEFINED in the base `.cpp`
  has internal linkage and is not callable from the subclass (undefined symbol
  at load time).  Remove the `inline` (the compiler still inlines within the
  base TU).
- Per-thread scratch: prefer thread-local `new[]`/`delete[]` INSIDE the parallel
  region, sized to the global-max neighbor count and clamped `>= 1`.  If the
  base keeps member scratch, refactor it into function arguments so the worker
  is reentrant (precedent: `ComputeOrientOrderAtom::calc_boop`).  `memory->grow`
  calls inside the per-atom loop must become thread-local allocation.
- Two-phase computes (build per-atom rows, then read neighbors' rows -- cna,
  cnp, entropy): use TWO `#pragma omp for` loops; the implicit barrier between
  them provides the ordering.  Per-thread error counters use
  `reduction(+:nerror)`.  All MPI collectives stay OUTSIDE the parallel region.
- Tallying to BOTH atoms of a pair (incl. ghosts, e.g. contact/atom) needs
  atomic increments; it stays bit-identical when every tally is an integer-
  valued `+= 1.0`.  Same for histogram computes (rdf): per-thread partial
  histograms, reduced into the shared array before the `MPI_Allreduce`.
- Known-hard cases, deferred deliberately: ML descriptor computes (sna/snad/
  snav/pace/pod -- need one descriptor engine per thread, memory-heavy, no
  template exists since SNAP has no /omp pair style; validate with a numeric
  tolerance, not bit-identity), `adf` (many member scratch arrays grown inside
  the loop), `fabric` (multi-pass floating-point tensor accumulation -- only
  tolerance-level reproducible).  SKIP: cheap reductions (`temp*`, `pressure`,
  `ke`, `pe`), delegating computes (`pe/atom`, `stress/atom`), external-library
  computes (`voronoi/atom`, `ptm/atom`).

## Validation recipes

- Compute `/omp` variants must be BIT-IDENTICAL: compare per-atom/array output
  at `-pk omp 1` vs `4`, at np=1 and np=2, and against the serial base style.
  Confirm the `/omp` variant actually resolved with `info computes` (occasional
  computes do not appear in the neighbor-list summary).
- When a threaded PAIR style would sit in the deck, keep the pair style and
  integrator UNTHREADED in all runs (identical trajectories) and vary only the
  compute under test (base vs `/omp`, threads via `package omp N`); comparing
  full `-sf omp` trajectories across thread counts is invalid (reduction order).
- Force styles are covered by the force-style YAML tests: the drivers
  automatically run `/omp` sub-tests off the base reference when the package is
  compiled in -- no new YAML needed, but remove stale `omp` entries from
  `skip_tests` lines when a port makes them runnable.
- PERFORMANCE RULE: measure before keeping any threading "optimization".
  OpenMP dispatch overhead makes fine-grained parallelism a net LOSS (per-cluster
  threading for constraint solvers was implemented, benchmarked, and rejected
  twice).  Benchmark or drop.

## Maintenance: copy-adapt drift and the 2026-08 audit bug classes

`/omp` files are copy-adapted snapshots of their base styles and DRIFT: base
bugfixes, new keywords, and lifted restrictions do not propagate by themselves.
Sync check method: side-by-side semantic diff of the `/omp` file against the
CURRENT base, plus `git log --follow` on the base file since the `/omp` file's
last real sync.  The full-package audit (PR #5163) found these recurring bug
classes -- check for each of them when touching or reviewing an `/omp` file:

1. **Group-mask drift**: per-atom output rows for atoms outside the compute
   group must be zeroed/skipped exactly like the base does (stale values were
   reported for non-group atoms).
2. **Inherited base bugs**: the copy faithfully reproduces base bugs (a
   division by zero shipped in both copies of one compute).  A bug found in
   either file is rarely alone -- fix base, `/omp`, and all other accelerator
   siblings together.
3. **Argument-parsing drift**: parsing code duplicated into the `/omp`
   constructor drifts worst (equal-style variables mis-bound, keywords
   silently ignored).
4. **Physics/feature drift**: `/omp` integrator fixes missed base fixes for
   years (2d disc moment-of-inertia, point-dipole rotation silently skipped).
5. **Orphaned OpenMP directives**: a `#pragma omp for` outside any ACTIVE
   parallel region runs serially without warning (two solver loops did).  But
   note the valid inverse pattern: `omp master`/`omp barrier` inside helpers
   that are CALLED from a parallel region are dynamically enclosed and fine --
   a textual scan for orphans produces false positives on those.
6. **Missing terms/masking**: contributions present in the base (special-bond
   masking in a matrix build) dropped from the threaded loop.
7. **Stale per-thread scratch tallied**: when the code path that FILLS
   per-thread scratch is conditional but the tally that READS it is not, stale
   data from a previous iteration enters energy/virial (sub-particle force
   buffers in line/tri pair styles).
8. **Stale restrictions**: `/omp` variants still refusing what the base has
   long supported (triclinic in one kspace style, ten years after the base
   lifted the restriction).
9. **`omp master` has NO implied barrier**: master-only sections that write
   shared state (reverse_comm results, size counters) race with worker threads
   unless followed by an explicit `#pragma omp barrier`; and thread-chunked
   allocation loops must span the SAME index range as the base's serial loop --
   chunking `[lfrom,lto)` over owned atoms only silently skips ghosts
   (the FixNeighHistoryOMP newton-on segfault class).

Related policy: a pair style that creates an internal fix (e.g. neighbor
history) should request it with `trysuffix=1` only if its own `/omp` children
actually work with the threaded fix variant; serial-only pair styles pass
`trysuffix=0` so `-sf omp` does not attach a threaded internal fix behind a
serial pair style.

Runtime race detection (beyond code review): a ThreadSanitizer build of the
OPENMP unit tests needs the OpenMP-aware tooling (libarcher / `OMP_TOOL`)
to avoid drowning in false positives from the OpenMP runtime itself.

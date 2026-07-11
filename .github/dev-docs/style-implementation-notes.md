# Style Implementation Notes (MPI, restart, buffers, parsing)

Hard-won implementation lessons for writing and modifying LAMMPS styles.  Read this
before writing a new pair/bond/fix/compute style or debugging parallel/restart
misbehavior in one.

## MPI parallelism patterns

**`newton on` vs `newton off` for bond/angle topology.** Under `newton on` (default),
each bond/angle is stored at only one of its participating atoms (typically the one
with the lowest local index).  That means a rank may see only *part* of a cluster
whose atoms span multiple ranks -- the other half of the cluster may live entirely in
topology arrays on other ranks.  Under `newton off`, every bond is stored at both
endpoints, so each rank sees all bonds touching its local atoms, but each bond must be
processed exactly once to avoid double counting.

**Per-atom accumulators in bond styles: reverse_comm is a newton-ON construct.**
For a bond (or angle/dihedral) style that builds per-atom accumulators by summing over
its topology: under `newton bond off`, a boundary-straddling bond appears in BOTH
ranks' bond lists and is computed twice, each rank handling its own local endpoint --
so each local atom's accumulator completes locally under an `if (i < nlocal)` guard
and needs only `forward_comm` to publish the finalized value to ghosts.  No
`reverse_comm`.  Folding ghost partial sums onto owners via `reverse_comm` is needed
only under `newton bond on` (canonical example: `BondBPMZero::calculate_manybody`
guards `comm->reverse_comm(this)` with `if (newton_bond)` but always calls
`comm->forward_comm(this)`).

**Determinism across rank counts.** LAMMPS trajectories are deterministic across MPI
rank counts (modulo ~1e-12 summation-order roundoff).  A large, immediate 1-vs-N
divergence almost always means decomposition-dependent random numbers in the setup
(`velocity create` without `loop geom`, `create_atoms ... random`, stochastic fixes),
NOT a force/communication bug.  Reproduce from an identical start (data/restart file,
or `velocity create ... loop geom`) before suspecting the style's forward/reverse comm.

**Warnings print on all ranks.** `error->warning()` writes to `lmp->screen`, which
defaults to stdout on EVERY rank; only the `log.lammps` file is rank-0-only.  Do not
"fix" per-rank warnings on the assumption they are suppressed -- guard with
`comm->me == 0` only when a single message is actually wanted.

## Restart and atom-migration mechanics in fixes

**`restart_peratom` and the `Atom::RESTART` callback are independent; gate BOTH.**
`atom->add_callback(Atom::RESTART)` registers the fix in `atom->extra_restart[]`,
whose entries are packed UNCONDITIONALLY on `write_restart` -- the callback, not the
flag, drives writing.  `restart_peratom` only feeds the read-back matching list used
to re-attach saved data.  A fix that must write no per-atom restart data has to leave
`restart_peratom = 0` AND skip the callback (plus the matching `delete_callback` in
the destructor); registering the callback with `restart_peratom = 0` writes data the
reader cannot re-assign -> silently corrupt restarts (read fails with err0016).

**A fix packing per-atom exchange data must set `maxexchange`.** The migration buffer
slack is `sum(fix->maxexchange) + BUFEXTRA` (1024 doubles); anything packed in
`pack_exchange()` beyond the declared `maxexchange` rides on that fixed cushion and
overflows the heap once a fix packs more than ~1024 doubles/atom.  Fixed upper bound
-> set `maxexchange` in the constructor; runtime-growing size -> set
`maxexchange_dynamic = 1` and update `maxexchange` as it grows (canonical example:
`fix_neigh_history` with `maxpartner`).

## Argument parsing and validity guards

**`ArgInfo` parses `[i][j]` only for prefixed references** (`c_ f_ v_ d_ i_ c2_ d2_`).
A bare keyword like `history[1][2]` falls through with `dim = 0` and unparsed indices.
For bare-keyword brackets: `ValueTokenizer(arg, "[]")` (consecutive separators
collapse), then convert with `utils::inumeric(FLERR, tok, false, lmp)` -- NOT
`next_int()`, which throws a non-LAMMPS exception (uncaught abort) on overflow.  Note
`utils::strmatch` supports `^ $ . * + ? [] ranges \s \S \w \W \d \D \i \f` and literal
`\[ \]`, but NOT groups `(...)`, branches `a|b`, or inverted classes `[^...]`.

**Periodic-output computes must self-guard timestep validity (err0007).** A `Fix`
advertises `peratom_freq`/`global_freq` and consumers validate at init; a `Compute`
has no freq member, so a compute whose output is only valid on certain steps must
check `update->ntimestep` inside its compute method and
`error->all(..., utils::errorurl(7))` (precedent: `compute_chunk_atom.cpp`).  The only
no-error alternative is the always-available pattern (recompute on demand).

**Box/ntypes-dependent guards belong in `settings()`, not the constructor.** The pair
constructor runs on BOTH the `pair_style` command path and the `read_restart` path;
only the former calls `settings()` afterward.  A `domain->box_exist` check (or any
`atom->ntypes`-derived sizing) in the constructor spuriously trips or under-allocates
during restarts and under `hybrid/overlay` (where the constructor commonly runs before
`create_box`, so `ntypes == 0` -> under-sized `map`, zero comm sizes, heap corruption).
Convention: fixed `comm_forward`/`comm_reverse` constants in the constructor; ntypes-
dependent checks and sizing in `settings()`/`coeff()`.

## Small conventions

- `src/pointers.h` includes and re-exports `<mpi.h>`, `<cstdio>`, `<string>`, and
  `<vector>`; virtually every class derives from `Pointers`, so do not add redundant
  includes for these -- only for headers it does not export (`<cmath>`, `<map>`,
  `<algorithm>`, ...).
- The `Memory` class does NO usage accounting; per-style `memory_usage()` is a rough
  estimate covering large (per-atom) allocations only.  Using `new[]` instead of
  `memory->create()` for small arrays does not "lose accounting" -- there is none.

## Recorded performance decisions

- **`fix ilves/omp` (per-cluster OpenMP threading) was implemented, benchmarked, and
  rejected twice** (2026): the constraint solver is only ~6% of step time even in the
  all-bonds best case, per-cluster work is too fine-grained for OpenMP dispatch
  overhead, and the net effect was a slowdown.  Do not re-propose it (same reasoning
  applies to SHAKE/RATTLE-style per-cluster threading).  Promising directions instead:
  single-rank cluster ownership in MPI mode; intra-cluster threading of the band
  Cholesky for the one-giant-cluster case -- benchmark first.

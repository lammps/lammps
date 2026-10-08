# oxDNA3 KOKKOS GPU performance commits (branch `oxdna3KK-kk-perf`)

Guide for benchmarking the commits on `oxdna3KK-kk-perf`, which sits on top of
`oxdna3KK-kk-fixes` at 2afff57fb5 (explicit fp32/fp64 conversions in the CG-DNA
styles).  How to build and run both codes and work through the options is in
`oxdna3-kokkos-gpu-handover.md`.  The item ids (P0, A1-A5, B1-B4, C1-C4, D1-D6, E1-E5,
F1) refer to the master list in section 7 of
`oxdna3-kokkos-gpu-performance.md` on branch
`claude/kokkos-oxdna3-perf-analysis-osblax`; ids S1-S13 refer to the second
search, listed at the end of this file.  Target: the oligomer benchmark
(dilute octamer duplexes), where KOKKOS CUDA is about 2.7x slower than the
standalone oxDNA CUDA code, while polybrick is within about 5%.

## How the commits were verified (no GPU was available)

- Kokkos Serial, OpenMPI with 4 ranks, `-D KOKKOS_DEBUG_RNG=on` (same Langevin
  random numbers as the CPU styles), 200 steps of `fix nve/asphere` +
  `fix langevin ... angmom`, compared with the CPU styles at every thermo step
  (total and per-sub-style energies via `compute pair`, bond energy).
- Decks: oxDNA1/oxDNA2/oxRNA2/oxDNA3 lj-unit duplexes, oxDNA3 unique base
  pairs, an oligomer proxy (duplex2 replicated 4x4x4), a dense proxy, a
  minimize + run deck, a `rerun` deck, trimmed lists (the default) with a 0.3
  skin, and oxdna3/xstk listed before oxdna3/hbond (with and without per-atom
  energies), and one deck with `pair_modify neigh/trim no`.
- Each deck ran with `neigh half newton on` and `neigh full newton off`, in a
  normal build (host kernels) and in a test build that forces the GPU
  (screened-pair) code paths of hbond, xstk and coaxstk onto the host.
- Exact commits were bit-identical to their parent on all decks;
  rounding-level ones agree with the CPU styles to 1e-12 or better after 200
  steps (the minimize deck to about 1e-9, as before the changes).
- Every commit compiles with nvcc 12.6 for sm_89 (`Kokkos_ARCH_ADA89`,
  `KOKKOS_PREC=mixed`).  Registers, stack and kernel parameter size from
  `cuobjdump -res-usage` are listed below.  At HEAD, the affected files also
  compile for sm_89 with each non-default switch below
  (`OXDNA_KK_FUSE_HBXSTK=0`, `OXDNA_KK_TWO_PHASE=1`,
  `OXDNA_KK_SCREENED_PER_ATOM=1`, nonzero launch bounds for all three kernel
  classes); combining the two exclusive options stops at the intended
  `#error`.
- ctest (Kokkos Serial + MPI build with CG-DNA, MOLECULE, ASPHERE; 1060
  tests) passes, including `AtomStyles`, `AtomStylesKokkos` and the
  `FixTimestep` tests.  Five tests fail only because the container runs as
  root: three MPI tests pass with Open MPI's run-as-root settings, and
  `Platform` and `TextFileReader` expect an unreadable file, which root can
  read.  There are no force-style tests for the oxDNA pair or bond styles.
- Float conversions: every commit was compiled with clang 18 in KOKKOS mixed
  precision with `-Wimplicit-float-conversion -Wdouble-promotion` (the check
  used for kk-fixes 2afff57fb5); the oxDNA sources and the other files changed
  on this branch have no such warnings in any commit, and none at HEAD in
  single precision either (apart from Kokkos-internal and existing
  `verlet_kokkos.cpp` lines).

## Commits (oldest first)

"Exact" means bit-identical results; "rounding" means a different summation
order (on GPUs the atomics are nondeterministic anyway).

| # | Hash | Id | Change | Result | What to measure |
|---|------|----|--------|--------|-----------------|
| 1 | 193529c43c | P0 | Monotonic `Neighbor::nbuild` counter replaces `ncalls` as the cache key of the 3'/5' tables (stale tables in `rerun`, which resets `ncalls`) | fixes rerun | correctness only |
| 2 | c3b61e9010 | B1 | Leave zero-weight special pairs (1-2 bonded) out of the screened pair list | exact | fewer threads in hbond/xstk/coaxstk |
| 3 | eb3a9acc2d | B3a | coaxstk: radial factor first, drop unused cosphi3 terms | exact | coaxstk kernel |
| 4 | f13e648f7f | B3b | hbond: return early for base pairs with `epsilon_hb == 0` | exact | hbond kernel |
| 5 | 8c3c97df40 | B2 | coaxstk GPU kernel runs over a separate list of strand-end pairs (built by the npair fix) | exact | coaxstk kernel (only strand-end pairs get a thread) |
| 6 | 24439021f3 | E5 | stk styles no longer request an unused neighbor list | exact | Neigh time |
| 7 | 95333e88c6 | B4 | stk: skip tetramer type loads when the tables are uniform; excv: hoist per-atom loads | exact | stk, excv kernels |
| 8 | 417388398b | C4 | excv looks up 3'/5' neighbors directly; removes the per-slot table (anum x maxneigh x 4 ints) | exact | excv kernel, rebuild time, memory |
| 9 | e82f01d4db | A5 | fene: `Kokkos::atomic_add` on the view instead of atomic-trait local views; STACK 64 -> 0 | exact | fene kernel |
| 10 | a61f378bc4 | A5/A4 | stk: accumulate, then one round of 12 atomics; STACK 96 -> 32 | rounding | stk kernel |
| 11 | 19c48f9fd0 | D4 | Tunable launch bounds: `-DOXDNA_KK_{ATOM,PAIR,BOND}_{MAXT,MINB}` | exact | sweep MAXT/MINB |
| 12 | a39288a5de | A4 | Screened kernels: r x f torque and pure torque summed, one torque atomic per atom | rounding | hbond, xstk, coaxstk |
| 13 | 12d142b791 | A4b | Optional `-DOXDNA_KK_SCREENED_PER_ATOM=1`: one thread per atom a over its screened pairs, a-side in registers | exact (default) / rounding | A/B the macro |
| 14 | 2862235dee | C1a | Frame vectors in LayoutRight | exact | all consumers |
| 15 | f217d44399 | C1b | LRF kernel writes a packed per-atom record (x, nx, ny, nz; 16 floats) read by all force kernels | exact | all force kernels |
| 16 | 4bdd797ab8 | D2a | LRF frame kernel also zeroes f and torque (Verlet skips its two zero kernels when safe) | exact | 2 fewer launches per step |
| 17 | 58c8373c2e | C3 | hbond coefficients packed into one struct per type pair; PARAM 16304 -> 3704 B | exact | hbond launch latency |
| 18 | 87f3af21a7 | C3 | Same for coaxstk, stk, excv; PARAM 14-15 KB -> 3.4-3.9 KB | exact | launch latency |
| 19 | 0831c23424 | A3 | Optional `-DOXDNA_KK_TWO_PHASE=1`: radial prefilter kernel compacts the hbond/xstk pairs, then the full kernels run over the survivors | exact | A/B the macro |
| 20 | e3213e282e | new | Fix OXDNA/NPAIR/kk requests its list at the screen cutoff (trimmed from the dh-sized list) | exact (pair order may change with trim) | npair count/fill kernels |
| 21 | 9c31bbf78e | new | hbond/xstk/coaxstk request neighbor lists only for their host kernels; oxdna3/xstk never | exact | Neigh time (with trim: 3 lists instead of 6) |
| 22 | 46acb30b6d | A2 | Fused hbond + oxdna3/xstk kernel (default on, `-DOXDNA_KK_FUSE_HBXSTK=0` disables) | rounding | A/B the macro |
| 23 | 26e529222a | S11 | LRF no longer writes the unused separate nx/ny/nz arrays (36 B per atom and step) | exact | LRF kernel |
| 24 | 2989879798 | S11 | LRF record carries type (column 3) and qeff (column 13); excv and dh read the neighbor's type/qeff there, next to its position | exact | excv, dh kernels |
| 25 | c2ec198b3c | S9 | excv: skip neighbors beyond the center-of-mass cutoff before building site vectors; topology loads only for 1-2 special neighbors; uniform coefficients read from the functor | exact (see note) | excv kernel |
| 26 | 80291b3036 | S10 | dh coefficients packed into one struct per type pair (was 7 views) | exact | dh kernel |
| 27 | 95ceee49c7 | S13 | npair fix: one count/scan/readback/fill for the screened and coax lists (4 launches + 1 readback per rebuild instead of 6 + 2) | exact | rebuild steps, small skin |
| 28 | 94f8461506 | S12 | stk angles from the cross-product norm and atan2 instead of acos and sin | rounding | stk kernel (no stack left) |
| 29 | 79070ecf2d | S12 | fene force with one division instead of three | rounding | fene kernel |

Note on #25: the bonded base-base terms are now evaluated only for 1-2 special
neighbors.  atom style oxdna sets the 3'/5' neighbors from the Bonds section, so
they always are 1-2 neighbors; results could only differ if bonds were deleted
while stale 3'/5' neighbors were kept.

## Compile-time switches (all in `src/KOKKOS/mf_oxdna_kokkos.h`)

Pass them with `-D CMAKE_CXX_FLAGS="-DOXDNA_KK_..."` (with nvcc_wrapper).

| Macro | Default | Effect |
|-------|---------|--------|
| `OXDNA_KK_ATOM_MAXT/MINB` | 64/1 (CUDA), 128/1 (HIP) | launch bounds of per-atom kernels (excv, dh, stk, ...) |
| `OXDNA_KK_PAIR_MAXT/MINB` | 0/0 (none) | launch bounds of the screened-pair kernels (hbond, xstk, coaxstk, fused) |
| `OXDNA_KK_BOND_MAXT/MINB` | 0/0 (none) | launch bounds of the fene kernel |
| `OXDNA_KK_SCREENED_PER_ATOM` | 0 | 1 = one thread per atom over its screened pairs |
| `OXDNA_KK_TWO_PHASE` | 0 | 1 = radial prefilter + compacted hbond/xstk kernels |
| `OXDNA_KK_FUSE_HBXSTK` | 1 | 0 = separate hbond and xstk kernels |

`OXDNA_KK_TWO_PHASE` and `OXDNA_KK_SCREENED_PER_ATOM` cannot be combined;
either one disables the fused kernel.

## Kernel resources (sm_89, mixed precision)

| Kernel | kk-fixes aed39b4499 REG/STACK/PARAM | HEAD REG/STACK/PARAM |
|--------|--------------------------|----------------------|
| oxdna3/xstk | 95 / 0 / 3904 | 98 / 0 / 4056 |
| hbond | 78 / 0 / 16216 | 83 / 0 / 3768 |
| coaxstk | 69 / 0 / 14584 | 67 / 0 / 3512 |
| excv | 94 / 0 / 15304 | 109 / 0 / 4168 |
| dh | 72 / 0 / 4096 | 64 / 0 / 3096 |
| stk | 72 / 96 / 14816 | 117 / 0 / 3440 |
| fene | 42 / 64 / 2432 | 48 / 0 / 2496 |
| LRF | 34 / 0 / 1512 | 38 / 0 / 1336 |
| fused hbond+xstk | - | 123 / 0 / 13872 |
| two-phase radial (hbond / xstk) | - | 32 / 39 registers |

No kernel uses a stack frame any more (the last one, in stk, was the slow
path of `sin()`, removed by #28).  stk (117) and excv (109) now use more
registers than before #28 and #25; if their kernels come out slower, try the
`OXDNA_KK_BOND_*` / `OXDNA_KK_ATOM_*` launch bounds.  The fused kernel carries
copies of both styles, hence its larger parameter block.  HEAD numbers are
for the branch on kk-fixes 2afff57fb5 (nvcc 12.6, KOKKOS_PREC=mixed).

## Suggested benchmark protocol

1. Baseline: `oxdna3KK-kk-fixes`, then HEAD of `oxdna3KK-kk-perf`, oligomer
   and polybrick, with `nsys`/`ncu` or the Kokkos simple kernel timer.
2. If HEAD is faster, bisect by groups: after #10 (exact cleanups + stack
   removal), after #18 (layout, launch count, packed coefficients), after #22
   (two-phase, neighbor lists, fused kernel), then #23-29 (second search).
3. A/B the macros at HEAD: `OXDNA_KK_FUSE_HBXSTK=0`,
   `OXDNA_KK_TWO_PHASE=1` (implies no fusion), `OXDNA_KK_SCREENED_PER_ATOM=1`,
   and a sweep of `OXDNA_KK_PAIR_MAXT` in {64, 128, 256} with `MINB` in {1, 2, 4}
   (the fused kernel uses 123 registers).
4. Neighbor lists: trimming is now the default (see below); time the Neigh
   section and the npair fix kernels, and try a smaller skin.

## Neighbor list trimming (kk-fixes commits e1a8c85f05, e2f233566c)

- The oxDNA pair cutoffs now include the distances of the interaction sites
  from the centers of mass, and the CG-DNA styles no longer turn trimming off,
  so the lists of the hybrid/overlay sub-styles are trimmed to their own
  cutoffs by default (`pair_modify neigh/trim no` restores the old behavior).
  kk-fixes reports a 1.2x to 1.7x speedup of oxDNA runs from this alone.
- The larger cutoffs grow the master list (e.g. from 5.64 to 6.60 sigma with a
  skin of 2.0 and oxdna3/dh).  With trimming, excv loops over its own, much
  shorter list instead of the dh-sized one.
- Trimming does not apply to the list of fix OXDNA/NPAIR/kk, which is not a
  pair sub-style: without #20 it is still a copy of the dh-sized list, and its
  count and fill passes (and those of the coaxstk list) loop over all of those
  neighbors.  With trimming on, the hbond, xstk and coaxstk lists would each
  cost a trim pass per rebuild on GPUs although only their host kernels read
  them; #21 removes them.  With the oxDNA3 styles a GPU run now builds 3 lists:
  dh (binned), the npair fix list (binned) and excv (trimmed from the fix
  list).
- Baselines: compare against kk-fixes at 2afff57fb5 (the base of this
  branch), so that the trimming speedup and the conversion fixes are not
  attributed to this branch.

## Not implemented

| Id | Item | Reason |
|----|------|--------|
| A1 | Per-atom fused kernel of all oxDNA3 terms (full list, no atomics) | New kernel architecture; much larger than the rest together.  Design: one thread per atom over a full list, bonded n3/n5 terms from the atom's own side, all nonbonded terms accumulated in registers, one write, energies by reduction; enabled only for the full oxDNA3 set under hybrid/overlay with a full list.  Worth it only if #22 plus the macros leave a large gap. |
| - | One kernel per edge as in the standalone code | not practical in LAMMPS |
| C2 | Cache the interaction site vectors in the LRF record | the sites are 1-3 FMAs from the frame vectors already loaded by #15; storing them adds loads instead |
| D1 | Slim functors for all hot kernels | the parameter size goal is mostly met by #17-18 (14-16 KB -> 3.4-4 KB) |
| D3 | Cheaper `check_distance` | the host synchronization is inherent (the standalone code also synchronizes to decide on rebuilds) |
| E1 | Build the bond 3'/5' table once for fene and stk | rebuild-only work, negligible per step |
| E2 | Skip the host copy of the bond list | rebuild-only; risky, since host readers of the bond list are hard to enumerate |
| E3 | Avoid reading back the screened pair count | rebuild-only, needed to size the list |
| E4, D5 | Input tuning (skin, bin size, trim) and benchmark fairness | input/benchmark notes, see above |
| D6 | Scatter-view duplication | OpenMP backend only; on GPUs the scatter views are atomic without duplication |
| F1 | Compile time of unused template instantiations | no runtime effect |
| - | Double precision constants in the kernels | handled separately |

## Second search: further candidates (not implemented yet)

From an audit of the per-step path outside the force kernels and of the
remaining kernels at #22.  Items S9-S13 are #23-29 above.

Input settings (no code change; check the benchmark inputs):
- `fix nve/dotc/langevin` has no KOKKOS version: two full host round trips of
  all atom data per step.  Use `fix nve/asphere` + `fix langevin ... angmom`.
- Atom sorting falls back to the host for this atom style (every 1000 steps by
  default): use `atom_modify sort 0 0.0`.
- `fix balance` on one GPU forces full host round trips on each rebuild;
  `fix print` disables the fused integrator every step; `neigh_modify every 1
  check yes` adds a host wait every step; `bond_style hybrid` with one
  sub-style adds a device-to-host copy every step (see S1).

Code:

| Id | Item | Expected effect | Risk |
|----|------|-----------------|------|
| S1 | `BondHybridKokkos::compute()` copies the bond counts to the host every step (`bond_hybrid_kokkos.cpp:140-141`); they are current from the rebuild | 1 blocking D2H per step | low |
| S2 | `Kokkos::fence()` before every reverse comm (`verlet_kokkos.cpp:609`), not needed on one rank | 1 fence per step | low |
| S3 | any END_OF_STEP fix disables the fused integrator on every step (`verlet_kokkos.cpp:838`); test only the steps where it fires | 1 launch per step | low |
| S4 | merge the two `fix langevin/kk` kernels (force and angmom), one RNG state per atom | 1 launch per step | medium |
| S5 | fold the quaternion forward-comm kernel into the position kernel | 1 launch per step | medium |
| S6 | fuse stk and fene over the bond list (same 3'/5' table; fall back with bond hybrid) | 1 launch, about half the bond atomics | medium |
| S7 | fuse excv into the dh neighbor loop (same backbone site distance, dh list covers excv) | 1 launch and one list pass | medium |
| S8 | fold coaxstk into the fused hbond+xstk launch | 1 launch | medium |
| - | intermediate force/torque math declared `KK_ACC_FLOAT` runs in double with `KOKKOS_PREC=mixed` (1/64 rate on consumer GPUs) | check the FP64 pipe with ncu first | precision decision |
| - | CUDA graph of the fixed per-step kernel sequence; fast-math flags for the oxDNA sources | only if nsys shows launch gaps / after measuring | high |
| - | oxrna2/stk has the same acos/sin pattern as stk before #28 | oxRNA2 only | low |

## H100 results and the follow-up commits (branch `oxdna3KK-kk-perf-h100`)

Measured on one node with 4x H100 (CUDA 12.8, `Kokkos_ARCH_HOPPER90`,
`KOKKOS_PREC=mixed`), run with `-k on g 1 -sf kk -pk kokkos neigh half`,
benchmark inputs unchanged.  Each commit was profiled before and after with
nsys (per-step kernel table, GPU idle time, memcpys, host syncs) and ncu
(occupancy, active threads per warp, sectors per request), checked against
the CPU build (step-0 energy of every sub-style for oligomer 512 and
polybrick 136K; duplex2 NVE trajectory; rerun of a duplex2 dump), and kept
only if no benchmark case got slower.

### Time per step (ms, 1 GPU)

| System | kk-fixes | kk-perf | this branch | standalone oxDNA | branch / standalone |
|--------|---------:|--------:|------------:|-----------------:|--------------------:|
| oligomer 128 nt   | 0.122 | 0.110 | 0.072 | 0.065 | 1.11 |
| oligomer 1k nt    | 0.128 | 0.122 | 0.077 | 0.066 | 1.17 |
| oligomer 8k nt    | 0.141 | 0.132 | 0.085 | 0.079 | 1.08 |
| oligomer 65k nt   | 0.279 | 0.255 | 0.191 | 0.163 | 1.17 |
| oligomer 524k nt  | 1.443 | 1.375 | 0.97  | 1.003 | 0.97 |
| oligomer 4.2M nt  | 11.44 | 11.44 | 6.71  | 10.10 | 0.66 |
| polybrick 136K nt | 0.722 | 0.715 | 0.485 | 0.787 | 0.62 |
| polybrick 272K nt | 1.168 | 1.163 | 0.777 | 1.621 | 0.48 |
| polybrick 543K nt | 2.103 | 2.107 | 1.377 | 3.234 | 0.43 |
| polybrick 1087K nt| 4.093 | 4.333 | 2.570 | 7.175 | 0.36 |
| polybrick 2174K nt| 8.17  | 8.450 | 5.000 | 14.53 | 0.34 |

Polybrick runs are sensitive to other jobs on the node (host work and small
copies in the reneighboring steps): the numbers above are from one job at a
time.  The oligomer input applies the Langevin thermostat every step, the
standalone input its thermostat every 103 steps (about 5-9 percent of the
LAMMPS time).

### Commits on top of `oxdna3KK-kk-perf` (oldest first)

| Commit | Change | Effect on H100 |
|--------|--------|----------------|
| cb44d93ba1, f7009fa385 | from `lewis/oxdna3KK`: oxRNA2 scatter-view simplification, review cleanups | neutral |
| 6879c18d96 | `Neighbor::morph_copy_trim()`: trim a per-type (`cut_fixed`) list from a list with a longer default cutoff, so pair excv is trimmed from the dh list instead of binned separately | polybrick ~6% |
| 50cbaa778e | VerletKokkos: drop the global fence before `reverse_comm()` | ~15 us/step |
| 05e8f60075 | excv and dh: share the neighbors of an atom between up to 4 threads when there are too few atoms to fill the GPU; excv sums the b-side force over the site pairs of a neighbor before one atomic update | 8k nt: excv 18->8.5 us, dh 16->9.3 us |
| 0c56256472 | fix OXDNA/LRF/kk: write the packed record with 16-byte stores | LRF 38->14 us at 136K |
| f0a244ac61 | bond hybrid/kk: no per-step copy of the sub-style bond counts to the host | ~25 us/step polybrick |
| de0f6a7ea9 | NPairKokkos: flat build kernel when most bins are empty (the team kernel sizes every team by the fullest bin) | polybrick 4-8% |
| 829165c694, c30ca8e620, 00fe7a5442 | dh, excv, stk, fene: load the neighbor's packed record with 16-byte loads (excv: position/type chunk first, frames only inside the COM cutoff) | dh -26%, stk -33%, excv -10% |
| caf73136a5 | VerletKokkos: fuse the integration also with end-of-step fixes on the steps on which none of them is invoked (`fix print`) | small oligomers ~3% |
| 61d13d0227 | NeighborKokkos::check_distance(): kernel writes a stamped flag to pinned host memory instead of a reduction plus copy | 8k nt ~6% |
| 0918675bd3 | fix langevin/kk: angmom thermostat in the force kernel, sharing the RNG state | 1 launch; small oligomers ~8% |
| 8738e5955c | fused hbond+xstk: first kernel compacts the screened pairs that pass a radial test (~15%), second kernel evaluates only those, so warps no longer idle in the angular terms (used when the list fills the GPU more than twice) | 524k nt+ and polybrick 10-15% |
| a4d040459a | lighter radial pre-test (one base-site distance, per-type-pair xstk range over all 3'/5' contexts) | 1-2% |
| 7e1feaa7fb | comm pbc/swap2, atom2bin, xstk 3'/5' table: grow with headroom instead of reallocating on every rebuild | removes cudaMalloc/Free from rebuilds |
| 689679a61f | dh: full neighbor list on GPUs when there is less than one neighbor per atom at the mean density (dilute solutions); dense systems keep the half list | oligomer 65k+ 8-13% |
| 7f05d92c0b | excv: 3'/5' context types from a per-atom table built per reneighboring, not tag->index map lookups in the kernel | 4.2M nt 15%, 2.2M nt polybrick 14% (hash atom map) |

### Tried on H100 and not kept

- full neighbor list for dh in dense systems (polybrick, ~40 neighbors/atom): +20%;
- more threads per atom (nsplit) at large N: excv 154/172/216 us for 1/2/4 threads;
- sizing the second fused kernel from the radial count read back to the host: the sync costs more than it saves;
- two-pass excv: 82-84% of the list pairs pass the COM cutoff, nothing to compact;
- bin sizes other than the cutoff (0.5x, 1.5x, 2x): 1-20% slower;
- the `OXDNA_KK_TWO_PHASE` and `OXDNA_KK_SCREENED_PER_ATOM` switches: equal or slower than the compacted fused kernel; launch-bound variants (64/128/256 threads) neutral.

### Remaining cost

- 8k-65k nt oligomer: launch bound.  14 kernel launches per step at ~4 us of
  host time each and one host sync for the reneighbor decision; the GPU is
  idle ~25-30 us of an 85 us step.  Fewer launches need fusions of kernels
  owned by different styles (stk+fene, coaxstk into hbond+xstk, langevin into
  the integrator, the two forward-comm pack kernels), or a reneighbor check
  that does not block every step.
- large systems: the fused hbond+xstk kernel (~25%), excv and dh (~15% each);
  the integration, frame and thermostat kernels run at ~60% of the memory
  bandwidth, partly because the quaternions are stored in double.
- polybrick reneighboring steps (~10-13% of the time): blocking small copies in
  exchange/borders, the host copy of the bond list (needed by non-KOKKOS
  readers of `neighbor->bondlist`), and a thrust temporary allocation in the
  atom map sort.

### Lewis' regression report (polybrick 136K, kk-fixes vs old `oxdna3KK`)

The slowdown comes from e1a8c85f05 (site offsets in the pair cutoffs): the dh
master list cutoff grows from 1.78 to 2.74 and the neighbor build cost about
doubles.  The larger cutoff is needed: the old code dropped interactions (step
0 dh energy at 136K: CPU 0.0031256, old branch 0.0028437, fixed 0.0031256).
The screened list of the old code also used half the skin, which can miss
pairs between rebuilds.  H100 polybrick 136K: old 0.60, kk-fixes 0.71, this
branch 0.485 ms/step.

## Tooling used for verification

Not part of the branch; described so it can be recreated:

- builds with `BUILD_SHARED_LIBS=on`, Ninja, ccache, `PKG_KOKKOS=on`,
  `Kokkos_ENABLE_SERIAL=on`, `-D KOKKOS_DEBUG_RNG=on`, plus a CPU-only
  reference build with the same packages (CG-DNA, MOLECULE, ASPHERE);
- the GPU-path test build replaces `execution_space != HostKK` by `true` (and
  `== HostKK` by `false`) in `pair_oxdna_hbond_kokkos.cpp`,
  `pair_oxdna2_coaxstk_kokkos.cpp` and `pair_oxdna_xstk_kokkos.cpp`, and sets
  `force_screening_all_backends = true` in the constructor of
  `FixOxdnaNpairKokkos`;
- CUDA compile checks without a GPU: NVIDIA's redistributable nvcc 12.6,
  `nvcc_wrapper`, `Kokkos_ARCH_ADA89`, building only the `KOKKOS/*oxdna*`
  objects, then `cuobjdump -res-usage`.

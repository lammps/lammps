# Testing and Verification Methodology

How to validate LAMMPS changes beyond the standard unit/regression suites: choosing
the right gate for refactors vs re-implementations, force verification, test pitfalls,
platform-dependent (ARM64/macOS) failures, GPU package debugging, benchmark
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
  shift with the source; strip them before comparing.  Likewise strip timing columns
  (`spcpu`, `cpu`, ...) -- they never reproduce.
- Runs using `fix balance` are not reproducible even at a fixed rank count:
  `Irregular` migration receives with `MPI_ANY_SOURCE`, so atom ordering -- and with
  it summation order and roundoff -- differs from run to run and trajectories diverge
  chaotically.  Gate such changes on analytic equivalence or disable balancing for
  the comparison; never on re-run matching.
- `-sf X` silently falls back to the base style when no `/X` variant exists; once a
  port adds the variant, those runs legitimately change -- exclude them from the
  byte-compare and validate numerically instead.
- GPU mixed precision: expect ~1e-5 relative deviation vs CPU; compare against an
  established style's gpu-vs-cpu deviation as the reference scale.
- A test program that segfaults AFTER `[ PASSED ]` is destroying globals in undefined
  order at exit (for Kokkos: `lammps_kokkos_finalize()` never called).
- OPENMP pair trajectories are NOT thread-count invariant: the per-thread force-array
  reduction changes the summation order, so `-sf omp` runs only byte-compare at a
  FIXED thread count.  Never gate on matching across different `OMP_NUM_THREADS`.
- A failure seen only against force-style references generated with the kspace
  contribution suppressed (`pair_modify compute no` configs) can be a harness
  artifact -- that setup leaves the pair style's energy/virial accumulators in an
  undefined lifecycle.  Cross-check with a realistic full deck before root-causing
  (an ELECTRODE/omp "garbage results" report evaporated this way, 2026-08).

**Restart/settings I/O -> behavioral assertions, not just round-trips.**  A
write/read/write byte-compare only proves the new code reads back what it writes; a
setting missing from BOTH write paths passes it silently.  Add one assertion per
setting that the restored value changes observable behavior (example: a granular
`limit_damping` restart test must show the contact force clamped to zero after the
read-back, not just identical restart bytes).

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

## Test-system pitfalls

**Never use a huge simulation box in tests.**  Memory scales with box VOLUME, not
just atom count: neighbor-list bins and atom-sorting bins are allocated over the
whole box, so an almost empty box of 20000^3 length units OOMs or dies with "Too
many atom sorting bins".  For large coordinate *values* (e.g. a single-precision
resolution check) prefer a unit test that calls the class directly and runs no
simulation at all.  If a huge box is unavoidable, disable both volume-scaling
allocations:

```
neighbor        2.0 nsq        # N^2 neighbor list, no binning grid
atom_modify     sort 0 0       # no atom sorting bins
```

The same caution applies to anything else that multiplies out over the box volume:
fine `fix ave/grid`/`dump grid` grids and kspace meshes.

**Portability and silent-skip pitfalls in unit tests.**

- gtest's `ContainsRegex`/`MatchesRegex` are NOT portable: on Windows gtest falls back
  to its own regex engine with no character classes, groups, or alternation, while the
  POSIX path lacks `\d`.  For content checks use `utils::strmatch()` (via the
  `ASSERT_MATCH` idiom in `unittest/commands/test_info.cpp`); note the bundled
  tiny-regex also has no `{n}` repetition.
- CTest reports "Passed" even when EVERY gtest case in the binary skipped, and a
  force-style YAML with a missing prerequisite or `input_coeffs` entry self-skips
  silently (one shipped test never ran for weeks that way).  After adding a test,
  verify from the gtest/junit output that its cases actually EXECUTED, not just that
  ctest went green (see the note in `doc/src/Developer_unittest.rst`).
- End states of iterative algorithms (minimizers, converged solvers) are not portable
  reference data -- assert energies/box dimensions with tolerances, never per-atom
  end-state positions.
- A CI failure confined to one Linux job usually means `LAMMPS_SIZES=bigbig` (64-bit
  `tagint`/`bigint`); a minimal bigbig build (`-D LAMMPS_SIZES=bigbig` plus only the
  needed packages) reproduces it in a few minutes.
- Multi-rank gtests (`add_mpi_test`): `TEST_FAILURE` matches screen output that only
  rank 0 prints, so the other ranks return early and deadlock in the next collective
  call.  Check error paths by the exception text without an early return (the
  `expect_error()` helper in `unittest/commands/test_compute_voronoi.cpp`).
- An exception thrown between `BEGIN_HIDE_OUTPUT` and `END_HIDE_OUTPUT` leaves the
  gtest stdout capturer active; the NEXT capture aborts with "Only one stdout
  capturer" and the real error stays hidden.  Rerun the test with `-v`
  (`TEST_ARGS=-v` under ctest) to see it.
- An intermittently failing ("flaky") test usually means uninitialized memory in the
  style under test, not a harness or ctest-parallelism problem: run the driver under
  `valgrind --track-origins=yes` first.  On glibc, `MALLOC_PERTURB_=165` turns reads
  of stale or uninitialized heap memory into visibly wrong values without valgrind's
  cost (a cheap way to detect per-atom values left unset for atoms outside a group).
- Grep-filtering gtest output for "OK" hides crashes; check the exit status or the
  `[  PASSED  ]` summary line.
- A crash in untouched core code after a large merge, seen only in a long-lived
  incremental build directory, is often stale objects (ABI drift, e.g. renumbered
  bitmask constants): rebuild cleanly before debugging it.

**Test executables do not inherit the lammps target's PRIVATE compile definitions.**
Feature macros like `LAMMPS_ZLIB` (from `WITH_ZLIB`) are PRIVATE to the library, so
a unit test that must mirror the library configuration adds the define itself under
the same CMake option (see `test_vtk_writer` in `unittest/formats/CMakeLists.txt`).

## Platform-dependent test failures (ARM64, macOS)

The unit tests run in CI on Linux x86-64, Linux aarch64 (`unittest-arm64.yml`),
macOS, and Windows.  Linux aarch64 differs from x86-64 in two ways that matter
numerically; both can be emulated on an x86-64 machine in a separate build directory:

| Linux aarch64 property | x86-64 emulation |
|---|---|
| GCC contracts `a*b + c` into fused multiply-add (FMA) by default | `-D CMAKE_CXX_FLAGS="-mfma -ffp-contract=fast"` |
| plain `char` is UNSIGNED (signed on x86-64 and macOS) | `-D CMAKE_CXX_FLAGS=-funsigned-char` |

The macOS runner (`macos-26-intel`) is x86-64: a macOS-only failure comes from Apple's
libm or clang (e.g. a last-bit difference in `sin()`), never from FMA or char
signedness.  A clang build on Linux (`cmake/presets/clang.cmake`) covers the compiler
part.

Triage a failure seen on only one platform in this order, and do NOT tag the test
`unstable` before doing so -- the tag hides real portability bugs:

1. Random numbers seeded from coordinates or hashed bytes (`velocity ... loop geom`,
   `displace_atoms ... random`, `set ... random`, `fix atom/swap`, ...).
2. OpenMP data race: fails only in the `.omp` sub-test, and intermittently; weakly
   ordered ARM memory exposes races that x86-64 hides.
3. FMA: rebuild with the flags above and compare the driver's `-s` error statistics.
4. Near-zero relative comparisons in the reference data (see below).
5. Only then the `epsilon` tolerance.

Code bugs of this kind found in the 2026-09 survey of all `unstable`-tagged tests:

- **Arithmetic on plain `char`.**  `RanPark::reset(seed, coord)` hashed the bytes of the
  coordinates as plain `char`, so every coordinate-seeded random stream differed
  between x86-64 and aarch64 (36 of 88 failing tests).  Never compute with, hash, or
  compare plain `char` values; use `signed char`, `unsigned char`, or `uint8_t`.
- **Relying on exact cancellation.**  Under FMA, `r*r - m*m` with `r == m` is not 0.0
  but the rounding error of one product; a Voro++ radical-tessellation cutoff test
  silently returned wrong cells on all FMA platforms this way.  Write differences of
  squares as `(r - m) * (r + m)` or avoid depending on an exact zero.
- **Sorted eigendecompositions with (near-)degenerate eigenvalues.**  A 1-ulp change can
  reorder the principal axes; planar bodies came out tilted by 90 degrees.  See
  `finite-size-particles.md` for the 2d frame rules.
- **Racy lazy caches in `/omp` styles.**  The TIP4P `/omp` pair styles filled a shared
  M-site cache on demand from whichever thread needed it first and published a
  "computed" flag without memory ordering -> torn data on ARM.  Precompute shared data
  before the parallel region (see `openmp-porting.md`).

Test-design fragility (not code bugs), fixed in the test inputs:

- **Near-zero relative comparisons.**  The harness compares relatively, so a reference
  value that is pure roundoff (1e-16 to 1e-21: out-of-plane forces in 2d, the third
  gyration eigenvalue of a planar molecule, velocities of a system started at rest)
  fails with O(1) "errors" wherever the last bits differ.  Make such values exact
  zeros: define `fix enforce2d` AFTER every fix that adds forces or torques, start
  from nonzero velocities, zero roundoff-level eigenvalues.
- **Coordinate-hashed velocities.**  `velocity ... loop geom` hashes coordinate bits, so
  one coordinate differing by 1 ulp (FMA inside `displace_atoms`) changes that atom's
  velocity and, through momentum removal and rescaling, everybody else's.  Use
  `loop all` in test inputs.
- Iterative solvers converged to roundoff (QEq/ReaxFF charge tolerances near 1e-20)
  amplify last-bit differences; such tests need a looser `epsilon`, not bit-level
  references.

Cheap fragility probe: perturb the initial state by one ulp through an untracked copy
of the coefficient or data file (NOT via `post_commands`, which are re-applied after
the restart leg) and rerun the driver with `-s`; a robust test barely changes.

## GPU package debugging

- OpenCL error -48 (`CL_INVALID_KERNEL`, reported from `geryon/ocl_kernel.h`) means the
  kernel program FAILED TO COMPILE.  GPU builds define `UCL_NO_EXIT`, so the compiler
  log goes only to the screen output, which the unit tests capture.  Run the input
  directly (`lmp -sf gpu -in ...`) or the driver with `-v
  --gtest_filter=PairStyle.gpu` on the affected machine to see it.
- Offline syntax check without the vendor's OpenCL: extract the kernel source string
  from `<build>/gpu/<name>_cl.h` and run `clang -x cl -cl-std=CL1.2 -Xclang
  -finclude-default-header -fsyntax-only --target=spir64 -D_SINGLE_DOUBLE
  -DFAST_MATH=1 -DSHUFFLE_AVAIL=0 -DBLOCK_PAIR=256 -DSIMD_SIZE=32
  -DMAX_SHARED_TYPES=11 -DEVFLAG=1` (exactly one precision define; add
  `-Wimplicit-float-conversion` to find double arguments passed to the float-only
  fast-math helpers).  Vendor compilers differ: AMD silently narrows what NVIDIA
  rejects.
- On shared-memory devices (APUs, integrated GPUs) the GPU package falls back to host
  neighbor lists and aliases some host buffers, so device-neighbor code (e.g. the
  special-bond kernels) and packed-buffer paths are NOT exercised there.  A GPU test
  passing on an APU says nothing about those paths on a discrete GPU.
- Every pair kernel must deliver its results through `store_answers()`, also when it
  computes nothing for an atom; otherwise the host accumulates uninitialized device
  memory, which happens to be zero on some devices and garbage on others.

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

`fun:` patterns match MANGLED C++ symbol names: `fun:amd::guessTlsSize*` never
matches, `fun:_ZN3amd*guessTlsSizeEv*` does.  Do not hand-write patterns -- run the
reproducer once with `valgrind --gen-suppressions=all` and generalize the emitted
blocks.

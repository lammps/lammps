# Testing and Verification Methodology

How to validate LAMMPS changes beyond the standard unit/regression suites: choosing
the right gate for refactors vs re-implementations, force verification, test pitfalls,
the CI jobs and platform-dependent failures, GPU package debugging, testing of MDI
coupling, benchmark construction, and memory-checking conventions.

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

**A minimization that fails with `fix shake` or `fix rattle` is not a bug.**  During a
minimization these fixes replace the constraints by stiff harmonic restraints
(keyword `kbond`), which can stall the line search of `min_style cg` ("linesearch
alpha is zero") or stop `fire` with "Shake atoms ... missing".  Use a smaller `kbond`
in such a test or input instead of reporting or working around the failure; for
rigid bodies `doc/src/minimize.rst` describes relaxation with damped dynamics.

**Portability and silent-skip pitfalls in unit tests.**

- gtest's `ContainsRegex`/`MatchesRegex` are NOT portable: on Windows gtest falls back
  to its own regex engine with no character classes, groups, or alternation, while the
  POSIX path lacks `\d`.  Such a test passes on Linux and macOS and fails only in the
  Windows jobs.  For content checks use `HasSubstr` for plain text or
  `utils::strmatch()` (via the `ASSERT_MATCH` idiom in
  `unittest/commands/test_info.cpp`); note the bundled tiny-regex also has no `{n}`
  repetition.
- CTest reports "Passed" even when EVERY gtest case in the binary skipped, and a
  force-style YAML with a missing prerequisite or `input_coeffs` entry self-skips
  silently (one shipped test never ran for weeks that way).  After adding a test,
  verify from the gtest/junit output that its cases actually EXECUTED, not just that
  ctest went green (see the note in `doc/src/Developer_unittest.rst`).
- End states of iterative algorithms (minimizers, converged solvers) are not portable
  reference data -- assert energies/box dimensions with tolerances, never per-atom
  end-state positions.
- A CI failure confined to one job usually comes from what distinguishes that job
  (see the table of CI jobs below).  The most common case is `LAMMPS_SIZES=bigbig`
  (64-bit `tagint`/`bigint`, builder `gcc-bigbig`); a minimal bigbig build
  (`-D LAMMPS_SIZES=bigbig` plus only the needed packages) reproduces it in a few
  minutes.
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

## Continuous integration jobs

Pull requests and merges into `develop` and `maintenance` are tested by GitHub Actions
(`.github/workflows/`) and by the buildbot installation at https://ci.lammps.org; both
report a status per job on GitHub.  The "Details" link of a buildbot status opens the
build: its `configure` step shows the CMake options and its `ctest` step the test
output.  Summaries and trends of all results, including the nightly regression runs,
are at https://dashboard.lammps.org.

- Buildbot builders for pull requests have the prefix `pr-` and merge the current
  `develop` into the branch before building, so a conflict with `develop` fails the
  build.
- They start by themselves for pull requests of members of the LAMMPS organization
  (for outside collaborators once the pull request is not a draft).  For everybody
  else a maintainer starts them with the label `cpu_unit_tests`.  The GPU builders
  always need the label `gpu_unit_tests`; the macOS, FreeBSD, and Windows builders
  need `extra_unit_tests` unless the author is a member.
- The Linux, macOS, and FreeBSD builders skip the tests whose names start with
  `Python`; the Linux CPU builders also skip `GeturlCommand`, because they have no
  network access while testing.  The two GitHub unit-test workflows on Linux run
  those tests.
- `maintenance` is built by the CPU and VM builders, not by the GPU builders.

What distinguishes the jobs (they build most packages unless noted):

| Job | Configuration | Check first when only this job fails |
|---|---|---|
| `gcc-bigbig` | GCC, `LAMMPS_SIZES=bigbig`, Open MPI | 64-bit `tagint`/`bigint`; no FMA |
| `gcc-fftsingle` | GCC, `FFT_SINGLE=on`, MPICH | single precision FFTs; no FMA |
| `gcc11-alma9` | GCC 11 and the oldest supported CMake on AlmaLinux 9, Open MPI 4.1 | old compiler and library; no FMA |
| `gcc-basic` | only the packages of `basic.cmake` | code or a test that needs another package without saying so |
| `clang-cxx17`, `cpp23-gcc`, `cpp23-clang`, `intel-mkl` | Clang; C++23 with GCC and Clang; Intel compilers with MKL | compiler and C++ standard |
| `kokkos-serial-*`, `kokkos-openmp-*` | KOKKOS in double, mixed, and single precision | reduced precision tolerances; `kokkos-serial-single` and `kokkos-openmp-mixed` have no FMA |
| `kokkos-cuda-*`, `gpu-cuda-double`, `gpu-opencl-mixed` | KOKKOS and GPU package on NVIDIA GPUs | device-only code paths, mixed precision |
| `macos-x86_64` | macOS 12, AppleClang 14, no OpenMP | old compiler; no FMA |
| `macos26-x86_64` | macOS 26, AppleClang 17, OpenMP | no FMA |
| `macos-arm64` | macOS 27 on Apple silicon, AppleClang 21, OpenMP, KOKKOS/OpenMP | FMA; Apple's libm |
| `freebsd-x86_64` | FreeBSD 15, clang with libc++ | libc++; the only job that compiles the bundled linalg library |
| `windows-msvc`, GitHub `compile-msvc.yml` | Visual Studio, packages of `windows.cmake` | MSVC: OpenMP 2.0, the regex engine of gtest; runs the `unstable` tests |
| GitHub `unittest-arm64.yml` | GCC on Linux aarch64 | FMA; plain `char` is unsigned |
| GitHub `unittest-fftw.yaml` | GCC, MPI, FFTW3, KOKKOS/OpenMP | FMA |
| `check-examples` | every example input with `-skiprun` | an input that no longer parses |
| `plugins`, `windows-installer` | plugins with stand-in libraries; MinGW cross-compilation of the Windows installers (`develop` only) | `examples/plugins`, plugin CMake files, `cmake/packaging` |

**FMA and the `unstable` label.**  On Rocky Linux 10 (the buildbot machines) and in
the current Ubuntu image of GitHub, GCC targets x86-64-v3 and therefore contracts
`a*b + c` into fused multiply-add (FMA) instructions by default, as it does on aarch64.
The tests with the CTest label `unstable` are numerically sensitive and their
reference data holds for results without FMA.  The Linux jobs marked "no FMA" above
(compiled with `-ffp-contract=off`, or for x86-64-v2 on AlmaLinux 9), the builders
on GPUs, and the Windows jobs run them; the other jobs use `ctest -LE unstable`.  A
test without that label which fails only in the "no FMA" jobs may have reference data
that was generated with an FMA build.  Even a binary compiled without FMA executes
FMA code inside the math library on a CPU that has it (glibc selects its `exp()`,
`pow()`, `log()` at run time), so last-bit differences between two machines running
the same binary are possible.

## Platform-dependent test failures (ARM64, macOS, FMA)

Linux aarch64 differs from a plain x86-64 build in two ways that matter numerically;
both can be emulated on an x86-64 machine in a separate build directory:

| Linux aarch64 property | x86-64 emulation |
|---|---|
| GCC contracts `a*b + c` into fused multiply-add (FMA) by default | `-D CMAKE_CXX_FLAGS="-mfma -ffp-contract=fast"` |
| plain `char` is UNSIGNED (signed on x86-64 and macOS) | `-D CMAKE_CXX_FLAGS=-funsigned-char` |

For the opposite case, a build without FMA on a machine whose compiler uses it by
default, add `-ffp-contract=off` to the compiler flags.  `macos-arm64` uses FMA as
well, but plain `char` is signed there: emulate it with a clang build
(`cmake/presets/clang.cmake`) and `-D CMAKE_CXX_FLAGS="-mfma -ffp-contract=on"`.  A
failure this does not reproduce can come from Apple's libm (e.g. a last-bit difference
in `sin()`) or the compiler version.  The two macOS jobs on x86-64 do not use FMA, so
a failure there has another cause.

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

## Testing MDI coupling

The MDI package has no unit tests at all (`ctest -R mdi` only matches unrelated
names), so changes to `src/MDI/` are tested with `examples/mdi` (commands in its
`README` and `Run.sh`):

- LAMMPS can be both sides of `fix mdi/qm`: driver and (as a stand-in for a quantum
  code) engine.  Three examples (`aimd`, `snapshot`, `series`) run over MPI, TCP, and
  with the engine as a plugin library (`-mdi "... -method LINK -plugin_path <folder
  with liblammps.so>"`).  The thermo output must be identical to that of the matching
  `in.<name>.alone` input without MDI (the forces of the first snapshot may differ in
  the last digits).
- The LAMMPS engine does not provide the commands that `fix mdi/qmmm` requires
  (`>POTENTIAL_AT_NUCLEI` and `<CHARGES`, or `>NLATTICE` and `>LATTICE`), so that fix
  cannot be tested against LAMMPS.  Use a small stand-in engine instead: a C++ program
  linked to the MDI library that CMake built (`<build>/mdi_build_ext/`) that registers
  those commands and returns a deterministic response, driven by
  `examples/QUANTUM/NWChem/in.water.qmmm`.

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

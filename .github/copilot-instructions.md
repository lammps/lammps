# LAMMPS Instructions for AI Coding Agents

This file is the compact, always-loaded core read by GitHub Copilot (coding/cloud agent,
code review, chat), via an import from `.claude/CLAUDE.md` by Claude Code, and via the
`AGENTS.md` stub by other coding agents.  Detailed,
task-specific guides live in `.github/instructions/` (auto-attached by path patterns) and
`.github/dev-docs/` (read on demand); see the index at the end of this file.

## Code Review

Apply the general contribution requirements in
https://docs.lammps.org/Modify_requirements.html and the programming style in
https://docs.lammps.org/Modify_style.html, plus the path-specific rules in
`.github/instructions/` for the changed files: C++ sources and tests
(`source-code`), documentation (`documentation`), KOKKOS (`kokkos`), the GPU package
(`gpu`), unit tests (`force-style-tests`), example inputs (`regression-tests`), and the
build system (`build-system`).  Do not report GCC's `-Wmaybe-uninitialized` warnings:
in this code they have not pointed to a real problem so far.  Uninitialized data
reported by valgrind or a sanitizer is a real finding.

## Build System

**Always use CMake for new builds; always build out-of-source.**  The `CMakeLists.txt`
is in `cmake/`, NOT the repository root.

```bash
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake \
      -D ENABLE_TESTING=on -D BUILD_SHARED_LIBS=on -D DOWNLOAD_POTENTIALS=off -G Ninja
cmake --build build -j 4
# Executable: build/lmp
```

- Use `-S cmake` (NOT `-S .`); never run cmake or make in the repository root or `src/`.
- Presets are in `cmake/presets/` (`basic.cmake`, `gcc.cmake`, `most.cmake`, ...); they
  can be combined with repeated `-C` options.
- Enable packages with `-D PKG_<NAME>=on` (e.g. `-D PKG_MOLECULE=on`); LAMMPS has 80+
  optional packages in `src/<PACKAGE-NAME>/` directories.
- Use `-D DOWNLOAD_POTENTIALS=off` to avoid network dependence in CI or restricted
  environments.
- With `-D ENABLE_TESTING=on`, ALWAYS also set `-D BUILD_SHARED_LIBS=on`: with the
  default static library, each of the 100+ test executables links its own copy of
  LAMMPS (several hundred MB each with the default RelWithDebInfo debug info), and a
  static testing build can exceed 70 GB even with few packages.
- If MPI is not found, install your distribution's MPI development package or set
  `-D MPI_CXX_COMPILER=mpicxx` explicitly.  LAMMPS uses its bundled KISS FFT by
  default; FFTW3 is optional, not required.
- Build times: basic preset ~3-5 minutes; most packages ~10-15 minutes.
- The legacy GNU make build (`cd src && make serial`) and switching between the two
  build systems: `.github/instructions/build-system.instructions.md`.

## Testing & Validation

**Style checks (run before every commit/PR):**
```bash
cd src && make check            # all checks
cd src && make check-whitespace # most common CI failure
cd src && make fix-whitespace   # auto-fix whitespace
cd src && make fix-permissions  # auto-fix file permissions
```
Further named targets: `make check-homepage` (verifies https://www.lammps.org URLs),
`make check-errordocs`, `make check-fmtlib`.

**Unit tests (CTest; requires `-D ENABLE_TESTING=on` and a completed build):**
```bash
ctest --test-dir build -V                # all tests
ctest --test-dir build -V -R <pattern>   # subset by regex
```

Regression tests (the `examples/` inputs) and the documentation build have their own
guides: `.github/instructions/regression-tests.instructions.md` and
`.github/instructions/documentation.instructions.md`.

## Continuous Integration

Two systems run the tests for pull requests and for merges into `develop` and
`maintenance`; both report their status on GitHub:

- **GitHub Actions** (`.github/workflows/`): style checks, the check for
  variable-length arrays, the legacy make build, and unit tests with MPI and FFTW3,
  on Linux ARM64, and with MSVC on Windows.
- **Buildbot** at https://ci.lammps.org: unit tests on Linux (GCC, Clang, and Intel
  compilers, KOKKOS, GPUs), macOS, FreeBSD, and Windows, plus a check that every
  example input can still be read and set up.  The name of a builder describes its
  configuration, and the "Details" link of a status opens the build with the logs of
  all its steps.

**Debugging CI failures:** style-check -> run the matching `make check-*` target in
`src/` and the corresponding `make fix-*`; build failures -> check for `-S cmake`,
package dependencies, and VLA usage; unit tests -> rerun the single test with
`ctest -V -R <name>` in a build that uses the options of the failing job (its
configuration step lists them); example inputs -> verify the Python environment and
whether example inputs were modified.  A unit-test failure in only ONE job usually
comes from what is special about that job: 64-bit `tagint` (`LAMMPS_SIZES=bigbig`,
builder `gcc-bigbig`), single precision FFTs, FMA instructions used or not, macOS,
ARM64, or Windows.  The testing guide lists the jobs with these properties and the
order of the platform triage (char signedness, OpenMP race, FMA contraction,
near-zero comparisons): follow it before loosening tolerances or tagging a test
`unstable`.

## General Conventions

- **7-bit US-ASCII everywhere** (sources, docs, scripts); Unicode is forbidden
  (security policy) and fails CI.
- **User-facing text** (error messages, docs) must avoid computer-science jargon;
  the audience is researchers, not software engineers.
- **File permissions:** `.cpp`/`.h` must NOT be executable; `.sh`/`.py` scripts SHOULD
  be (checked by `make check-permissions`).
- **Python 3.9 is the oldest supported version:** Python code in the tree (the
  `python/` module, `tools/`, test and example scripts) must run with it, so do not
  use newer syntax or library features (e.g. `match` statements, `X | Y` type hints).
- **Library interface with `-DLAMMPS_BIGBIG`:** the per-atom scatter/gather functions
  (`lammps_gather_atoms()`, `lammps_gather()`, `lammps_scatter*()`, their `_concat`
  and `_subset` variants, and the Python and Fortran methods of the same names) stop
  with an error in builds with 64-bit atom IDs.  Code and tests that call them need a
  fallback or a skip (idioms in `.github/instructions/force-style-tests.instructions.md`).
- Root `README` has no extension; subdirectories may use `.md`.
- C++ coding rules and the steps for adding a new style:
  `.github/instructions/source-code.instructions.md`.

## Development Workflow

- Feature branches; PRs target `develop` (NOT `master` or `release`).  The `develop`
  branch is always kept functional (continuous release model).
- Run `cd src && make check` before committing; watch CI on the PR.
- Never rewrite history that has been pushed: no `git rebase`, `git commit --amend`,
  or `git push --force` on a published branch.  Bring a branch up to date by merging
  `develop` into it.
- A bug found in any style is rarely alone: styles and their accelerator variants are
  created by copy-adapt, so defects propagate in both directions.  After root-causing
  a bug, check the base style, all suffix variants (`/omp`, `/kk`, `/gpu`, `/opt`,
  `/intel`), and sibling styles cloned from the same template for the same code shape,
  and fix all occurrences together.
- A report produced by an automated tool (static analyzer, AI-written issue) is a
  lead, not a finding: reproduce it, look for the same defect in sibling code, and
  describe the fix on its own merits (defect, cause, test).  Refer to such a report
  in a commit message only when a maintainer asks for it.  Input files attached to
  such reports can differ from the example they were derived from in ways that are
  hard to see (e.g. a control character inside a number): compare them first.
- The core developers no longer develop the INTEL package: their own work on it is
  limited to bug fixes and adjustments to API changes.  Do not add or propose new
  `/intel` variants yourself.  New INTEL code from contributors is acceptable, since
  the contributor takes on its maintenance: review it on its merits, and do not ask
  for its removal only because of this rule.
- The PR template contains a mandatory **AI Tools Usage** section whose default text
  states no AI was used; when AI tools generated code, edit that section to disclose it
  honestly (tools and models used, which parts they generated).  This section is the
  ONLY place for AI attribution: do NOT add `Co-Authored-By:`, `Claude-Session:`,
  `Generated with ...`, or similar AI-attribution trailer lines to commit messages or
  PR descriptions, and do NOT name an AI tool or model as contributing author in
  source file headers.  This applies to Claude Code, GitHub Copilot, and any other
  coding agent alike.  The policy details are in `doc/src/Developer_agents.rst`.

## Task-Specific Guides

Path-scoped instructions in `.github/instructions/` are attached automatically when
matching files are touched.  The deep dives in `.github/dev-docs/` are NOT loaded
automatically: read them before starting the corresponding kind of work.

| Working on ... | Read |
|---|---|
| C++ code in `src/`, `lib/`, `unittest/`; adding styles | `.github/instructions/source-code.instructions.md` (auto) |
| build system (`cmake/`, legacy make) | `.github/instructions/build-system.instructions.md` (auto) |
| `src/KOKKOS/` styles (rules, policies) | `.github/instructions/kokkos.instructions.md` (auto) |
| porting a style to KOKKOS | `.github/dev-docs/kokkos-porting-guide.md` + `accelerator-porting-backlog.md` |
| porting/auditing `src/OPENMP/` (`/omp`) styles | `.github/dev-docs/openmp-porting.md` + `accelerator-porting-backlog.md` |
| GPU package code in `src/GPU/`, `lib/gpu/` (rules) | `.github/instructions/gpu.instructions.md` (auto) |
| porting a pair style to the GPU package | `.github/dev-docs/gpu-porting-guide.md` + `accelerator-porting-backlog.md` |
| granular/DEM code or tests | `.github/instructions/granular-tests.instructions.md` (auto) |
| documentation (`doc/`) | `.github/instructions/documentation.instructions.md` (auto) |
| force-style YAML tests (`unittest/`) | `.github/instructions/force-style-tests.instructions.md` (auto) |
| example inputs, regression tests | `.github/instructions/regression-tests.instructions.md` (auto) |
| rRESPA support in a fix | `.github/dev-docs/respa-integration.md` |
| finite-size particles, inertia/angmom | `.github/dev-docs/finite-size-particles.md` |
| new/changed styles: MPI, restart, buffers | `.github/dev-docs/style-implementation-notes.md` |
| refactor validation, CI jobs and platform failures, benchmarks, debugging | `.github/dev-docs/testing-and-verification.md` |
| backporting fixes to the `maintenance` branch | `.github/dev-docs/maintenance-backports.md` |

## Trust These Instructions

These instructions are tested and validated.  Only search for additional information
if a specific command fails, a package has special requirements, or the instructions
appear outdated based on error messages.  For package-specific documentation, build
options, and advanced features, refer to https://docs.lammps.org

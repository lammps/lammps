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
(`source-code`), documentation (`documentation`), KOKKOS (`kokkos`), unit tests
(`force-style-tests`), example inputs (`regression-tests`), and the build system
(`build-system`).

## Build System

**Always use CMake for new builds; always build out-of-source.**  The `CMakeLists.txt`
is in `cmake/`, NOT the repository root.

```bash
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake \
      -D ENABLE_TESTING=on -D DOWNLOAD_POTENTIALS=off -G Ninja
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

**Debugging CI failures:** style-check -> run the matching `make check-*` target in
`src/` and the corresponding `make fix-*`; build failures -> check for `-S cmake`,
package dependencies, and VLA usage; unit tests -> rerun the single test with
`ctest -V -R <name>`; regression tests -> verify the Python environment and whether
example inputs were modified.  A unit-test failure on only ONE Linux CI job is
usually the `LAMMPS_SIZES=bigbig` configuration (64-bit `tagint`) -- reproduce with
a minimal bigbig build first (see the testing guide) before suspecting anything else.
A failure only on ARM64 or macOS: follow the platform triage in the testing guide
(char signedness, OpenMP race, FMA contraction, near-zero comparisons) before
loosening tolerances or tagging the test `unstable`.

## General Conventions

- **7-bit US-ASCII everywhere** (sources, docs, scripts); Unicode is forbidden
  (security policy) and fails CI.
- **User-facing text** (error messages, docs) must avoid computer-science jargon;
  the audience is researchers, not software engineers.
- **File permissions:** `.cpp`/`.h` must NOT be executable; `.sh`/`.py` scripts SHOULD
  be (checked by `make check-permissions`).
- Root `README` has no extension; subdirectories may use `.md`.
- C++ coding rules and the steps for adding a new style:
  `.github/instructions/source-code.instructions.md`.

## Development Workflow

- Feature branches; PRs target `develop` (NOT `master` or `release`).  The `develop`
  branch is always kept functional (continuous release model).
- Run `cd src && make check` before committing; watch CI on the PR.
- A bug found in any style is rarely alone: styles and their accelerator variants are
  created by copy-adapt, so defects propagate in both directions.  After root-causing
  a bug, check the base style, all suffix variants (`/omp`, `/kk`, `/gpu`, `/opt`,
  `/intel`), and sibling styles cloned from the same template for the same code shape,
  and fix all occurrences together.
- The INTEL package is unmaintained: it receives only bug fixes and adjustments to
  API changes.  Do not add or propose new `/intel` variants.
- The PR template contains a mandatory **AI Tools Usage** section whose default text
  states no AI was used; when AI tools generated code, edit that section to disclose it
  honestly.  This section is the ONLY place for AI attribution: do NOT add
  `Co-Authored-By:`, `Claude-Session:`, `Generated with ...`, or similar AI-attribution
  trailer lines to commit messages or PR descriptions.  This applies to Claude Code,
  GitHub Copilot, and any other coding agent alike.

## Task-Specific Guides

Path-scoped instructions in `.github/instructions/` are attached automatically when
matching files are touched.  The deep dives in `.github/dev-docs/` are NOT loaded
automatically: read them before starting the corresponding kind of work.

| Working on ... | Read |
|---|---|
| C++ code in `src/`, `lib/`, `unittest/`; adding styles | `.github/instructions/source-code.instructions.md` (auto) |
| build system (`cmake/`, legacy make) | `.github/instructions/build-system.instructions.md` (auto) |
| `src/KOKKOS/` styles (rules, policies) | `.github/instructions/kokkos.instructions.md` (auto) |
| porting a style to KOKKOS | `.github/dev-docs/kokkos-porting-guide.md` + `kokkos-porting-backlog.md` |
| porting/auditing `src/OPENMP/` (`/omp`) styles | `.github/dev-docs/openmp-porting.md` |
| granular/DEM code or tests | `.github/instructions/granular-tests.instructions.md` (auto) |
| documentation (`doc/`) | `.github/instructions/documentation.instructions.md` (auto) |
| force-style YAML tests (`unittest/`) | `.github/instructions/force-style-tests.instructions.md` (auto) |
| example inputs, regression tests | `.github/instructions/regression-tests.instructions.md` (auto) |
| rRESPA support in a fix | `.github/dev-docs/respa-integration.md` |
| finite-size particles, inertia/angmom | `.github/dev-docs/finite-size-particles.md` |
| new/changed styles: MPI, restart, buffers | `.github/dev-docs/style-implementation-notes.md` |
| refactor validation, platform failures, benchmarks, debugging | `.github/dev-docs/testing-and-verification.md` |

## Trust These Instructions

These instructions are tested and validated.  Only search for additional information
if a specific command fails, a package has special requirements, or the instructions
appear outdated based on error messages.  For package-specific documentation, build
options, and advanced features, refer to https://docs.lammps.org

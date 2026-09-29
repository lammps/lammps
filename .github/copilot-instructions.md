# LAMMPS Instructions for AI Coding Agents

This file is the compact, always-loaded core read by GitHub Copilot (coding/cloud agent,
code review, chat) and, via an import from `.claude/CLAUDE.md`, by Claude Code.  Detailed,
task-specific guides live in `.github/instructions/` (auto-attached by path patterns) and
`.github/dev-docs/` (read on demand); see the index at the end of this file.

## Repository Overview

**LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical
molecular dynamics simulation code for parallel computers: a large, mature C++ code base
(~600MB, ~4,000 C++ files in `src/`) maintained by an international team of developers led by
staff at Sandia National Laboratories, open-source under GPL v2.

**Primary languages:** C++17 (core), C, Fortran, Python (interfaces)
**Build systems:** CMake (primary), traditional Make (legacy, subset of packages)
**Key frameworks:** MPI (parallelization), OpenMP (threading), Kokkos (GPU/many-core)

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

**Traditional Make (legacy):** `cd src && make serial` (or `make mpi`); the executables
are `lmp_serial` / `lmp_mpi`.  Enable/disable packages first with `make yes-<package>` /
`make no-<package>` or preset bundles like `make yes-basic` (MANYBODY, MOLECULE, KSPACE,
RIGID); `make pi` shows package status.  Packages needing external libraries or
downloads are CMake-only.

**Switching build systems:** Make -> CMake requires `make -C src purge` first;
CMake -> Make requires `make -C src clean-all` first.  CMake errors out if it detects
make-generated header files in `src/`.

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
cd build && ctest -V                # all tests
cd build && ctest -V -R <pattern>   # subset by regex
```
Tests live in `unittest/` by category: `c-library/`, `commands/`, `force-styles/`,
`formats/`, `fortran/`, `python/`, `utils/`, `granular/`.

**Regression tests** (CI runs them after code review; local runs rarely needed):
```bash
python3 -m venv testenv && source testenv/bin/activate
pip install numpy pyyaml junit_xml
python3 tools/regression-tests/run_tests.py --lmp-bin=build/lmp \
    --config-file=tools/regression-tests/config_quick.yaml --examples-top-level=examples
```

**Documentation build:** `cd doc && make html` must complete without new warnings;
`make spelling` must not report issues (see the documentation guide in the index below).

## Continuous Integration

GitHub Actions workflows in `.github/workflows/`.  On every PR to `develop`:
`style-check.yml` (coding standards), `unittest-linux.yml` (CTest), and
`quick-regression.yml` (regression subset); see that directory for the further
platform, style, and regression workflows.

**Debugging CI failures:** style-check -> run the matching `make check-*` target in
`src/` and the corresponding `make fix-*`; build failures -> check for `-S cmake`,
package dependencies, and VLA usage; unit tests -> rerun the single test with
`ctest -V -R <name>`; regression tests -> verify the Python environment and whether
example inputs were modified.  A unit-test failure on only ONE Linux CI job is
usually the `LAMMPS_SIZES=bigbig` configuration (64-bit `tagint`) -- reproduce with
a minimal bigbig build first (see the testing guide) before suspecting anything else.

## Repository Structure

The top-level `LAMMPS` class (`src/lammps.h`) owns pointers to all subsystems (`atom`,
`force`, `neighbor`, `comm`, `domain`, `modify`, `update`, `output`, `error`, `memory`).
Almost all physics is implemented as named "styles" inheriting from abstract base
classes (`Pair`, `Fix`, `Compute`, `Bond`, `Angle`, `Dihedral`, `Improper`, `Command`),
mapped to keywords via macros (`PairStyle`, `FixStyle`, ...) in the style headers.

## Coding Standards

- **C++17**; follow `.clang-format` in `src/`; keep code ASCII-only.
- **7-bit US-ASCII everywhere** (sources, docs, scripts); Unicode is forbidden
  (security policy) and fails CI.
- **No variable-length arrays** (checked by CI); use `memory->create()` or
  `std::vector`.
- **No alternative logical-operator tokens:** use `&&`, `||`, `!`, `^` -- never `and`,
  `or`, `not`, `xor` (breaks MSVC).
- **Parenthesize each operand of chained `&&`/`||` conditionals** for readability.
- **No two-trip loops for trivial initialization:** assign pairs directly
  (`xstyle[0] = xstyle[1] = NONE;`) instead of a `for` loop with only two trivial
  trips; keep the loop when the body is substantial (unrolling would duplicate code).
- **String formatting with fmtlib** (`fmt::format()`), not `sprintf`.
- **Error handling:** `error->all()` when all MPI ranks hit the error, `error->one()`
  for a single rank; `error->warning()` prints on every rank, so guard with
  `comm->me == 0` where a single message is wanted.
- **User-facing text** (error messages, docs) must avoid computer-science jargon;
  the audience is researchers, not software engineers.
- **RAII for C resources:** prefer `SafeFilePtr` (`src/safe_pointers.h`) over raw
  `FILE *`/`fopen` when touching such code.
- **`delete[]` before `utils::strdup()`:** when storing a copied name (variable,
  region, group ID, ...) in a class member, always `delete[]` the member immediately
  before re-assigning it with `utils::strdup()` -- even when it is provably still
  `nullptr`.  Static analysis (Coverity) flags the bare assignment as a leak, and the
  idiom is defensive against keywords being parsed twice.
- **MPI stubs:** if a serial build misses an MPI symbol, add it to `src/STUBS/mpi.h`
  instead of special-casing the caller.
- **Block comments:** inside `/* ... */`, an embedded `*/` (e.g. in a glob like
  `gb_*/ga_*`) silently terminates the comment; reword or use `//` comments.
- **File permissions:** `.cpp`/`.h` must NOT be executable; `.sh`/`.py` scripts SHOULD
  be (checked by `make check-permissions`).
- Root `README` has no extension; subdirectories may use `.md`.

## Adding New Styles

1. Place `style_name.cpp`/`.h` in `src/` or the appropriate package directory; use a
   similar existing style as template (see https://docs.lammps.org/Modify_style.html).
2. Add new package files to `src/.gitignore`; add renamed/removed file names to
   `src/Purge.list`.
3. Create/update the matching `doc/src/*.rst` file; new publicly visible commands and
   keywords need `.. versionadded:: TBD` (see the documentation guide).
4. Internal styles (upper-case style names) need no documentation.

## Development Workflow

- Feature branches; PRs target `develop` (NOT `master` or `release`).  The `develop`
  branch is always kept functional (continuous release model).
- Run `cd src && make check` before committing; watch CI on the PR.
- A bug found in any style is rarely alone: styles and their accelerator variants are
  created by copy-adapt, so defects propagate in both directions.  After root-causing
  a bug, check the base style, all suffix variants (`/omp`, `/kk`, `/gpu`, `/opt`,
  `/intel`), and sibling styles cloned from the same template for the same code shape,
  and fix all occurrences together.
- The PR template contains a mandatory **AI Tools Usage** section whose default text
  states no AI was used; when AI tools generated code, edit that section to disclose it
  honestly.  This section is the ONLY place for AI attribution: do NOT add
  `Co-Authored-By:`, `Claude-Session:`, `Generated with ...`, or similar AI-attribution
  trailer lines to commit messages or PR descriptions.  This applies to Claude Code,
  GitHub Copilot, and any other coding agent alike.

## Code Review

When performing a code review, apply the general instructions for
contributions to LAMMPS in
https://docs.lammps.org/Modify_requirements.html and the programming
style instructions in https://docs.lammps.org/Modify_style.html

When performing a code review, check any changes to the documentation
(in the `doc/src/` folder) to be written in American English and with
plain ASCII characters.

When performing a code review, ensure that the documentation for any new
commands or added keywords to existing commands contains a
`.. versionadded:: TBD` directive.  For completely new commands, the
`.. versionadded:: TBD` statement should be added after the
"Description" header.  For new keywords to an existing command, the
statement should be added before the paragraph introducing the new
keyword.  For any modified commands or keywords a `.. versionchanged::
TBD` directive should be included in the documentation.  This does not
apply to internal commands (style names written in upper case) or when
the change only adds an accelerated variant of an existing style (then
add the code letter to the respective `Commands_*.rst` file instead).
Check if any examples use the new or modified commands and whether they
need updating.

When reviewing C++ code, ensure that no alternative tokens are used for
logical operators (`&&` not `and`, `||` not `or`, `!` not `not`, `^` not
`xor`); alternative tokens cause compilation failures with some
compilers, most prominently Microsoft Visual C++.

There should not be any printf() statements or fprintf(screen,...) /
fprintf(logfile,...)  in new code. Those should be either removed or
replaced by utils::logmesg().  C++ iostreams (std::cout, std::cerr)
should be replaced by using C-style stdio, if needed using
utils::print() or use utils::logmesg(); both of which support
std::format style formatting.  Any output statements that were added
for the obvious purpose of aiding in debugging should be removed
entirely.

There should not be any new error messages of the error->all(FLERR,
"Illegal XXX command") kind, instead utils::missing_cmd_args() should be
used or more specific error messages with the error pointer argument to
highlight the location of the error as described in
https://docs.lammps.org/Developer_notes.html#errors-warnings-and-informational-messages
If needed compare with similar code that has already been fully
converted to this new style of error messages.  Also note the use of
utils::errorurl() to directing users to more detailed explanations on
the https://docs.lammps.org/Errors_details.html page.

When parsing text files, there should be no use of strtok(), sscanf(),
atoi(), atof() and similar, but the Tokenizer or ValueTokenizer classes
be used and - where possible - also on of the file reader classes. For
converting arguments to numbers, there are also the utils::numeric(),
utils::inumeric(), utils::bnumeric(), and utils::tnumeric() classes.
Same as with error messages, there are likely code block that have been
modernized and are sufficiently similar to serve as an instructive
example.

Commented out code that was apparently added for debugging purposes or
represents disabled features or unused alternative implementations
should be removed entirely, unless a pull request explicitly explains
that those are placeholders for future added features.

When new files are added to package directories in `src`, make sure they
are added to the `src/.gitignore` file, so that copies made in `src` by
the traditional make build are not accidentally committed.  When files
are renamed or removed in package directories, make sure the old names
are added to `src/Purge.list` so stale copies are removed by `make
purge`.


## Task-Specific Guides

Path-scoped instructions in `.github/instructions/` are attached automatically when
matching files are touched.  The deep dives in `.github/dev-docs/` are NOT loaded
automatically: read them before starting the corresponding kind of work.

| Working on ... | Read |
|---|---|
| `src/KOKKOS/` styles (rules, policies) | `.github/instructions/kokkos.instructions.md` (auto) |
| porting a style to KOKKOS | `.github/dev-docs/kokkos-porting-guide.md` + `kokkos-porting-backlog.md` |
| porting/auditing `src/OPENMP/` (`/omp`) styles | `.github/dev-docs/openmp-porting.md` |
| granular/DEM code or tests | `.github/instructions/granular-tests.instructions.md` (auto) |
| documentation (`doc/`) | `.github/instructions/documentation.instructions.md` (auto) |
| force-style YAML tests (`unittest/`) | `.github/instructions/force-style-tests.instructions.md` (auto) |
| rRESPA support in a fix | `.github/dev-docs/respa-integration.md` |
| finite-size particles, inertia/angmom | `.github/dev-docs/finite-size-particles.md` |
| new/changed styles: MPI, restart, buffers | `.github/dev-docs/style-implementation-notes.md` |
| refactor validation, benchmarks, debugging | `.github/dev-docs/testing-and-verification.md` |

## Trust These Instructions

These instructions are tested and validated.  Only search for additional information
if a specific command fails, a package has special requirements, or the instructions
appear outdated based on error messages.  For package-specific documentation, build
options, and advanced features, refer to https://docs.lammps.org

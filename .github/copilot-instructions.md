# LAMMPS Copilot Instructions

## Repository Overview

**LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics simulation code designed for parallel computers. This is a large, mature C++ codebase (~600MB, ~4,000 C++ files in src/) maintained by an international team of developers lead by staff at Sandia National Laboratories as open-source software under GPL v2.

**Primary Languages:** C++17 (core), C, Fortran, Python (interfaces)
**Build Systems:** CMake (primary, modern), Make (traditional, still supported)
**Key Frameworks:** MPI (parallel execution), OpenMP (threading), Kokkos (performance portability)

## Build System & Workflow

### CMake Build (Recommended)

**ALWAYS use CMake for new builds.** The traditional Make system is maintained and only supports a subset of packages. Thus CMake is the primary build system.

**Basic build sequence:**
```bash
# 1. Create build directory (REQUIRED - out-of-source builds only)
mkdir build

# 2. Configure with CMake
cmake -S cmake -B build -C cmake/presets/basic.cmake

# 3. Build
cmake --build build -j 4

# 4. The executable will be: build/lmp
```

**Important CMake details:**
- CMake configuration files are in `cmake/` directory (NOT at repo root)
- Use `-S cmake` to specify the source directory (this is NOT standard - most projects use `-S .`)
- Presets are in `cmake/presets/` - use `-C` to load them
- Common presets: `basic.cmake`, `gcc.cmake`, `most.cmake`
- Combine presets: `-C cmake/presets/gcc.cmake -C cmake/presets/most.cmake`
- Standard CMake options work: `-D BUILD_SHARED_LIBS=on`, `-D ENABLE_TESTING=on`

**Typical CI build configuration (from workflows):**
```bash
cmake -S cmake -B build \
      -C cmake/presets/gcc.cmake \
      -C cmake/presets/most.cmake \
      -D CMAKE_CXX_COMPILER_LAUNCHER=ccache \
      -D CMAKE_C_COMPILER_LAUNCHER=ccache \
      -D BUILD_SHARED_LIBS=off \
      -D DOWNLOAD_POTENTIALS=off \
      -D ENABLE_TESTING=on \
      -G Ninja
cmake --build build
```

**IMPORTANT:** Use `-D DOWNLOAD_POTENTIALS=off` by default to avoid network dependency issues in CI/restricted environments. Only omit this flag if you specifically need LAMMPS to download potential files during the build.

**Build time:** Basic build: ~3-5 minutes, Full build with most packages: ~10-15 minutes

### Traditional Make Build (Legacy, Still Supported)

```bash
cd src
make serial     # Serial build (no MPI)
make mpi        # MPI parallel build
# Executable will be: lmp_serial or lmp_mpi
```

**Available Make targets in src/:**
- `make serial` or `make mpi` - basic builds
- Machine-specific makefiles in `src/MAKE/MACHINES/`
- Options in `src/MAKE/OPTIONS/` for different compiler/feature combinations

### Package Management

LAMMPS has 80+ optional packages. Packages are in `src/[PACKAGE-NAME]/` directories.

**With CMake:** Use `-D PKG_[NAME]=on` (e.g., `-D PKG_MOLECULE=on`, `-D PKG_PYTHON=on`)

**With Make:** Use `make yes-[package]` or `make no-[package]` before building; also
some presets exist: `make yes-basic` or `make-yes-most`; packages requiring extra libraries
or downloads are only supported by CMake.
```bash
cd src
make yes-basic      # Enable MANYBODY, MOLECULE, KSPACE, and RIGID packages
make yes-openmp     # Enable OPENMP package
make yes-misc       # Enable MISC
make serial         # Then build
```

**View package status:** `cd src && make pi` (shows which packages are installed)

## Testing & Validation

### Unit Tests (via CMake + CTest)

```bash
# Configure with testing enabled
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake -D ENABLE_TESTING=on -G Ninja

# Build
cmake --build build

# Run all tests
cd build && ctest -V

# Note: Tests require the executable to be built first
```

**Test organization** (in `unittest/`):
- `c-library/` - C library interface tests
- `commands/` - Input command tests
- `force-styles/` - Pair, bond, angle, kspace style tests
- `formats/` - File format tests
- `fortran/` - Fortran module tests
- `python/` - Python module tests
- `utils/` - Utility function tests

### Regression Tests

```bash
# Setup Python environment (REQUIRED)
python3 -m venv testenv
source testenv/bin/activate
pip install numpy pyyaml junit_xml

# Run regression tests
python3 tools/regression-tests/run_tests.py \
    --lmp-bin=build/lmp \
    --config-file=tools/regression-tests/config_quick.yaml \
    --examples-top-level=examples
```

### Style/Coding Standard Checks

**ALWAYS run these before submitting PRs:**
```bash
cd src
make check-whitespace    # Check whitespace/formatting issues
make check-permissions   # Check file permissions
make check-homepage      # Verify homepage URLs
make check-errordocs     # Check error documentation
make check-fmtlib        # Check fmtlib formatting library
```

**To auto-fix issues:**
```bash
make fix-whitespace
make fix-permissions
```

**All style checks:**
```bash
cd src && make check
```

## GitHub Workflows / CI

The repository uses GitHub Actions for CI with multiple workflows:

**Pull Request Checks (run on all PRs to develop branch):**
1. **style-check.yml** - Runs coding standard checks (`make check-*` targets in src/)
2. **quick-regression.yml** - Builds with most packages, runs regression tests on modified code
3. **unittest-linux.yml** - Builds with LAMMPS_BIGBIG, runs unit tests via CTest

**Additional Checks:**
- **codeql-analysis.yml** - Security scanning
- **check-vla.yml** - Checks for variable-length arrays (not allowed)
- **check-cpp23.yml** - C++23 compatibility check
- **compile-msvc.yml** - Windows MSVC compilation
- **full-regression.yml** - Complete regression suite (on push to develop)

**Workflow dependencies:**
- Ubuntu latest with `ccache`, `ninja-build`, `libeigen3-dev`, `libcurl4-openssl-dev`, `python3-dev`, `mpi-default-bin`, `mpi-default-dev`
- Python packages: `numpy`, `pyyaml`, `junit_xml`
- Build typically uses Ninja generator for speed

## Key Repository Structure

```
lammps/
├── .github/          # GitHub workflows, templates, CodeQL config
│   ├── workflows/    # CI/CD workflow files (12 workflows)
│   ├── CONTRIBUTING.md
│   └── CODEOWNERS
├── cmake/            # CMake build system (USE -S cmake for CMake!)
│   ├── CMakeLists.txt          # Main CMake file
│   ├── presets/                # CMake preset files
│   ├── Modules/                # CMake modules
│   └── packaging/              # Packaging scripts
├── src/              # Source code (~3,777 C++/H files)
│   ├── main.cpp              # Main entry point
│   ├── lammps.cpp/h          # Main LAMMPS class
│   ├── Makefile              # Traditional make system
│   ├── MAKE/                 # Make configurations
│   │   ├── Makefile.serial   # Serial build
│   │   ├── Makefile.mpi      # MPI build
│   │   ├── OPTIONS/          # Compiler/feature options
│   │   └── MACHINES/         # Machine-specific configs
│   ├── [PACKAGE]/            # 80+ optional package directories
│   │   ├── MOLECULE/         # Molecular systems
│   │   ├── KSPACE/           # Long-range electrostatics
│   │   ├── RIGID/            # Rigid body dynamics
│   │   ├── KOKKOS/           # Kokkos acceleration
│   │   └── ...
│   ├── Package.sh            # Package management script
│   └── .clang-format         # Code formatting rules
├── unittest/         # Unit test suite (CTest-based)
├── examples/         # Example input files
├── bench/            # Benchmark inputs
├── tools/            # Pre/post-processing tools
│   ├── coding_standard/      # Style checking scripts
│   ├── regression-tests/     # Regression test framework
│   └── ...
├── doc/              # Documentation source
├── lib/              # External libraries (colvars, kokkos, etc.)
├── python/           # Python interface
├── potentials/       # Potential files
└── README            # Main readme (not .md!)
```

## Common Pitfalls & Important Notes

### Build Issues

1. **CMake source directory:** Use `-S cmake` NOT `-S .` (CMakeLists.txt is in cmake/, not root)

2. **Out-of-source builds only:** ALWAYS create a separate build directory. Never build in source tree.

3. **Switching between Make and CMake:** If you previously used `make` to build, you MUST run `make -C src purge` before using CMake. CMake will error if it detects make-generated header files. Similarly, run `make clean-all` in src/ before switching from CMake to Make.

4. **Package dependencies:** Some packages require others. CMake will warn you; check console output.

5. **MPI detection:** If MPI is not found, install `mpi-default-dev` or set `MPI_CXX_COMPILER=mpicxx` explicitly.

6. **FFTW not required:** LAMMPS uses KISS FFT by default. FFTW3 is optional.

### Code Style

1. **All source must be ASCII:** Unicode characters are not allowed (security policy).

2. **Whitespace matters:** Run `make check-whitespace` and `make fix-whitespace` before committing.

3. **Use .clang-format:** Code should follow .clang-format rules in src/.

4. **No VLAs:** Variable-length arrays are not allowed (checked by CI).

5. **Documentation:** All new commands or features must be documented.  This does *not* apply to
   *internal* commands which can be recognized by their style name written in upper case.  The
   following applies only to publicly visible commands that have style names in lower case.
   Put `.. versionadded:: TBD` or `.. versionchanged:: TBD` in front of paragraphs documenting
   the new or changed functionality or in front of the "Description" headline for completely
   new commands. The `.. versionadded:: TBD` directive should be used with new features or added
   keywords.  The `.. versionchanged:: TBD` directive should be used when the behavior of a
   keyword changes.  The `TBD` will be manually replaced with the release version string during
   the release preparation.  This does not apply when the change is only adding an accelerated
   version of an existing style.  Instead the corresponding code letter should be added to the
   respective Commands_\*.rst file.  Also the documentation should pass running "make check"
   in the doc folder without any output indicating non-ASCII characters, missing entries, or
   duplicate anchors.

### Testing

1. **Build before test:** CTest requires the executable to be built first. If tests fail to find executable, run `cmake --build build` first.

2. **Python environment:** Regression tests require a virtual environment with numpy, pyyaml, junit_xml.

3. **Test selection:** Unit tests are optional packages. Use `-D ENABLE_TESTING=on` with CMake.

### File Permissions & Naming

1. **Source files:** `.cpp` and `.h` files, no executable permission
2. **Scripts:** `.sh` and `.py` files should have executable permission
3. **README vs README.md:** Root README has no extension; subdirs may use .md

## Development Workflow

1. **Branch:** Work on feature branches, submit PRs to `develop` (NOT `master` or `release`)

2. **Style check first:** Run `cd src && make check` before committing

3. **Build locally:** Test with CMake using gcc.cmake + most.cmake presets to match CI

4. **Test changes:** Run relevant unit tests if touching core code, regression tests if modifying examples

5. **Watch CI:** All PR checks must pass. Review CI logs if failures occur.

6. **Continuous release model:** The `develop` branch is always functional. All changes go through PRs with mandatory CI checks.

## Quick Reference Commands

```bash
# Standard development build and test cycle
mkdir build
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake -D ENABLE_TESTING=on -D DOWNLOAD_POTENTIALS=off
cmake --build build -j 4
cd src && make check  # Style checks
cd ../build && ctest -V  # Unit tests

# Minimal build for quick testing
cmake -S cmake -B build -C cmake/presets/basic.cmake
cmake --build build -j 4

# Add a package
cmake -S cmake -B build -C cmake/presets/basic.cmake -D PKG_MOLECULE=on
cmake --build build -j 4

# Traditional make (if needed)
cd src
make serial  # or 'make mpi'
./lmp_serial -in input_file

# Clean everything
rm -rf build
cd src && make clean-all

# Switch from Make to CMake (purge make-generated files)
cd src && make purge
cd .. && mkdir build && cmake -S cmake -B build -C cmake/presets/basic.cmake
```

## Code Review

When performing a code review, apply the general instructions for contributions
to LAMMPS in https://docs.lammps.org/Modify_requirements.html

When performing a code review, apply the programming style instructions
for LAMMPS in https://docs.lammps.org/Modify_style.html

When performing a code review, check any changes to the documentation
(in the `doc/src/` folder) to be written in American English and with
plain ASCII characters.

When performing a code review, ensure that the documentation for any new
commands or added keywords to existing commands contains a
`.. versionadded:: TBD` directive.  For any modified commands or
keywords a `.. versionchanged:: TBD` directive should be included in the
documentation. Check if any examples use the new or modified commands
and check if they need updating.

When reviewing C++ code, ensure that no alternative tokens are used for
logical operators.  That is, use `&&` instead of `and`, `||` instead of
`or`, `!` instead of 'not', `^` instead of `xor` and so on.  These
alternative tokens are only required for ASCII text in some non
US-English characters sets, but the LAMMPS sources code are *supposed*
be in US-English 7-bit ASCII.  Using alternative tokens causes
compilation failures with some compilers by default, most prominently
Microsoft Visual C++.

When new files are added to package directories in `src`, make sure
they are added to the `src/.gitignore` file, so that the copies in
`src` make by the traditional make build system are not accidentally
added to the git repository.

Whe files are renamed or removed from package directories in `src`,
make sure the old file names are added the `src/Purge.list` file so
their copies from the traditional make build system are properly
removed with "make purge".

## Documentation Changes

When modifying documentation files in `doc/src/`:

**Build and validate documentation:**
```bash
cd doc
make html          # Build HTML, check for warnings
make pdf           # Build PDF (requires pdflatex)
make spelling      # Check spelling
make anchor_check  # Check for duplicate anchors
make style_check   # Verify style lists are complete
```

Ensure that building the documentation with "make html" and "make pdf"
can complete and does *NOT* produce any errors or warnings about sphinx
or docutils syntax issues or incorrect references.

Also make certain that "make spelling" does not report in any spelling
issues; those need to be either remedied or exceptions (e.g. for code
examples or author names of cited references) should be added to the
file "doc/utils/sphinx-config/false-positives.txt".

**Documentation conventions:**
- Use reStructuredText format (`.rst` files)
- Use American English spelling
- Use ASCII characters only
- Wrap code examples in `.. code-block::` with appropriate language (LAMMPS, bash, c++, python)
- Use `.. note::` for important remarks and `.. warning::` for critical warnings
- New commands require `.. versionadded:: TBD`
- Modified commands require `.. versionchanged:: TBD`

## Debugging CI Failures

When a CI check fails, diagnose using these steps:

**1. Style check failures (`style-check.yml`):**
```bash
cd src
make check-whitespace    # Most common - fix with: make fix-whitespace
make check-permissions   # Fix with: make fix-permissions
make check-homepage      # Verify https://www.lammps.org URLs
make check-errordocs     # Check error documentation
make check-fmtlib        # Verify fmtlib formatting
```

**2. Build failures:**
- Check CMake output for missing dependencies
- Ensure `-S cmake` (not `-S .`) is used
- Verify package dependencies are met
- Check for VLA (variable-length array) usage - not allowed

**3. Unit test failures:**
- Run specific failing test: `cd build && ctest -V -R <test_name>`
- Check if test requires specific packages to be enabled
- Verify the executable was built before running tests

**4. Regression test failures:**
- Ensure Python environment has numpy, pyyaml, junit_xml
- Check if example inputs were modified correctly

## Short-Circuit Instructions

**STOP and check these common mistakes:**

1. **Wrong CMake source directory:**
   - WRONG: `cmake -S . -B build`
   - CORRECT: `cmake -S cmake -B build`

2. **Building in source tree:**
   - NEVER run cmake or make in the repository root
   - ALWAYS create a separate `build/` directory

3. **Mixed build systems:**
   - If switching from Make to CMake: run `make -C src purge` first
   - If switching from CMake to Make: run `make -C src clean-all` first

4. **Unicode in source files:**
   - All source code must be ASCII only
   - Unicode characters will cause CI to fail

5. **Missing whitespace fixes:**
   - Always run `cd src && make fix-whitespace` before committing

6. **Incorrect file permissions:**
   - `.cpp` and `.h` files must NOT be executable
   - `.sh` and `.py` scripts SHOULD be executable

## Sample Prompts

**Adding a new pair style:**
> "Create a new pair style called `pair_example` that implements [description]. Follow the pattern in `src/pair_lj_cut.cpp` and add documentation in `doc/src/pair_example.rst`."

**Fixing a bug in a compute:**
> "Fix bug in `compute_temp.cpp` where [description]. Add a unit test in `unittest/` to prevent regression."

**Adding a new package:**
> "Create a new package called MYPACKAGE with [features]. Include CMakeLists.txt entries, documentation, and example inputs."

**Updating documentation:**
> "Update the documentation for `fix_nve.rst` to include the new `keyword` option. Use `.. versionchanged:: TBD` directive."

**Debugging build failure:**
> "The CI build is failing with [error]. Diagnose and fix the issue."

# Porting Pair Styles to KOKKOS

Below is a structured plan and reference for porting some
pair styles to the KOKKOS package.  It covers the required code changes,
documentation updates, and groups the work into Copilot-session-sized
batches ordered from simplest to most complex.

---

## Background and Motivation

The KOKKOS package provides GPU-accelerated (CUDA/HIP) and many-core CPU
(OpenMP/threads) variants of LAMMPS styles using Kokkos abstractions.
Porting a pair style to KOKKOS means:

1. Creating `src/KOKKOS/pair_<name>_kokkos.h` and
   `src/KOKKOS/pair_<name>_kokkos.cpp`.
2. Registering the new files in `src/KOKKOS/Install.sh`.
3. Making small changes to the **base class** header (adding `virtual` to
   `allocate()` and adding `if (copymode) return;` to the destructor of the
   base class if not already present).
4. Updating the RST documentation in `doc/src/pair_<name>.rst`.
5. Updating `doc/src/Commands_pair.rst`.

---

## Key Rules for KOKKOS Porting

### 1.  `if (copymode) return;` in the base-class destructor

When a Kokkos pair style is copied for use on-device (via `copymode = 1`),
the copy's destructor must not free memory that belongs to the original object.
Every **base class** that a `_kokkos` style will inherit from must have:

```cpp
PairFoo::~PairFoo()
{
  if (copymode) return;   // MUST be the first line
  // ... rest of destructor ...
}
```

If the base-class destructor already has this line, no change is needed.
Check with:

```bash
grep -n "copymode" src/<PKG>/pair_<name>.cpp
```

### 2.  `virtual void allocate()` in the base-class header

The KOKKOS style overrides `allocate()` to replace plain arrays with dual
views.  The override only works correctly in C++ when the base-class
declaration is `virtual`.  Change:

```cpp
void allocate();          // old
virtual void allocate();  // required
```

Check with:

```bash
grep -n "allocate" src/<PKG>/pair_<name>.h
```

### 3.  Per-atom view members: `typename AT::` not `DAT::`

Every member that stores a per-atom view obtained from
`atomKK->k_<field>.view<DeviceType>()` **must** be declared with the
template-parameterized alias `typename AT::t_<type>`, **never** with the
hardcoded alias `DAT::t_<type>`.

```cpp
// WRONG — kk/host instantiation fails with Kokkos ViewMapping assertion:
DAT::t_tagint_1d_randomread molecule;

// CORRECT — adapts to DeviceType (LMPDeviceType or LMPHostType):
typename AT::t_tagint_1d_randomread molecule;
```

`AT` is `ArrayTypes<DeviceType>`.  `DAT` is `ArrayTypes<LMPDeviceType>`.
In GPU builds, the pair style is instantiated for *both* `LMPDeviceType` and
`LMPHostType`.  Views with `DAT::` types live in GPU memory; assigning a
host-side `atomKK` view to them triggers a static assertion in Kokkos.

**Safe to use `DAT::`:** dual views that the pair style allocates itself and
never assigns from an `atomKK` view:
`k_eatom`, `k_vatom`, `k_cutsq`, `k_cut_ljsq`, `k_cut_coulsq`.

**Must use `typename AT::`:** all per-atom views set via
`atomKK->k_<field>.view<DeviceType>()`:
`x`, `c_x`, `f`, `type`, `q`, `molecule`, `radius`, etc.

Quick audit command:
```bash
grep "DAT::t_" src/KOKKOS/pair_<name>_kokkos.h
```
Every hit should be a dual-view (type contains `tdual` or `ttransform`), not
a plain per-atom view (type contains `t_`).

### 4.  `src/KOKKOS/Install.sh` entries

Add two lines per new file:

```sh
action pair_<name>_kokkos.cpp pair_<name>.cpp   # if base is in another package
action pair_<name>_kokkos.h   pair_<name>.h     # idem
```

If the base class is in the **core** (`src/`), omit the second argument:

```sh
action pair_<name>_kokkos.cpp
action pair_<name>_kokkos.h
```

### 5.  Documentation changes

**In `doc/src/pair_<name>.rst`:**

* Add `.. index:: pair_style <name>/kk` near the top.
* Add `*<name>/kk*` to the `Accelerator Variants:` line.

**In `doc/src/Commands_pair.rst`:**

* Add `k` to the accelerator letter string in the entry for this style
  (letters: `g`=GPU, `i`=INTEL, `k`=KOKKOS, `o`=OPENMP, `t`=OPT).

---

## Standard KOKKOS Template (pair_kokkos.h framework)

Most simple pairwise styles (no three-body, no many-body embedding) fit the
`pair_kokkos.h` template.  The KOKKOS class must implement:

| Function | COUL_FLAG=0 | COUL_FLAG=1 |
|---|---|---|
| `compute_fpair` | required | required (VdW part) |
| `compute_evdwl` | required | required |
| `compute_fcoul` | return 0 | required |
| `compute_ecoul` | return 0 | required |

The template automatically handles half/half-thread/full neighbor list
dispatch, reduction of forces, and energy/virial accumulation.

See `src/KOKKOS/pair_born_kokkos.{h,cpp}` and
`src/KOKKOS/pair_lj_cut_kokkos.{h,cpp}` as canonical examples.

### Typical header skeleton (COUL_FLAG=0)

```cpp
#ifdef PAIR_CLASS
PairStyle(<name>/kk,     Pair<Name>Kokkos<LMPDeviceType>);
PairStyle(<name>/kk/device,Pair<Name>Kokkos<LMPDeviceType>);
PairStyle(<name>/kk/host,  Pair<Name>Kokkos<LMPHostType>);
#else
#ifndef LMP_PAIR_<NAME>_KOKKOS_H
#define LMP_PAIR_<NAME>_KOKKOS_H
#include "pair_kokkos.h"
#include "pair_<name>.h"
#include "neigh_list_kokkos.h"
namespace LAMMPS_NS {
template<class DeviceType>
class Pair<Name>Kokkos : public Pair<Name> {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  Pair<Name>Kokkos(class LAMMPS *);
  ~Pair<Name>Kokkos() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_<name> {
    KOKKOS_INLINE_FUNCTION params_<name>() { /* zero all */ }
    KOKKOS_INLINE_FUNCTION params_<name>(int) { /* zero all */ }
    KK_FLOAT cutsq, /* ... pair-specific params ... */;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(...) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(...) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(...) const { return 0; }

  Kokkos::DualView<params_<name>**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_<name>**,...>::t_dev_const_um params;
  params_<name> m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  Kokkos::DualView<KK_FLOAT**,DeviceType> k_cutsq;
  typename Kokkos::DualView<KK_FLOAT**,...>::t_dev d_cutsq;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  Kokkos::DualView<KK_FLOAT*,DeviceType> k_eatom, k_vatom;
  int neighflag, newton_pair;
  int nlocal, nall, eflag, vflag;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,FULL,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,FULL,false>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALF,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALF,false>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALFTHREAD,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    FULL,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    HALF,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    HALFTHREAD,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<Pair<Name>Kokkos<DeviceType>,
    CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<Pair<Name>Kokkos<DeviceType>>(
    Pair<Name>Kokkos<DeviceType>*);
};
}
#endif
#endif
```

For `COUL_FLAG=1`, add `compute_fcoul` and `compute_ecoul`, add
`special_coul`, `qqrd2e`, `q` array, and change `CoulLongTable<0>` to
`CoulLongTable<1>` in the friend declarations.

### Typical `.cpp` skeleton

```cpp
template<class DeviceType>
Pair<Name>Kokkos<DeviceType>::Pair<Name>Kokkos(LAMMPS *lmp) : Pair<Name>(lmp)
{
  respa_enable = 0;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

template<class DeviceType>
Pair<Name>Kokkos<DeviceType>::~Pair<Name>Kokkos()
{
  if (copymode) return;
  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in; vflag = vflag_in;
  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);
  if (eflag_atom) { /* reallocate k_eatom */ }
  if (vflag_atom) { /* reallocate k_vatom */ }
  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);
  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal; nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0..3] = ...;
  NeighListKokkos<DeviceType>* k_list = ...;
  ...
  EV_FLOAT ev = pair_compute<Pair<Name>Kokkos<DeviceType>,
    CoulLongTable<COUL_FLAG>>(this, k_list);
  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) { virial[0..5] += ev.v[0..5]; }
  if (eflag_atom) { k_eatom.template modify<DeviceType>(); k_eatom.sync_host(); }
  if (vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }
  if (vflag_fdotr) pair_virial_fdotr_compute(this);
  copymode = 0;
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::allocate()
{
  Pair<Name>::allocate();
  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_<name>**,...>("Pair<Name>::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::init_style()
{
  Pair<Name>::init_style();
  // Request full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
}

template<class DeviceType>
double Pair<Name>Kokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = Pair<Name>::init_one(i,j);
  k_params.view_host()(i,j).<field> = static_cast<KK_FLOAT>(<field>[i][j]);
  /* ... fill all params fields ... */
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i < MAX_TYPES_STACKPARAMS+1 && j < MAX_TYPES_STACKPARAMS+1)
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
  m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();
  return cutone;
}

// Explicit instantiation at the bottom of the .cpp file
namespace LAMMPS_NS {
template class Pair<Name>Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class Pair<Name>Kokkos<LMPHostType>;
#endif
}
```

---

## Candidate Pair Styles to Port

52 pair styles are possible candidates for porting to KOKKOS.  They are
grouped below by porting complexity, from easiest to hardest.

---

## Group 1 — Simple Pairwise, No Coulomb (EXTRA-PAIR; 13 styles)

**Status: COMPLETED**

**Complexity:** Low.  Use the `pair_kokkos.h` template with `COUL_FLAG=0`.
Implement `compute_fpair` and `compute_evdwl`; `compute_ecoul` returns 0.
Pattern to follow: `pair_born_kokkos` or `pair_beck_kokkos`.

**Session size recommendation:** 4–6 styles per session.

| Pair style | Package | Base header | Status |
|---|---|---|---|
| `mie/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_mie_cut.h` | **done** |
| `nm/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_nm_cut.h` | **done** |
| `gauss/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_gauss_cut.h` | **done** |
| `harmonic/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_harmonic_cut.h` | **done** |
| `cosine/squared` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_cosine_squared.h` | **done** |
| `lj/smooth` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lj_smooth.h` | **done** |
| `lj/smooth/linear` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lj_smooth_linear.h` | **done** |
| `morse/smooth/linear` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_morse_smooth_linear.h` | **done** |
| `ufm` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_ufm.h` | **done** |
| `wf/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_wf_cut.h` | **done** |
| `born/gauss` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_born_gauss.h` | **done** |
| `lj/mdf` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lj_mdf.h` | **done** |
| `lennard/mdf` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lennard_mdf.h` | **done** |

**Required base-class changes:**
- Add `virtual void allocate();` in each header (if not already `virtual`).
- Add `if (copymode) return;` as first line of each destructor.

**Documentation:** Add `.. index:: pair_style <name>/kk` and add `*<name>/kk*`
to `Accelerator Variants`.  Add `k` to each entry in `Commands_pair.rst`.

---

## Group 2 — Simple Pairwise, No Coulomb (continued; 5 styles)

**Status: COMPLETED**

**Complexity:** Low.  Same as Group 1 but with more parameters per pair-type.

| Pair style | Package | Notes | Status |
|---|---|---|---|
| `momb` | EXTRA-PAIR | London + Buckingham, extra exponent params | **done** |
| `buck/mdf` | EXTRA-PAIR | Buckingham with switching, see `pair_buck_kokkos` | **done** |
| `pedone` | EXTRA-PAIR | Morse-like with r^-6 correction, no Coulomb | **done** |
| `lj/pirani` | EXTRA-PAIR | Generalized LJ with Pirani potential, per-type m,n | **done** |
| `ylz` | ASPHERE | Orientation-dependent potential, requires custom kernel | **done** |

**Note:** `ylz` is in the ASPHERE package but computes only short-range VdW
with fixed charge; it should use `COUL_FLAG=0`.

---

## Group 3 — Pairwise with Short-Range Coulomb (Wolf/DSF/cut; 8 styles)

**Status: COMPLETED**

**Complexity:** Moderate.  Use the `pair_kokkos.h` template with `COUL_FLAG=1`.
Add `compute_fcoul` and `compute_ecoul`.  The Coulomb part uses no long-range
tables (Wolf damping, DSF cutoff, or plain cut).
Pattern to follow: `pair_lj_cut_coul_dsf_kokkos`.

| Pair style | Package | Coulomb type | Status |
|---|---|---|---|
| `born/coul/wolf` | EXTRA-PAIR | Wolf damped | **done** |
| `lj/cut/coul/wolf` | EXTRA-PAIR | Wolf damped | **done** |
| `nm/cut/coul/cut` | EXTRA-PAIR | cut | **done** |
| `coul/diel` | EXTRA-PAIR | dielectric screening | **done** |
| `coul/shield` | INTERLAYER | exponential screening | **done** |
| `buck6d/coul/gauss/dsf` | MOFFF | Gaussian-damped DSF | **done** |
| `coul/cut/global` | EXTRA-PAIR | tiny wrapper of `coul/cut` | **done** |
| `lj/cut/sphere` | EXTRA-PAIR | radius from per-atom property; no Coulomb (COUL_FLAG=0) | **done** |

**Notes:**
- `coul/diel`: stores scalar dielectric parameters (`a_eps`, `b_eps`, `eps_s`) as `KK_FLOAT`
  members and copies them before each kernel launch.
- `coul/shield`: includes per-atom `molecule` view to skip same-layer pairs; the taper
  function (Tap, dTap) is inlined using Horner's method with hardcoded coefficients.
  **Important:** the `molecule` member must be declared `typename AT::t_tagint_1d_randomread`
  (not `DAT::t_tagint_1d_randomread`), otherwise the `kk/host` instantiation fails to
  compile in GPU builds.  See *Lessons Learned* for the general rule.
- `coul/cut/global`: inherits `PairCoulCutKokkos<DeviceType>` directly; only overrides
  `coeff()` (enforces 2 args) and `extract()`.
- `lj/cut/sphere`: COUL_FLAG=0; the per-atom `radius` view (declared `typename AT::
  t_kkfloat_1d_randomread`) is synced; cutsq stores the maximum possible cutoff per
  type pair; the actual per-atom cutoff and sigma are computed inside
  `compute_fpair`/`compute_evdwl`.

---

## Group 4 — Pairwise with Long-Range Coulomb (Ewald/PPPM; 4 styles)

**Status: COMPLETED**

**Complexity:** Moderate.  `COUL_FLAG=1`, use `init_tables()` for Coulomb.
Pattern to follow: `pair_born_kokkos` extended with `pair_lj_cut_coul_long_kokkos`.

| Pair style | Package | Notes | Status |
|---|---|---|---|
| `nm/cut/coul/long` | EXTRA-PAIR | NM with long-range Ewald | **done** |
| `lj/switch3/coulgauss/long` | YAFF | YAFF + Gaussian charge Ewald | **done** |
| `mm3/switch3/coulgauss/long` | YAFF | YAFF MM3 variant | **done** |
| `buck6d/coul/gauss/long` | MOFFF | MOFFF Buckingham + Ewald | **done** |

---

## Group 5 — TIP4P Water Models (4 styles)

**Complexity:** High.  The TIP4P geometry correction (virtual oxygen site M)
requires special handling in the neighbor loop.  Study the existing
`pair_lj_cut_tip4p_long_kokkos` carefully before attempting.

| Pair style | Package | Notes |
|---|---|---|
| `tip4p/cut` | MOLECULE | TIP4P cut |
| `tip4p/long` | KSPACE | TIP4P long Coulomb |
| `lj/cut/tip4p/cut` | MOLECULE | LJ + TIP4P cut |
| `lj/cut/tip4p/long` | KSPACE | LJ + TIP4P long |

---

## Group 6 — Granular, Colloidal, Lubrication (2 styles)

**Complexity:** High.  These involve contact mechanics (overlap detection,
friction history) or hydrodynamic lubrication requiring per-atom state.

| Pair style | Package | Notes |
|---|---|---|
| `gran/hooke` | GRANULAR | frictionless Hertz; no history arrays needed |
| `gran/hertz/history` | GRANULAR | requires shear history Kokkos view |

---

## Group 7 — MANYBODY Three-body (SW-like, Tersoff tables; 3 styles)

**Complexity:** High.  These use a short neighbor list in addition to the
regular neighbor list.  Follow `pair_sw_kokkos` and `pair_tersoff_kokkos`
as patterns; the short-neigh build kernel and the three-body loop kernel
must be implemented manually (the `pair_kokkos.h` template does not support
three-body).

| Pair style | Package | Notes |
|---|---|---|
| `sw/mod` | MANYBODY | SW variant; inherits from PairSW; follow `pair_sw_kokkos` |
| `nb3b/harmonic` | MANYBODY | three-body harmonic, independent from SW hierarchy |
| `tersoff/mod/c` | MANYBODY | Tersoff mod + carbon correction |

---

## Group 8 — Moderate Many-Body / Geometry-Dependent (6 styles)

**Complexity:** High.  These require non-trivial per-atom or angle-dependent
loops that do not map onto the `pair_kokkos.h` template.

| Pair style | Package | Notes |
|---|---|---|
| `atm` | MANYBODY | three-body Axilrod-Teller-Muto; triple loop |
| `local/density` | MANYBODY | local density embedding; two-pass loop |
| `gayberne` | ASPHERE | ellipsoid-ellipsoid; per-atom quaternions on device |
| `resquared` | ASPHERE | re-squared ellipsoid; similar to gayberne |
| `line/lj` | ASPHERE | line-segment LJ; per-atom shape/orientation |
| `tri/lj` | ASPHERE | triangle LJ; per-atom shape/orientation |

---

## Group 9 — H-Bond, Dipole, Peridynamics (6 styles)

**Complexity:** High.  Angle-dependent (H-bond) or orientation-dependent
(dipole) or continuum mechanics (peri) interactions.

| Pair style | Package | Notes |
|---|---|---|
| `hbond/dreiding/lj` | MOLECULE | angle-dependent H-bond; triple loop |
| `hbond/dreiding/morse` | MOLECULE | Morse variant of above |
| `hbond/dreiding/lj/angleoffset` | EXTRA-MOLECULE | angle-offset extension |
| `hbond/dreiding/morse/angleoffset` | EXTRA-MOLECULE | angle-offset extension |
| `lj/sf/dipole/sf` | DIPOLE | shifted-force dipole + LJ |
| `lj/expand/sphere` | EXTRA-PAIR | per-atom sphere radii required |

---

## Documentation Checklist per Style

For each pair style ported, update the following:

1. **`src/KOKKOS/pair_<name>_kokkos.h`** — new file
2. **`src/KOKKOS/pair_<name>_kokkos.cpp`** — new file
3. **`src/KOKKOS/Install.sh`** — add `action` lines
4. **Base class `.h` in `src/<PKG>/pair_<name>.h`**:
   - Make `allocate()` virtual
5. **Base class `.cpp` in `src/<PKG>/pair_<name>.cpp`**:
   - Add `if (copymode) return;` as first line of destructor
6. **`doc/src/pair_<name>.rst`**:
   - Add `.. index:: pair_style <name>/kk` near the top and so that they are ordered in groups by the base pair style and then for the same base style alphabetically according to the accelerator suffix. So `/kk` follows `/gpu` but comes before `/omp`.
   - Add `*<name>/kk*` to `Accelerator Variants:`
   - Order the list of accelerator variants alphabetically
7. **`doc/src/Commands_pair.rst`**:
   - Add letter `k` to the entry for this style
   - When adding to existing letter use (ko) and *not* (o,k), that is order alphabetically and do not use a comma.

---

## Build and Test Instructions

### Build for testing

```bash
mkdir -p build
cmake -S cmake -B build \
  -C cmake/presets/gcc.cmake \
  -C cmake/presets/most.cmake \
  -D PKG_KOKKOS=on \
  -D Kokkos_ENABLE_OPENMP=on \
  -D DOWNLOAD_POTENTIALS=off \
  -D ENABLE_TESTING=on \
  -G Ninja
cmake --build build -j 4
```

### Style checks before committing

```bash
cd src && make check-whitespace && make check-permissions
```

### Running unit tests

```bash
cd build && ctest -V -R pair
```

---

## Summary of Common Mistakes to Avoid

1. **Forgetting `if (copymode) return;`** — causes double-free crash on GPU.
2. **Forgetting `virtual void allocate()`** — silently calls the wrong allocate,
   missing KK dual-view setup, leading to incorrect results or crashes.
3. **Wrong `CoulLongTable` template parameter** — use `<0>` when
   `COUL_FLAG=0`, `<1>` when `COUL_FLAG=1`.
4. **Missing friend declarations** — `pair_compute*` helpers are templated
   free functions; the class must declare them as friends.
5. **Not using `KK_FLOAT`** — all floating-point fields in `params_*` structs
   must use `KK_FLOAT` (not `double`).
6. **Not adding `// NOLINTNEXTLINE` before each `KOKKOS_INLINE_FUNCTION`** —
   required to suppress linter warnings in LAMMPS.
7. **Wrong Install.sh action format** — when base class is in another package
   (e.g., EXTRA-PAIR), pass both filenames; when base is in `src/`, omit
   second filename.
8. **Not rebuilding with `make purge` after switching build systems.**
9. **Using `DAT::t_<type>` for per-atom views assigned from `atomKK`** —
   every per-atom view assigned via `atomKK->k_<field>.view<DeviceType>()`
   must be declared `typename AT::t_<type>` (not `DAT::t_<type>`).  `AT` is
   `ArrayTypes<DeviceType>`, so the view type adapts to the template parameter.
   `DAT` is hardcoded to `ArrayTypes<LMPDeviceType>`, which compiles fine for
   the `kk/device` instantiation but triggers a Kokkos `ViewMapping` static
   assertion when the class is instantiated for `LMPHostType` (the `kk/host`
   variant required in GPU builds).  Standard views that hold their own data
   (e.g., `DAT::ttransform_kkacc_1d k_eatom`, `DAT::ttransform_kkfloat_2d
   k_cutsq`) are correctly declared with `DAT::` because they are allocated
   locally and are not assigned from an `atomKK->k_*.view<DeviceType>()` call.

---

## Lessons Learned (updated after Group 2)

### Global scalar parameters in pair styles

Some pair styles (e.g., `momb`) have global scalar parameters set by
`settings()` rather than per-pair parameters set by `pair_coeff`.  These
scalars are stored as plain `double` members of the base class (e.g.,
`sscale`, `dscale`).  For the KOKKOS version, copy them to `KK_FLOAT`
members (e.g., `m_sscale`, `m_dscale`) in `compute()` before the kernel
launch.  They can then be captured by the `KOKKOS_INLINE_FUNCTION` methods
through `this` (which is captured by the Kokkos functor).

### MDF (Mei-Davenport-Fernando) tapering

The `buck/mdf` style has a smooth taper applied when the inter-atomic
distance is in the range `[cut_inner, cut_outer]`.  In the standard
`pair_kokkos.h` template, the `compute_fpair` function receives `rsq` and
must compute `r = sqrt(rsq)` itself.  The tapering logic (computing `d`,
`tt`, `dt`) can be inlined directly in `compute_fpair` and `compute_evdwl`
using a conditional `if (rsq > cut_inner_sq)`.  Store `cut_inner`,
`cut_inner_sq`, and the outer `cutsq` in the per-pair `params_*` struct so
the kernel can branch correctly.

### Orientation-dependent pair styles (ylz, gayberne, etc.)

The `ylz` style computes orientation-dependent forces and torques based on
per-atom quaternion data.  It does not fit the `pair_kokkos.h` template
because:
1. It needs to read quaternion data from `AtomVecEllipsoidKokkos`.
2. It accumulates torques as well as forces.
3. The inner loop is not a simple `compute_fpair` scalar.

For such styles a **custom kernel** is required, following the pattern of
`pair_lj_cut_dipole_cut_kokkos`:

- Use `TagPairXXXKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>` tags.
- Override `compute()` to launch the appropriate template instantiation with
  `Kokkos::parallel_for` or `Kokkos::parallel_reduce`.
- Accumulate torques with atomic operations via a view with
  `Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>`.
- Access the ellipsoid bonus data (quaternions) via
  `avecKK->k_bonus.view<DeviceType>()` where `avecKK` is a pointer to
  `AtomVecEllipsoidKokkos` cast from `atom->style_match("ellipsoid")`.
- Use `MathExtraKokkos::quat_to_mat()` to convert quaternion to rotation
  matrix on-device.
- Set `datamask_read` to include `TORQUE_MASK | ELLIPSOID_MASK`.
- Set `datamask_modify` to include `TORQUE_MASK`.
- Store the neighbor-list arrays explicitly (`d_neighbors`, `d_numneigh`,
  `d_ilist`) rather than relying on the `pair_compute` helper.

### lj/pirani: Kokkos::pow in device code

The ILJ/Pirani potential uses `pow(rx, n(rx))` where `n(rx)` is itself a
function of the inter-atomic distance.  `Kokkos::pow(base, exp)` is the
device-portable equivalent of `std::pow`.  Using `std::pow` inside a
`KOKKOS_INLINE_FUNCTION` causes a compiler error on GPU targets.  Always
prefer `Kokkos::` math functions in device kernels:
`Kokkos::sqrt`, `Kokkos::exp`, `Kokkos::pow`, `Kokkos::log`,
`Kokkos::sin`, `Kokkos::cos`, etc.

### Force sign convention

Be careful with force sign conventions.  In the base CPU code, forces are
often accumulated as `f[i][0] += delx * fpair` where `fpair = F/r` and
`delx = xj - xi`.  The KOKKOS kernel must use the same convention.  For the
`ylz` style the force on atom `i` due to atom `j` is `-dU/dr * rhat` where
`rhat` points from `i` to `j` (i.e., `r12 = xj - xi`); the force on `j`
is the Newton-3 partner `-force_on_i`.

---

## Lessons Learned (updated after Group 3, session May 2026)

### Per-atom views: always use `typename AT::` not `DAT::`

The most dangerous portability mistake for COUL_FLAG=1 styles (and any style
accessing extra per-atom arrays) is declaring a per-atom view with the
hardcoded `DAT::` type alias instead of the template-parameterized `AT::`:

```cpp
// WRONG — compiles for kk/device but fails at kk/host instantiation:
DAT::t_tagint_1d_randomread molecule;

// CORRECT — adapts to DeviceType, works for both kk/device and kk/host:
typename AT::t_tagint_1d_randomread molecule;
```

The rule is simple: any member that is assigned from
`atomKK->k_<field>.view<DeviceType>()` must be declared with `typename AT::`.
Views that manage their own data (dual views allocated in `allocate()`, e.g.,
`k_eatom`, `k_vatom`, `k_cutsq`) are correctly declared with `DAT::` because
they do not hold per-atom data obtained from `atomKK`.

This bug manifests as a Kokkos `ViewMapping` static assertion at compile time
when targeting GPU architectures (because GPU builds instantiate both
`Pair<X>Kokkos<LMPDeviceType>` and `Pair<X>Kokkos<LMPHostType>`).  On CPU-only
builds both `LMPDeviceType` and `LMPHostType` share the same memory space, so
the bug is silent.

The `pair_coul_shield_kokkos` style contains a `molecule` array (to skip
same-layer pairs) that was originally declared as `DAT::t_tagint_1d_randomread`
and was fixed to `typename AT::t_tagint_1d_randomread` in this session.

**Checklist:** when writing any new KK pair style, audit every per-atom
array member with `grep "DAT::" pair_<name>_kokkos.h`.  Only these are safe
as `DAT::`:

- `DAT::ttransform_kkacc_1d k_eatom`
- `DAT::ttransform_kkacc_1d_6 k_vatom`
- `DAT::ttransform_kkfloat_2d k_cutsq`
- `DAT::tdual_kkfloat_2d k_cut_ljsq` / `k_cut_coulsq`

Everything else that is set to `atomKK->k_<field>.view<DeviceType>()` must
use `typename AT::`.

### `coul/shield`: per-molecule filtering via the `molecule` view

`pair_coul_shield` skips interactions between atoms in the same layer (same
molecule id).  In the KOKKOS kernel this is implemented with a per-atom
`molecule` view:

```cpp
// in compute_fcoul / compute_ecoul:
if (molecule(i) == molecule(j)) return static_cast<KK_FLOAT>(0.0);
```

The view is synced in `compute()` before the kernel launch:

```cpp
molecule = atomKK->k_molecule.view<DeviceType>();
```

Because the `molecule` view's type must match `DeviceType`'s memory space,
it must be declared `typename AT::t_tagint_1d_randomread`, not
`DAT::t_tagint_1d_randomread`.  This is the concrete example that motivated
the general rule described above.

### `coul/shield`: inline taper function

The `coul/shield` potential uses a polynomial taper function to smooth the
force to zero at the cutoff.  In the KOKKOS kernel, all math must be
device-portable.  The taper and its derivative are computed via Horner's
method with hardcoded coefficients (exactly the same coefficients as the CPU
base class), called directly inside `compute_fcoul` and `compute_ecoul`.
No LAMMPS utility function is needed; the 8 polynomial coefficients are
embedded as literal constants in the kernel.

### `coul/cut/global`: thin wrapper inheriting from `PairCoulCutKokkos`

`pair_coul_cut_global_kokkos` does not need its own full implementation.
It is sufficient to:
1. Create a subclass `PairCoulCutGlobalKokkos<DeviceType>` that inherits
   from `PairCoulCutKokkos<DeviceType>`.
2. Override `coeff()` to enforce exactly two arguments (as in the base CPU
   class).
3. Override `extract()` to return the `epsilon` pointer.
4. Redirect the `PairStyle` macro entries to the new class.

No new `compute_fpair`/`compute_fcoul` is needed; the parent handles all
kernel logic.

### `lj/cut/sphere`: per-atom radius with a global `cutsq` ceiling

`pair_lj_cut_sphere_kokkos` reads per-atom radii and computes per-pair
cutoffs on-the-fly inside `compute_fpair`/`compute_evdwl`.  `cutsq` stores
the **maximum possible** squared cutoff per atom-type pair (set in
`init_one()`), so the neighbor list is built conservatively.  Inside the
kernel, the actual cutoff for atoms `i`–`j` is
`(sigma_ij + radius[i] + radius[j])^2`, and pairs beyond that distance
return zero force regardless of the neighbor list entry.

The `radius` per-atom view is declared `typename AT::t_kkfloat_1d_randomread
radius` (using `typename AT::`, not `DAT::`) and is synced in `compute()`:

```cpp
radius = atomKK->k_radius.view<DeviceType>();
```

The same `typename AT::` rule applies here as to `molecule`.

---

## Candidate Fix Styles to Port

This section covers porting of **fix styles** from `src/` and the `EXTRA-FIX` package.
Unlike pair styles, fix styles implement a wide range of lifecycle hooks
(`initial_integrate`, `final_integrate`, `post_force`, `end_of_step`, etc.) and
have very different complexity profiles.  The groups below are labeled **A through H**
and ordered from simplest to most complex (or least to most GPU-unfriendly).

### Key differences from pair-style porting

1. **No `pair_kokkos.h` template** — fix styles have no equivalent of the
   pairwise-dispatch framework.  Each fix implements one or more lifecycle methods
   directly with `Kokkos::parallel_for` / `Kokkos::parallel_reduce`.

2. **Fix lifecycle hooks** — identify which hooks need acceleration.  Common
   target hooks are `post_force`, `initial_integrate`, `final_integrate`, and
   `end_of_step`.  Hooks that run rarely (e.g., `write_restart`) stay on the
   CPU and do not need porting.

3. **`datamask_read` / `datamask_modify`** — every KOKKOS fix must declare
   which per-atom arrays it reads and which it modifies, using the bitmask
   constants from `atom.h` (e.g., `F_MASK`, `V_MASK`, `X_MASK`).
   Override `get_compute_flag()` if needed.

4. **`atomKK` sync/modified** — before a kernel reads a per-atom array, call
   `atomKK->sync<DeviceType>(flag)`.  After writing, call
   `atomKK->modified<DeviceType>(flag)`.  Use the same `typename AT::` rule
   for per-atom view members (see Key Rule 3 in the pair-styles section).

5. **KokkosBase mixin** — fixes that use per-atom stored data (those that
   implement `grow_arrays`, `copy_arrays`, `pack_exchange`, `unpack_exchange`)
   must also inherit from `KokkosBase` so the Kokkos atom-map machinery
   correctly manages the dual views.  See `fix_spring_self_kokkos` for the
   pattern.

6. **`Install.sh` entries** — same format as pair styles.  When the base fix
   lives in `src/`, omit the second argument; when it lives in an optional
   package (e.g., EXTRA-FIX), pass both filenames:

   ```sh
   # Base in src/:
   action fix_lineforce_kokkos.cpp
   action fix_lineforce_kokkos.h

   # Base in EXTRA-FIX:
   action fix_drag_kokkos.cpp fix_drag.cpp
   action fix_drag_kokkos.h   fix_drag.h
   ```

7. **Documentation** — update `doc/src/fix_<name>.rst` (add `/kk` index entry
   and `Accelerator Variants:` line) and `doc/src/Commands_fix.rst` (add letter
   `k` to the accelerator string).

---

## Group A — Trivial per-atom force/velocity loops (src/ and EXTRA-FIX) ✓ DONE

**Complexity:** Very low.  A single `Kokkos::parallel_for` over `[0, nlocal)`
reading `f[]` (or `v[]`) and writing back.  No per-atom stored state, no global
reduction, no communication.  Pattern: inherit from the base fix, override
`post_force()` (and `min_post_force()` if present) with a Kokkos functor.
Follow `fix_viscous_kokkos` as the primary template.

| Fix style | Package | Operation | Notes |
|---|---|---|---|
| `nve/noforce` | `src/` | Advance positions: `x += dt*v` | No force reads; pattern: `fix_nve_kokkos` |
| `aveforce` | `src/` | Replace per-atom `f` with group-average force | Needs a global `Kokkos::parallel_reduce` + MPI allreduce before applying; similar to `fix_addforce_kokkos` |
| `viscous/sphere` | `EXTRA-FIX` | Viscous drag + angular drag (torque) | Reads `omega`, writes `torque`; follow `fix_viscous_kokkos` |
| `drag` | `EXTRA-FIX` | Drag force toward target point/line/plane | Simple per-atom `f` update |
| `oneway` | `EXTRA-FIX` | Zero velocity component when moving the wrong way | Simple per-atom `v` check + zero |

---

## Group B — Wall styles extending FixWall (src/ and EXTRA-FIX) ✓ DONE

**Complexity:** Low.  All are subclasses of `FixWall` and need only
`precompute(int)` (precompute per-wall coefficients) and `wall_particle(int,int,double)`
(compute force/energy for one atom at given wall distance) to be ported.
Follow `fix_wall_lj93_kokkos` exactly.

| Fix style | Package | Potential form | Notes |
|---|---|---|---|
| `wall/harmonic` | `src/` | Harmonic spring | Simplest: `F = k*(r0-r)` |
| `wall/lj126` | `src/` | LJ 12-6 | Follow `wall/lj93/kk` directly |
| `wall/lj1043` | `src/` | LJ 10-4-3 integrated | Slightly more complex `wall_particle` |
| `wall/morse` | `src/` | Morse potential | Two exponential terms; follow `wall/lj93/kk` |

**Note:** `wall/reflect/stochastic` (EXTRA-FIX) inherits `FixWallReflect` (already ported
as `wall/reflect/kk`), adds a Marsaglia random-number generator for thermostatting.  The
random-number generator cannot be called from device code without extra work (see Group F).

---

## Group C — Sphere-body NVE/NVT/NPT integrators extending FixNH ✓ DONE

**Complexity:** Low-to-moderate.  `FixNHSphere` extends `FixNH` (already
ported as `fix_nh_kokkos`) by overriding `nve_v()`, `nve_x()`, and `nh_v_temp()`
to operate on both linear and angular velocity.  The KOKKOS variants
`fix_nph_kokkos` / `fix_npt_kokkos` / `fix_nvt_kokkos` already exist; the sphere
variants merely add the angular DOF loops.  Follow `fix_nve_sphere_kokkos` for the
angular loop pattern (or `fix_nve_asphere_kokkos` for the rigid-body variant).

| Fix style | Package | Base class chain | Notes |
|---|---|---|---|
| `nh/sphere` | `src/` | `FixNHSphere` → `FixNH` → `Fix` | Adds `nve_v_sphere`, `nve_x_sphere`, `nh_v_temp_sphere` |
| `nph/sphere` | `src/` | `FixNPHSphere` → `FixNHSphere` | Thin wrapper; most work in `nh/sphere/kk` |
| `npt/sphere` | `src/` | `FixNPTSphere` → `FixNHSphere` | Same as above |
| `nvt/sphere` | `src/` | `FixNVTSphere` → `FixNHSphere` | Same as above |

## compute temp/sphere ✓ DONE

**Complexity:** Low.  Mirrors `compute_temp_kokkos` but adds angular velocity
terms.  Sphere particles always use `rmass` (no per-type mass path), so there
is no `RMASS` template parameter — only a `MODE` parameter (ALL vs ROTATE-only).
The `dof_compute()` and `init()` methods are kept CPU-side (called infrequently
and already optimised).  Thin-wrapper constructors for the sphere NH fixes are
updated to create `temp/sphere/kk` instead of `temp/sphere`.

| Compute style | Package | Notes |
|---|---|---|
| `temp/sphere/kk` | `src/KOKKOS/` | `compute_temp_sphere_kokkos.{h,cpp}` |

---

## Group D — Per-atom stored-state fixes (grow_arrays + KokkosBase)

**Complexity:** Moderate.  These fixes maintain per-atom arrays that survive
atom migration (they implement `grow_arrays`, `copy_arrays`, `pack_exchange`,
`unpack_exchange`).  The KOKKOS version must:

1. Inherit from both the base fix and `KokkosBase`.
2. Replace CPU arrays with `Kokkos::DualView` members.
3. Override `grow_arrays` to resize dual views.
4. Override `pack_exchange_kokkos` / `unpack_exchange_kokkos`.

Use `fix_spring_self_kokkos` as the primary template.

| Fix style | Package | Per-atom data | Notes |
|---|---|---|---|
| `spring/rg` | `EXTRA-FIX` | Radius-of-gyration spring; needs global reduce | Moderate; follow `spring/self/kk` |
| `ti/spring` | `EXTRA-FIX` | Reference positions `xoriginal[nmax][3]` | Similar to `spring/self` |
| `addtorque/group` | `EXTRA-FIX` | Adds torque to group; needs COM reduce | Reads `x`, `torque`; pattern like `addforce/kk` |

---

## Group E — Thermostat/barostat coupling with global statistics

**Complexity:** Moderate.  These fixes compute a global temperature or pressure,
then rescale velocities/forces.  They call `Temperature->compute()` or
`Pressure->compute()` on the CPU, then launch a Kokkos kernel to rescale.
Pattern: `fix_temp_rescale_kokkos`, `fix_temp_berendsen_kokkos`.

| Fix style | Package | Coupling type | Notes |
|---|---|---|---|
| `press/berendsen` | `src/` | Isotropic/anisotropic box rescaling | Temperature + pressure compute → box + velocity rescale |
| `press/langevin` | `EXTRA-FIX` | Stochastic pressure coupling | Similar to `press/berendsen` with random term |
| `temp/csvr` | `EXTRA-FIX` | Canonical sampling; velocity rescaling | Stochastic factor; needs host-side random then device rescale |
| `temp/csld` | `EXTRA-FIX` | Canonical sampling; Lowe-Denbigh-Andersen | Per-pair velocity reassignment; moderate complexity |
| `gjf` | `EXTRA-FIX` | Grønbech-Jensen/Farago Langevin integrator | Has per-atom `vhalf[nmax][3]`; inherits Langevin structure (`fix_langevin_kokkos` pattern) |

---

## Group F — Moderate complexity: geometry/indenter/move/restraints

**Complexity:** Moderate-to-high.  These styles have variable run-time geometry,
tabulated data, or complex per-atom logic that requires careful structuring of
device kernels.

| Fix style | Package | Operation | Notes |
|---|---|---|---|
| `indent` | `src/` | Force from sphere/cylinder/plane indenter | Runtime shape dispatch; evaluate distance from shape |
| `move` | `src/` | Prescribe atom trajectories (linear/wiggle/rotate/variable) | Four motion modes; complex per-atom state; `initial_integrate` |
| `heat` | `src/` | Add kinetic energy to two groups | Needs two-group velocity scale; similar to `temp/rescale` |
| `restrain` | `src/` | Harmonic restraints on bonds/angles/dihedrals | Needs per-restraint neighbor lookup; moderate loop |
| `deform/pressure` | `EXTRA-FIX` | Deformation with pressure control | Extends `fix_deform_kokkos` (ported); adds pressure coupling |

---

## Documentation Checklist per Fix Style

For each fix style ported, update the following:

1. **`src/KOKKOS/fix_<name>_kokkos.h`** — new file; see structural notes below.
2. **`src/KOKKOS/fix_<name>_kokkos.cpp`** — new file.
3. **`src/KOKKOS/Install.sh`** — add `action` lines (omit second argument if base
   is in `src/`, include both filenames if base is in an optional package).
4. **Base class `.cpp`** — add `if (copymode) return;` as first line of destructor
   if not already present.
5. **`doc/src/fix_<name>.rst`**:
   - Add `.. index:: fix_style <name>/kk` near the top.
   - Add `*<name>/kk*` to `Accelerator Variants:`.
   - Order accelerator variants alphabetically.
6. **`doc/src/Commands_fix.rst`**:
   - Add letter `k` to the accelerator string for this style.
   - Use alphabetical order within the letters string (e.g., `ko` not `ok` or `o,k`).

### Typical fix KOKKOS header skeleton

```cpp
#ifdef FIX_CLASS
// clang-format off
FixStyle(<name>/kk,Fix<Name>Kokkos<LMPDeviceType>);
FixStyle(<name>/kk/device,Fix<Name>Kokkos<LMPDeviceType>);
FixStyle(<name>/kk/host,Fix<Name>Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_<NAME>_KOKKOS_H
#define LMP_FIX_<NAME>_KOKKOS_H

#include "fix_<name>.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFix<Name>{};

template<class DeviceType>
class Fix<Name>Kokkos : public Fix<Name> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  Fix<Name>Kokkos(class LAMMPS *, int, char **);
  ~Fix<Name>Kokkos() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFix<Name>, const int &) const;

 private:
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif
```

### `datamask_read` / `datamask_modify` convention

Declare in the constructor body which per-atom arrays the fix touches:

```cpp
Fix<Name>Kokkos<DeviceType>::Fix<Name>Kokkos(LAMMPS *lmp, int narg, char **arg)
    : Fix<Name>(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read   = F_MASK | MASK_MASK;
  datamask_modify = F_MASK;
}
```

### Syncing per-atom arrays before and after a kernel

```cpp
void Fix<Name>Kokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync<DeviceType>(datamask_read);

  x    = atomKK->k_x.view<DeviceType>();
  f    = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFix<Name>>(0, nlocal), *this);

  atomKK->modified<DeviceType>(datamask_modify);
}
```

---

## Lessons Learned from Fix Porting

### Fix styles vs pair styles: key structural difference

Pair styles always implement `compute(int eflag, int vflag)` and fit the
`pair_kokkos.h` dispatch template.  Fix styles implement whichever lifecycle
hooks are relevant.  Before porting, identify the set of hooks used:

```bash
grep "void Fix[A-Z]" src/<PKG>/fix_<name>.cpp | grep -v "^//"
```

Only those hooks need Kokkos kernels; others can call the base-class implementation.

### Wall styles (FixWall subclasses): use the `wall_particle` override pattern

All `FixWall` subclasses share a generic `post_force` loop in `FixWall::post_force`.
The KOKKOS port of `fix_wall_lj93_kokkos` replaces that loop with a Kokkos functor
and calls `wall_particle` inline.  Every other `FixWall` subclass can follow the
same override, needing only a new `precompute` (compute per-wall coefficients) and
a new `wall_particle` (compute single-atom force from wall distance):

```cpp
template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWall<Name>Kokkos<DeviceType>::wall_particle(int m, int which, double coord) const
{
  // compute fwall and ewall for atom at position coord from wall m
}
```

### Wall styles: `private:` must become `protected:` for coeff arrays

The base-class coeff arrays (`coeff1[]`, `coeff2[]`, ..., `offset[]`) are often
declared `private:`.  Because the KK subclass inherits from the base and reads
those arrays in `precompute()`, they must be changed to `protected:`.  Always
check with:

```bash
grep -n "private:\|coeff\|offset" src/fix_wall_<name>.h
```

### Wall styles: `fix_wall_harmonic` has no per-wall DualViews

`wall/harmonic` does not precompute coefficients (it computes `dr = cutoff - delta`
inline).  The `precompute()` override is empty in the KK class.  The base-class
`epsilon[]` and `cutoff[]` (protected in `FixWall`) are accessed directly as scalar
casts in the functor since they are small fixed arrays (6 elements) that remain on
the host.  For consistency with other wall KK styles the DualView pattern is still
used for `k_vatom` but no coefficient DualViews are needed.

### Wall styles: `fix_wall_lj1043` — avoid `powint` in device kernels

`powint(x, -N)` is a **host-only** function from `math_special.h`.  In the base-class
`wall_particle`, integer powers like `powint(delta + coeff4[m], -3)` must be replaced
in the device functor with explicit repeated multiplication:

```cpp
KK_FLOAT dc4inv = 1.0 / (delta + d_coeff4(m));
KK_FLOAT dc4inv3 = dc4inv * dc4inv * dc4inv;      // replaces powint(..., -3)
KK_FLOAT dc4inv4 = dc4inv3 * dc4inv;              // replaces powint(..., -4)
```

Store `coeff4[]` in its own DualView so the device kernel can read it.

### Wall styles: `fix_wall_morse` — `exp` is available on device

Morse potential uses `exp()`, which is available in CUDA/HIP device code as
`Kokkos::exp()`.  Simply replace `exp(...)` with `Kokkos::exp(...)` in the
device functor.  The parameters `alpha[]`, `sigma[]`, and `epsilon[]` are defined
in the `FixWall` base (`protected:`), so copy them to DualViews in `precompute`.

### Wall styles: `post_force` wraps base and handles virial

The KK `post_force` override only:
1. Reallocates `k_vatom` if `vflag_atom` is set.
2. Calls `FixWall<Name>::post_force(vflag)` — which calls the overridden
   `wall_particle` (the KK version).
3. Syncs `k_vatom` back to host if needed.

The virial accumulation inside `wall_particle` uses the `result[7..12]` array
and is flushed to `virial[0..5]` after `parallel_reduce` returns.

### Sphere integrators: extending FixNHKokkos

`FixNHKokkos` already exposes virtual `nve_v`, `nve_x`, `nh_v_temp` hooks.
The sphere variants only add angular velocity updates.  The Kokkos port of
`FixNHSphereKokkos` should add a second functor that updates `omega` alongside
`v`, dispatched from the same `initial_integrate` / `final_integrate` call.

### Sphere NH integrators: inherit from FixNHKokkos, not FixNHSphere

The instinct is to inherit `FixNHSphereKokkos` from `FixNHSphere` (like
`FixWallLJ93Kokkos` inherits from `FixWallLJ93`).  But for NH sphere, this is
wrong: `FixNHSphere::nve_v()` calls `FixNH::nve_v()` internally — the CPU
version — losing the Kokkos parallelism.

**Correct approach:** inherit `FixNHSphereKokkos<DeviceType>` from
`FixNHKokkos<DeviceType>`.  This brings in the Kokkos implementations of
`nve_v`, `nve_x`, `nh_v_temp`, and `nh_v_press`.  Override just the three
sphere-specific methods:

```
nve_v()      → call FixNHKokkos::nve_v(),      then launch omega update kernel
nve_x()      → call FixNHKokkos::nve_x(),      then launch dipole update kernel (if needed)
nh_v_temp()  → call FixNHKokkos::nh_v_temp(),  then launch omega scaling kernel
```

The constructor logic from `FixNHSphere` (omega/radius flag checks, `inertia`,
`disc` keyword) is duplicated directly in `FixNHSphereKokkos`.

### NH sphere: FixNH::omega[6] vs per-atom omega view

`FixNH` (inherited via `FixNHKokkos`) has a `protected: double omega[6]`
member for the barostat degrees of freedom.  Adding a per-atom Kokkos view
also named `omega` in `FixNHSphereKokkos` would shadow this.  Use a distinct
name — `omega_kk` — for the per-atom Kokkos view to avoid confusion.

### NH sphere: base class views remain valid after virtual-call chain

When `FixNHSphereKokkos::nve_v()` calls `FixNHKokkos<DeviceType>::nve_v()`
first, the base sets `this->mask` and `this->rmass` to the current device
views before its `parallel_for`.  After it returns, these views are still
valid (they are reference-counted handles).  The sphere omega kernel can
therefore reuse `this->mask` and `this->rmass` directly without resyncing,
and only needs to additionally sync `OMEGA_MASK`, `TORQUE_MASK`, and
`RADIUS_MASK`.

### NH sphere: thin-wrapper constructors use temp/sphere/kk, not temp/sphere

The non-KK thin wrappers (`nph/sphere`, `npt/sphere`, `nvt/sphere`) create
`temp/sphere` computes.  The KK thin wrappers (`nph/sphere/kk` etc.) should
instead create `temp/sphere/kk` computes (once `compute_temp_sphere_kokkos`
exists) so that the temperature calculation also runs on the device.
Using `temp/sphere` would still give the correct temperature but would require
a host synchronisation on every thermostat step, defeating the purpose of the
KOKKOS port.  Using `temp/kk` would give the *wrong* temperature (ignores
angular velocity).

### NH sphere: DLM dipole integrator not supported in KK mode

The DLM (Dullweber-Leimkuhler-Maclachlan) dipole integrator in
`FixNHSphere::nve_x()` requires multiple matrix-vector products
(`BuildRxMatrix`, `matvec`, `transpose_times3`, etc.).  These are feasible on
device in principle but complex and rarely used.  The KK port supports only
the simple dipole integrator (`dlm_flag == 0`); requesting `update dipole dlm`
raises an error.


### compute temp/sphere: sphere particles always use rmass, so no RMASS template

`compute temp` has an `RMASS` template parameter to handle both
per-atom mass (`rmass`) and per-type mass (`mass[type[i]]`).
`compute temp/sphere` does NOT need this because `atom_style sphere` always
provides `rmass` — the per-type mass path is never taken.  The KK variant
therefore uses a simpler `MODE` template parameter (0 = ROTATE-only,
1 = ALL, i.e. translational + rotational), with `rmass` always used.

### compute temp/sphere: view name collision with FixNH::omega

`FixNH` has a `protected: double omega[6]` member for barostat DOF.
`ComputeTempSphereKokkos` is a `Compute`, not a fix, so there is no
collision here — but for consistency and clarity, the per-atom Kokkos
omega view is named `omega_kk` in the compute too.

### Per-atom storage: dual-view pattern from fix_spring_self_kokkos

When a fix stores per-atom data that migrates with atoms (`grow_arrays` /
`pack_exchange` / `unpack_exchange`), replace the plain C arrays with
`Kokkos::DualView` members and inherit from `KokkosBase`.  The
`pack_exchange_kokkos` override packs atoms directly from the device view
without falling back to host memory, preserving GPU locality.

## Trust These Instructions

These instructions are tested and validated. Only search for additional information if:
- A specific command fails with an error
- You need details about a specific package's requirements
- Instructions appear outdated based on error messages
- Working with advanced features not covered here (GPU, Kokkos backends, etc.)

For package-specific documentation, build options, and advanced features, refer to https://docs.lammps.org

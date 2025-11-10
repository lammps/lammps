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

5. **Documentation:** All new commands or features must be documented. Put `.. versionadded:: TBD` or
   `.. versionchanged:: TBD` in front of paragraphs documenting the new or changed functionality.
   The `TBD` will be manually replaced with the release version string during the release preparation.

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

When performing a code review, check any changes to the documentation (in the
`doc/src/` folder) to be written in American English and with plain ASCII characters.

When performing a code review, ensure that the documentation for any new
commands or added keywords to existing commands contains a `.. versionadded:: TBD`
directive.  For any modified commands or keywords a `.. versionchanged:: TBD`
directive should be included in the documentation. Check if any examples use
the new or modified commands and check if they need updating. Ensure that
building the documentation with "make html", "make pdf", and "make spelling"
can complete and does *NOT* produce any *NEW* warnings or errors.

## Trust These Instructions

These instructions are tested and validated. Only search for additional information if:
- A specific command fails with an error
- You need details about a specific package's requirements
- Instructions appear outdated based on error messages
- Working with advanced features not covered here (GPU, Kokkos backends, etc.)

For package-specific documentation, build options, and advanced features, refer to https://docs.lammps.org

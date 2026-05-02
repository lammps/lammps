# GEMINI.md - LAMMPS Development Guide

This document provides specialized instructions and context for the Google Gemini CLI when working on the LAMMPS codebase. LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is a mature, high-performance molecular dynamics code.

## Core Mandates & Engineering Standards

### Security & System Integrity
- **7-bit ASCII only:** All source files, documentation, and scripts MUST use 7-bit US-English ASCII characters only. Unicode is strictly forbidden.
- **No VLAs:** Variable-length arrays are not allowed. Use `memory->create()` or `std::vector` instead.

### Build System
- **CMake is Primary:** ALWAYS use CMake for new builds. Traditional Make is legacy.
- **CMake Source Directory:** Use `-S cmake`, NOT `-S .`. The `CMakeLists.txt` is in the `cmake/` directory.
- **Out-of-source builds:** ALWAYS create a separate `build/` directory. Never build in the root or `src/`.
- **Switching builds:** Run `make -C src purge` before switching from Make to CMake.

### Coding Standards
- **Standard:** C++17.
- **Logical Operators:** Use `&&`, `||`, `!`, NOT `and`, `or`, `not`.
- **Formatting:** Adhere to `.clang-format` in `src/`. Use `fmtlib` for string formatting.
- **Inheritance:** Most features are implemented as "styles" inheriting from base classes (`Pair`, `Fix`, `Compute`, `Command`, etc.).
- **Error Handling:** Use `error->all()` or `error->one()` for reporting errors.

## Repository Structure

- `src/`: Core source code and optional packages in subdirectories.
- `cmake/`: CMake build system configuration.
- `doc/src/`: Documentation in reStructuredText (`.rst`).
- `lib/`: Bundled external libraries (e.g., Kokkos, Colvars).
- `unittest/`: CTest-based unit tests.
- `tools/`: Development utilities, including `coding_standard` and `regression-tests`.
- `examples/`: Input script examples for testing.

## Development Workflow

### 1. Research & Strategy
- When adding a new style, identify a similar existing style in `src/` or a package as a template.
- Check `https://docs.lammps.org/latest/` for architectural requirements: [Modify Style](https://docs.lammps.org/Modify_style.html).

### 2. Implementation
- Place new style files (`style_name.cpp`, `style_name.h`) in the appropriate `src/` or package directory.
- Add new files added to packages to `src/.gitignore`.
- If renaming/removing files, update `src/Purge.list`.

### 3. Documentation
- Create/update the corresponding `.rst` file in `doc/src/`.
- Use `.. versionadded:: TBD` for new features or `.. versionchanged:: TBD` for modified behavior.
- Ensure all documentation is ASCII and American English.

### 4. Validation
- **Style Checks:** Run `make check` in `src/` to verify whitespace and formatting.
- **Build:** Build with a preset like `cmake -C cmake/presets/most.cmake`.
- **Unit Tests:** Run `ctest` in the build directory (requires `-D ENABLE_TESTING=on`).
- **Regression Tests:** Are not necessary during development. They will be run automatically by the available CI infrastructure after code review.

## Common Commands

```bash
# Style and formatting checks
cd src && make check

# Recommended build sequence
mkdir -p build-gemini
cmake -S cmake -B build-gemini -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake -D ENABLE_TESTING=on -G Ninja
cmake --build build-gemini

# Run unit tests
cd build && ctest -V

# Fix whitespace issues
cd src && make fix-whitespace
```

## Documentation Validation
```bash
cd doc
make html          # Build HTML and check for syntax errors
make spelling      # Check spelling (exceptions in doc/utils/sphinx-config/false-positives.txt)
make anchor_check  # Check for duplicate anchors
```

Refer to `.github/copilot-instructions.md` for additional technical details and common pitfalls.

---
applyTo: "cmake/**,src/Makefile"
---

# LAMMPS Build System (CMake modules and src/Makefile)

The CMake quick start is in `.github/copilot-instructions.md`.  This file covers how
CMake collects the sources and conventions for changing the build machinery.

## CMake is the only build system

The legacy GNU make build was removed after the 30Sep2026 stable release
(`stable_22Jul2025` is the last version with the complete make build).  `src/Makefile`
only provides the `check-*`/`fix-*` style check targets, `clean`, and `purge`; all other
former targets print a message pointing to CMake.  `make -C src purge` removes files
that the make build had copied or generated in `src/` (package sources with the same
name as a file in a package folder, `style_*`/`packages_*.h` files, `lib/gpu/*_cubin.h`);
CMake refuses to configure while such files are present (`check_for_autogen_files()` in
`LAMMPSUtils.cmake`, `DetectBuildSystemConflict()` in `StyleHeaderUtils.cmake`).  Never
add new files or instructions for the make build (`Install.sh`, `Makefile.lammps`,
`src/.gitignore` or `Purge.list` entries).

## How CMake collects sources

- `cmake/CMakeLists.txt` globs all `*.cpp`/`*.h` files in `src/` and in the folder of
  every enabled package and detects styles from the `XxxStyle()` macros in the headers.
  New files in a package folder are picked up automatically after re-running the CMake
  configuration (`cmake -S cmake -B <build>`); no build file needs to be edited.
- Suffix packages (GPU, INTEL, KOKKOS, OPENMP, OPT, and CORESHELL) use
  `RegisterStylesExt()`: `<base>_<suffix>.{h,cpp}` is only compiled when `<base>.h` is a
  style header of the build (core or an enabled package), so variants of styles from
  optional packages follow that package automatically.  Files of a suffix package
  without a style macro (base classes, helpers) must be listed explicitly in
  `cmake/Modules/Packages/<PKG>.cmake`, inside an `if(PKG_<NAME>)` block if they depend
  on an optional package.
- A file in an ordinary package that depends on ANOTHER optional package is removed
  with `list(REMOVE_ITEM ...)` of the style header and source list in
  `cmake/Modules/Packages/<PKG>.cmake` when that package is disabled (see `MISC.cmake`,
  `MC.cmake`).

## CMake module conventions

- Declare the URL and checksum of downloaded libraries with
  `SetDownloadSettings(<prefix> <name> <url> <sha256>)` from
  `cmake/Modules/LAMMPSUtils.cmake` (keep the copy in
  `cmake/Modules/LAMMPSInterfacePlugin.cmake` in sync), not with hand-written cached
  `<PREFIX>_URL`/`<PREFIX>_SHA256` variables.
- An autotools-based `ExternalProject_Add()` needs `AutotoolsTouch.cmake` as its
  `PATCH_COMMAND` (see the MBX, PLUMED, and SCAFACOS modules): with policy CMP0135 NEW,
  extracted files get extraction-time stamps and automake rebuild rules fire
  (`missing aclocal-1.16`).
- Build workarounds (source patches for downloaded libraries, timestamp fixes) must be
  transparent to users: explain them in CMake comments and commit messages, not in the
  user manual.

---
applyTo: "cmake/**,src/Makefile*,src/MAKE/**,src/*.sh,src/*/Install.sh,lib/*/Makefile*,lib/*/Install.py"
---

# LAMMPS Build System (CMake modules and the legacy make build)

The CMake quick start is in `.github/copilot-instructions.md`.  This file covers the
legacy GNU make build and conventions for changing the build machinery.

## Legacy GNU make build

`cd src && make serial` (or `make mpi`) builds `lmp_serial`/`lmp_mpi`.  Enable or disable
packages first with `make yes-<package>`/`make no-<package>` or bundles like
`make yes-basic` (MANYBODY, MOLECULE, KSPACE, RIGID); `make pi` shows the package status.
Packages that need external libraries or downloads are CMake-only.  The make build
copies package files into `src/`, which is why package files need `src/.gitignore` and
`src/Purge.list` entries.

**Switching build systems:** make -> CMake requires `make -C src purge` first; CMake ->
make requires `make -C src clean-all` first.  CMake errors out if it detects
make-generated header files in `src/`.

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

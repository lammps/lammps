# Testing the plugins of packages with external libraries

The plugins in `examples/PACKAGES/*/plugin` (ML-PACE, PLUMED, MBX) are separate CMake projects
that users build after LAMMPS, and the Windows installers bundle them.  This folder has what the
continuous integration uses to build them quickly and to check that they load:

- `mock-libs.cmake`: initial cache file (`cmake -C`) for configuring a plugin.  It packs the
  stand-ins `mock-plumed` and `mock-mbx` into archives in the build folder and sets the download
  URLs and checksums of PLUMED and MBX to them, so that the plugins build them instead of the real
  libraries (CMake warns that the download settings differ from the defaults).
- `mock-plumed`, `mock-mbx`: stand-ins for the source archives of PLUMED and MBX.  Each has a
  `configure` script and make rules that take the options and produce the files (libraries,
  headers, programs, folders) that the LAMMPS CMake scripts expect, for Linux and for the
  cross-compilation for Windows.  The headers declare only what the LAMMPS sources use, and the
  functions do nothing: the plugins compile, link and load, but do not work.  When a package starts
  to use another function of its library, it must be added to the stand-in.  ML-PACE needs no
  stand-in: its library compiles in a few seconds.
- `check-plugins.sh`: loads plugins into LAMMPS and compares the styles that LAMMPS lists with
  those that the plugin sources register.  LAMMPS only prints a warning when a plugin cannot be
  loaded, so this comparison is the test.

To run the same check on Linux:

```bash
cmake -S cmake -B build -C cmake/presets/basic.cmake -D PKG_PLUGIN=on \
      -D BUILD_SHARED_LIBS=on -D BUILD_MPI=off
cmake --build build
for p in pace plumed mbx; do
  cmake -S examples/PACKAGES/$p/plugin -B build/${p}plugin -C unittest/plugins/mock-libs.cmake \
        -D LAMMPS_SOURCE_DIR=$PWD/src -D BUILD_MPI=off
  cmake --build build/${p}plugin
done
sh unittest/plugins/check-plugins.sh build/lmp \
   examples/PACKAGES/pace/plugin/paceplugin.cpp build/paceplugin/paceplugin.so \
   examples/PACKAGES/plumed/plugin/plumedplugin.cpp build/plumedplugin/plumedplugin.so \
   examples/PACKAGES/mbx/plugin/mbxplugin.cpp build/mbxplugin/mbxplugin.so
```

The plugins must be configured with the same MPI setting as LAMMPS.  For the Windows installers,
`cmake/packaging/cmake-win-on-linux.py -c unittest/plugins/mock-libs.cmake` uses the stand-ins.

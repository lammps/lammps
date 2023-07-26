LAMMPS GUI

# Overview

This is essentially a small graphical text editor that is linked to the LAMMPS
library and thus can run LAMMPS using the contents of the text buffer as input.
The output of the run is captured and displayed in a separate window dialog.

# Prerequisites and portability

The GUI is programmed with the Qt cross-platform GUI toolkit, currently using
version 5.15LTS. It has been compiled and tested on:
- Fedora Linux 38 x86\_64
- Apple macOS 12 (Monterey) x86\_64
- Windows 10 x86_64 (via MinGW Linux to Windows cross-compiler)

To compile you need to have LAMMPS compiled (so it can be linked to the LAMMPS
C-library interface) and Qt 5.15LTS installed.

# Compilation

The source is bundled with the LAMMPS source code distribution and thus can be
built as part of the regular LAMMPS build by setting `-D BUILD_LAMMPS_GUI=on`.
All other settings are then inherited from the regular LAMMPS build.

It is also possible to build the GUI as a standalone compilation, then the
CMake configuration needs to be told where to find the LAMMPS headers and
the LAMMPS library, via `-D LAMMPS_SOURCE_DIR=/path/to/lammps/src`. CMake
will try to guess a build folder with the LAMMPS library from that path,
but it can also be set with `-D LAMMPS_LIB_DIR=/path/to/lammps/lib`.

It is also possible to compile the GUI without embedding the LAMMPS library
directly, but to dynamically load it later as a shared library. This has
the advantage that you can update LAMMPS without having to recompile the GUI.
This is triggered by setting `-D LAMMPS_GUI_USE_PLUGIN=on` and then
`-D LAMMPS_PLUGINLIB_DIR=/path/to/lammps/plugin/loader`. Typically, this
would be the `examples/COUPLE/plugin` folder of the LAMMPS distribution.


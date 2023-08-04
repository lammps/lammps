LAMMPS GUI

# Overview

LAMMPS GUI is essentially a small graphical text editor that is linked
to the LAMMPS library and thus can run LAMMPS using the contents of the
text buffer as input. This is similar to what people usually would do
using a text editor to edit the input and then a command line terminal
window to run the input commands.  The main benefit is that this
integrates very well with graphical desktop environments and that it is
easier to use for beginners in running computations and thus very
suitable for tutorials on LAMMPS.

# Features

The main window of the LAMMPS GUI is a text editor window with syntax
highlighting set up for LAMMPS input files.  It can be used to edit any
kind of text file, but then trying to run those files will cause errors.
The output of a run is captured and displayed in a separate window
dialog.  The log window is updated during the run and a progress bar for
each run command shown in the main window.  Starting a new run will open
another log windows.  After the simulation is finished, an image of the
simulated system can be created and shown in another window.  Ongoing
runs can be stopped at the next iteration via triggering a timeout.

When opening a file, the editor will determine the directory where the
file resides and switch its current working directory to that folder.
Many LAMMPS inputs contain commands that read other files, typically
from the folder with the input file.  The GUI will show the current
working directory.  The editor window can also receive (entire) files
via drag-n-drop from a file manager GUI or a desktop environment.

Almost all commands are accessible via hotkeys. Which those hotkeys are,
is shown next to the entries in the menu.  Log and image viewer windows
can be closed with CTRL-W (or Command-W on macOS).  The "About LAMMPS"
dialog will show the LAMMPS version and the features included into the
LAMMPS library linked to the LAMMPS GUI.

Due to its nature as a graphical application, it is not possible to use
the LAMMPS GUI in parallel with MPI, but OpenMP multi-threading is
available.

# Prerequisites and portability

LAMMPS GUI is programmed using the Qt cross-platform GUI toolkit,
currently using Qt version 5.15LTS for better compatibility with older
compilers. It has been successfully compiled and tested on:

- Fedora Linux 38 x86\_64 using GCC 13 and Clang 16
- Apple macOS 12 (Monterey) and macOS 13 (Ventura) with Xcode on arm64 and x86\_64
- Windows 10 and 11 x86_64 with Visual Studio 2022 and Visual C++ 14.36

# Compilation

The source for the LAMMPS GUI is included with the LAMMPS source code
distribution in the folder `tools/lammps-gui` and thus it can be can be
built as part of a regular LAMMPS compilation.  Using CMake is required.
To enable its compilation the CMake variable `-D BUILD_LAMMPS_GUI=on`
must be set when creating the CMake configuration.  All other settings
(compiler, flags, compile type) for LAMMPS GUI are then inherited from
the regular LAMMPS build.  If the Qt library is packaged for Linux
distributions, then its location is typically auto-detected since the
required CMake configuration files are stored in a location where CMake
can find them without additional help.  Otherwise, the location of the
Qt library installation must be indicated by setting `-D
Qt5_DIR=/path/to/qt5/lib/cmake/Qt5`, which is a path to a folder inside
the Qt installation that contains the file `Qt5Config.cmake`.

It is also possible to build the LAMMPS GUI as a standalone executable
(e.g. when LAMMPS has been compiled with traditional make), then the
CMake configuration needs to be told where to find the LAMMPS headers
and the LAMMPS library, via `-D LAMMPS_SOURCE_DIR=/path/to/lammps/src`.
CMake will try to guess a build folder with the LAMMPS library from that
path, but it can also be set with `-D LAMMPS_LIB_DIR=/path/to/lammps/lib`.

Rather than linking to the LAMMPS library during compilation, it is also
possible to compile the GUI with a plugin loader library that will load
the LAMMPS library at runtime during startup of the GUI from a shared
library; e.g. `liblammps.so` or `liblammps.dylib` or `liblammps.dll`
depending on the operating system.  This has the advantage that the
LAMMPS library can be updated LAMMPS without having to recompile the
GUI.  The ABI of the LAMMPS C-library interface is very stable and
generally backward compatible.  This feature is enabled by setting
`-D LAMMPS_GUI_USE_PLUGIN=on` and then
`-D LAMMPS_PLUGINLIB_DIR=/path/to/lammps/plugin/loader`. Typically, this
would be the `examples/COUPLE/plugin` folder of the LAMMPS distribution.

# Platform notes

## macOS

When building on macOS, the build procedure will try to manufacture a
drag-n-drop installer, LAMMPS-macOS-multiarch.dmg.  To build multi-arch
executables that will run on both, arm64 and x86_64 architectures
natively, it is necessary to set the CMake variable
`-D CMAKE_OSX_ARCHITECTURES=arm64;x86_64`.  To achieve wide compatibility
with different macOS versions, you can also set
`-D CMAKE_OSX_DEPLOYMENT_TARGET=11.0` which will set compatibility to macOS
11 (Big Sur) and later, even if you are compiling on a more recent macOS
version.

--------

updated by Axel Kohlmeyer, 08/1023

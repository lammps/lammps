LAMMPS GUI

# Overview

LAMMPS GUI is essentially a small graphical text editor that is linked
to the LAMMPS library and thus can run LAMMPS using the contents of the
text buffer as input.

This is similar to what people usually would do using a text editor to
edit the input and then open a command line terminal window to run the
necessary commands.  The main benefit of a GUI is that this integrates
very well with graphical desktop environments and many basic tasks can
be done directly from within the GUI without switching to a text
console.  This makes it easier for beginners to get started running
computations and thus very suitable for tutorials on LAMMPS.

# Features

The main window of the LAMMPS GUI is a text editor window with line
numbers syntax and highlighting set up for LAMMPS input files.  When
starting a run the output to the console is captured and displayed in a
log window.  Also, generated thermodynamic data is collected and shown
in chart window.  An ongoing run can be stopped at the next iteration
and after a run is completed, a snapshot image can be generated and
displayed in an image viewer window.

When opening a file, the editor will switch to the directory where the
file was found.  The current working directory is shown in the status
bar at the bottom of the window. The editor window can also receive
(entire) files via drag-n-drop from a file manager GUI or a desktop
environment.

Almost all commands are accessible via hotkeys. Which those hotkeys are,
is shown next to the entries in the menu.  Log and image viewer windows
can be closed with CTRL-W (or Command-W on macOS).
A number of settings can be adjusted via a Preferences dialog.

The "About LAMMPS" dialog will show the LAMMPS version and the features
included into the LAMMPS library linked to the LAMMPS GUI.

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
- Windows 10 and 11 x86_64 with a MinGW GCC cross-compiler on Linux

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
generally backward compatible.  This feature is enabled by setting `-D
LAMMPS_GUI_USE_PLUGIN=on` and then `-D
LAMMPS_PLUGINLIB_DIR=/path/to/lammps/plugin/loader`. Typically, this
would be the `examples/COUPLE/plugin` folder of the LAMMPS distribution.

# Platform notes

## macOS

When building on macOS, the build procedure will try to manufacture a
drag-n-drop installer, LAMMPS-macOS-multiarch.dmg, when using the 'dmg'
target (i.e. `cmake --build <build dir> --target dmg` or `make dmg`.

To build multi-arch executables that will run on both, arm64 and x86_64
architectures natively, it is necessary to set the CMake variable `-D
CMAKE_OSX_ARCHITECTURES=arm64;x86_64`.  To achieve wide compatibility
with different macOS versions, you can also set `-D
CMAKE_OSX_DEPLOYMENT_TARGET=11.0` which will set compatibility to macOS
11 (Big Sur) and later, even if you are compiling on a more recent macOS
version.


## Windows

On Windows either native compilation from within Visual Studio 2022 with
Visual C++ is supported and tested, or compilation with the MinGW / GCC
cross-compiler environment on Fedora Linux.

### Visual Studio

Using CMake and Ninja as build system are required. Qt needs to be
installed, tested was a binary package downloaded from
https://www.qt.io, which installs into the `C:\\Qt` folder by default.
There is a custom `x64-GUI-MSVC` build configuration provided in the
`CMakeSettings.json` file that Visual Studio uses to store different
compilation settings for project.  Choosing this configuration will
activate building the `lammps-gui.exe` executable in addition to LAMMPS
through importing package selection from the `windows.cmake` preset
file and enabling building the LAMMPS GUI and disable building with MPI.
When requesting an installation from the `Build` menu in Visual Studio,
it will create a compressed `LAMMPS-Win10-amd64.zip` zip file with the
executables and required dependent .dll files.  This zip file can be
uncompressed and `lammps-gui.exe` run directly from there.  The
uncompressed folder can be added to the `PATH` environment and LAMMPS
and LAMMPS GUI can be launched from anywhere from the command line.

### MinGW64 Cross-compiler

The standard CMake build procedure can be applied and the
`mingw-cross.cmake` preset used. By using `mingw64-cmake` the CMake
command will automatically include a suitable CMake toolset file (the
regular cmake command can be used after that).  After building the
libraries and executables, you can build the target 'zip'
(i.e. `cmake --build <build dir> --target zip` or `make zip`
to stage all installable files into a LAMMPS_GUI folder and then
run a script to copy all required dependencies, some other files,
and create a zip file from it.

Linux
"""""

Version 5.12 or later of the Qt library and CMake version 3.16 are
required and those are provided by, e.g., Ubuntu 20.04LTS.  Thus older
Linux distributions are not likely to be supported, while more recent
ones will work, even for pre-compiled executables (see above).  After
compiling with `cmake --build <build folder>`, use
`cmake --build <build folder> --target tgz` or `make tgz` to build
a `LAMMPS-Linux-amd64.tar.gz` file with the executables and their
support libraries.

--------

updated by Axel Kohlmeyer, 2023-08-05

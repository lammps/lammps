Using CMake with LAMMPS tutorial
================================

The support for building LAMMPS with CMake is a recent addition to
LAMMPS thanks to the efforts of Christoph Junghans (LANL) and Richard
Berger (Temple U).  One of the key strengths of CMake is that it is not
tied to a specific platform or build system and thus generate the files
necessary to build and develop for different build systems and on
different platforms.  Note, that this applies to the build system itself
not the LAMMPS code. In other words, without additional porting effort,
it is not possible - for example - to compile LAMMPS with Visual C++ on
Windows.  The build system output can also include support files
necessary to programm LAMMPS as a project in integrated development
environments (IDE) like Eclipse, Visual Studio, QtCreator, Xcode,
CodeBlocks, Kate and others.

A second important feature of CMake is, that it can detect and validate
available libraries, optimal settings, available support tools and so
on, so that by default LAMMPS will take advantage of available tools
without requiring to provide the details about how to enable/integrate
them.

The downside of this approach is, that there is some complexity
associated with running CMake itself and how to achieve desired
customizations and modifications to the LAMMPS configuration and
compilation.  And for as long as this facility is relatively new and
not as widely used as the traditional build process, there are chances
that the scripts that CMake processes may have bugs or are missing
options, despite the best efforts to test and verify its functionality.

This tutorial will show how to manage this through some selected
examples.  Please see the chapter about :doc:`building LAMMPS <Build>`
for descriptions of specific flags and options for LAMMPS in general and
for specific packages.

CMake can be used through either the command-line interface (CLI)
program ``cmake`` (or ``cmake3``), a text mode interactive user
interface (TUI) program ``ccmake`` (or ``ccmake3``), or a graphical user
interface (GUI) program ``cmake-gui``.  All of them are portable
software available on all supported platforms and can be used
interchangeably.  The minimum supported CMake version is 3.10 (3.12 or
later is recommended).


Prerequisites
-------------

This tutorial assumes that you are operating in a command-line environment
using a shell like Bash.

- Linux: any Terminal window will work
- MacOS X: launch the Terminal application.
- Windows 10: install and run the :doc:`Windows subsystem for Linux <Howto_bash>`

We also assume that you have downloaded and unpacked a recent LAMMPS source code package
or used Git to create a clone of the LAMMPS sources on your compilation machine.

You should change into the top level directory of the LAMMPS source tree all
paths mentioned in the tutorial are relative to that.  Immediately after downloading
it should look like this:

.. code-block:: bash

    $ ls
    bench  doc       lib      potentials  README  tools
    cmake  examples  LICENSE  python      src

Build directory versus source directory
---------------------------------------

When using CMake the build procedure is separated into multiple distinct phases:

  #. **Configuration:** detect or define which features and settings
     should be enable and used and how LAMMPS should be compiled
  #. **Compilation:** generate and compile all necessary source files
     and build libraries and executables.
  #. **Installation:** copy selected files from the compilation into
     your file system, so they can be used without having to keep the
     source and build tree around.

The configuration and compilation of LAMMPS has to happen in a dedicated
*build directory* which must be different from the source directory.
Also the source directory (``src``) must remain pristine, so it is not
allowed to "install" packages using the traditional make process and
after an compilation attempt all created source files must be removed.
This can be achieved with ``make no-all purge``.

You can pick **any** folder outside the source tree. We recommend to
create a folder ``build`` in the top-level directory, or multiple
folders in case you want to compile LAMMPS in different configurations
(``build-parallel``, ``build-serial``) or with different compilers
(``build-gnu``, ``build-clang``, ``build-intel``) and so on.


Running CMake
-------------

CLI version
^^^^^^^^^^^

In the (empty) ``build`` directory, we now run the command ``cmake
../cmake``, which will start the configuration phase and you will see
the progress of the configuration printed to the screen followed by a
summary of the enabled features, options and compiler settings. A typical
summary screen will look like this:

.. code-block::

   $ cmake ../cmake/
   -- The CXX compiler identification is GNU 8.2.0
   -- Check for working CXX compiler: /opt/tools/gcc-8.2.0/bin/c++
   -- Check for working CXX compiler: /opt/tools/gcc-8.2.0/bin/c++ - works
   -- Detecting CXX compiler ABI info
   -- Detecting CXX compiler ABI info - done
   -- Detecting CXX compile features
   -- Detecting CXX compile features - done
   -- Found Git: /usr/bin/git (found version "2.25.2") 
   -- Running check for auto-generated files from make-based build system
   -- Found MPI_CXX: /usr/lib64/mpich/lib/libmpicxx.so (found version "3.1") 
   -- Found MPI: TRUE (found version "3.1")  
   -- Looking for C++ include omp.h
   -- Looking for C++ include omp.h - found
   -- Found OpenMP_CXX: -fopenmp (found version "4.5") 
   -- Found OpenMP: TRUE (found version "4.5")  
   -- Found JPEG: /usr/lib64/libjpeg.so (found version "62") 
   -- Found PNG: /usr/lib64/libpng.so (found version "1.6.37") 
   -- Found ZLIB: /usr/lib64/libz.so (found version "1.2.11") 
   -- Found GZIP: /usr/bin/gzip  
   -- Found FFMPEG: /usr/bin/ffmpeg  
   -- Performing Test COMPILER_SUPPORTS-ffast-math
   -- Performing Test COMPILER_SUPPORTS-ffast-math - Success
   -- Performing Test COMPILER_SUPPORTS-march=native
   -- Performing Test COMPILER_SUPPORTS-march=native - Success
   -- Looking for C++ include cmath
   -- Looking for C++ include cmath - found
   -- Generating style_angle.h...
   [...]
   -- Generating lmpinstalledpkgs.h...
   -- The following tools and libraries have been found and configured:
    * Git
    * MPI
    * OpenMP
    * JPEG
    * PNG
    * ZLIB
   
   -- <<< Build configuration >>>
      Build type:       RelWithDebInfo
      Install path:     /home/akohlmey/.local
      Generator:        Unix Makefiles using /usr/bin/gmake
   -- <<< Compilers and Flags: >>>
   -- C++ Compiler:     /opt/tools/gcc-8.2.0/bin/c++
         Type:          GNU
         Version:       8.2.0
         C++ Flags:     -O2 -g -DNDEBUG
         Defines:       LAMMPS_SMALLBIG;LAMMPS_MEMALIGN=64;LAMMPS_JPEG;LAMMPS_PNG;LAMMPS_GZIP;LAMMPS_FFMPEG
         Options:       -ffast-math;-march=native
   -- <<< Linker flags: >>>
   -- Executable name:  lmp
   -- Static library flags:    
   -- <<< MPI flags >>>
   -- MPI includes:     /usr/include/mpich-x86_64
   -- MPI libraries:    /usr/lib64/mpich/lib/libmpicxx.so;/usr/lib64/mpich/lib/libmpi.so;
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /home/akohlmey/compile/lammps/build

The ``cmake`` command has one mandatory argument, and that is a folder
with either the file ``CMakeLists.txt`` or ``CMakeCache.txt``. The
``CMakeCache.txt`` file is created during the CMake configuration run
and contains all active settings, thus after a first run of CMake
all future runs in the build folder can use the folder ``.`` and CMake
will know where to find the CMake scripts and reload the settings
from the previous step.  This means, that one can modify an existing
configuration by re-running CMake, but only needs to provide flags
indicating the desired change, everything else will be retained. One
can also mix compilation and configuration, i.e. start with a minimal
configuration and then, if needed, enable additional features and
recompile.

The steps above **will NOT compile the code**\ . The compilation can be
started in a portable fashion with ``cmake --build .``, or you use the
selected built tool, e.g. ``make``.

TUI version
^^^^^^^^^^^

For the text mode UI CMake program the basical principle is the same.
You start the command ``ccmake ../cmake`` in the ``build`` folder.

.. list-table::

   * - .. figure:: JPG/ccmake-initial.png
          :target: JPG/ccmake-initial.png
          :align: center

          Initial ``ccmake`` screen

     - .. figure:: JPG/ccmake-config.png
          :target: JPG/ccmake-config.png
          :align: center

          Configure output of ``ccmake``

     - .. figure:: JPG/ccmake-options.png
          :target: JPG/ccmake-options.png
          :align: center

          Options screen of ``ccmake``

This will show you the initial screen (left image) with the empty
configuration cache. Now you type the 'c' key to run the configuration
step. That will do a first configuration run and show the summary
(center image). You exit the summary screen with 'e' and see now the
main screen with detected options and settings. You can now make changes
by moving and down with the arrow keys of the keyboard and modify
entries. For on/off settings, the enter key will toggle the state.
For others, hitting enter will allow you to modify the value and
you commit the change by hitting the enter key again or cancel using
the escape key.  All "new" settings will be marked with a star '\*'
and for as long as one setting is marked like this, you have to
re-run the configuration by hitting the 'c' key again, sometimes
multiple times unless the TUI shows the word "generate" next to the
letter 'g' and by hitting the 'g' key the build files will be written
to the folder and the TUI exits.  You can quit without generating
build files by hitting 'q'.

GUI version
^^^^^^^^^^^

For the graphical CMake program the steps are similar to the TUI
version.  You can type the command ``cmake-gui ../cmake`` in the
``build`` folder.  In this case the path to the CMake script folder is
not required, it can also be entered from the GUI.

.. list-table::

   * - .. figure:: JPG/cmake-gui-initial.png
          :target: JPG/cmake-gui-initial.png
          :align: center

          Initial ``cmake-gui`` screen

     - .. figure:: JPG/cmake-gui-popup.png
          :target: JPG/cmake-gui-popup.png
          :align: center

          Generator selection in ``cmake-gui``

     - .. figure:: JPG/cmake-gui-options.png
          :target: JPG/cmake-gui-options.png
          :align: center

          Options screen of ``cmake-gui``

Again, you start with an empty configuration cache (left image) and need
to start the configuration step.  For the very first configuration in a
folder, you will have a popup dialog (center image) asking to select the
desired build tool and some configuration settings (stick with the
default) and then you get the option screen with all new settings
highlighted in red.  You can modify them (or not) and click on the
"configure" button again until satisfied and click on the "generate"
button to write out the build files. You can exit the GUI from the
"File" menu or hit "ctrl-q".


Setting options
---------------




Using presets
-------------

Since LAMMPS has a lot of optional features specifying them all
on the command line, or - when selecting a different compiler toolchain -
multiple options have to be changed 


Choosing generators
-------------------


After the initial build, if you edit LAMMPS source files, or add your
own new files to the source directory, you can just re-type make from
your build directory and it will re-compile only the files that have
changed.  If you want to change CMake options you can run cmake (or
ccmake or cmake-gui) again from the same build directory and alter
various options; see details below.  Or you can remove the entire build
folder, recreate the directory and start over.

. Further details about features and settings for CMake
are in the `CMake online documentation <https://cmake.org/documentation/>`_
For the rest of this manual we will assume that the build environment
is generated for "Unix Makefiles" and thus the ``cmake --build .`` will
call the ``make`` command or you can use it directly.



There are 3 variants of the CMake command itself: a command-line version
(``cmake`` or ``cmake3``), a text mode UI version (``ccmake`` or ``ccmake3``),
and a graphical GUI version (``cmake-gui``).  You can use any of them
interchangeably to configure and create the LAMMPS build environment.
On Linux all the versions produce a Makefile as their output by default.
See more details on each below.

You can specify a variety of options with any of the 3 versions, which
affect how the build is performed and what is included in the LAMMPS
executable.  Links to pages explaining all the options are listed on
the :doc:`Build <Build>` doc page.

You must perform the CMake build system generation and compilation in
a new directory you create.  It can be anywhere on your local machine.
In these Build pages we assume that you are building in a directory
called ``lammps/build``.  You can perform separate builds independently
with different options, so long as you perform each of them in a
separate directory you create.  All the auxiliary files created by one
build process (executable, object files, log files, etc) are stored in
this directory or sub-directories within it that CMake creates.

----------

**Command-line version of CMake**\ :

.. code-block:: bash

   cmake  [options ...] /path/to/lammps/cmake  # build from any dir
   cmake  [options ...] ../cmake               # build from lammps/build
   cmake3 [options ...] ../cmake               # build from lammps/build

The cmake command takes one required argument, which is the LAMMPS
cmake directory which contains the CMakeLists.txt file.

The argument can be prefixed or followed by various CMake
command-line options.  Several useful ones are:

.. code-block:: bash

   -D CMAKE_INSTALL_PREFIX=path  # where to install LAMMPS executable/lib if desired
   -D CMAKE_BUILD_TYPE=type      # type = RelWithDebInfo (default), Release, MinSizeRel, or Debug
   -G output                     # style of output CMake generates (e.g. "Unix Makefiles" or "Ninja")
   -D CMAKE_MAKE_PROGRAM=builder # name of the builder executable (e.g. when using "gmake" instead of "make")
   -DVARIABLE=value              # setting for a LAMMPS feature to enable
   -D VARIABLE=value             # ditto, but cannot come after CMakeLists.txt dir
   -C path/to/preset/file        # load some CMake settings before configuring

All the LAMMPS-specific -D variables that a LAMMPS build supports are
described on the pages linked to from the :doc:`Build <Build>` doc page.
All of these variable names are upper-case and their values are
lower-case, e.g. -D LAMMPS_SIZES=smallbig.  For boolean values, any of
these forms can be used: yes/no, on/off, 1/0.

On Unix/Linux machines, CMake generates a Makefile by default to
perform the LAMMPS build.  Alternate forms of build info can be
generated via the -G switch, e.g. Visual Studio on a Windows machine,
Xcode on MacOS, or KDevelop on Linux.  Type ``cmake --help`` to see the
"Generator" styles of output your system supports.

.. note::

   When CMake runs, it prints configuration info to the screen.
   You should review this to verify all the features you requested were
   enabled, including packages.  You can also see what compilers and
   compile options will be used for the build.  Any errors in CMake
   variable syntax will also be flagged, e.g. mis-typed variable names or
   variable values.

CMake creates a CMakeCache.txt file when it runs.  This stores all the
settings, so that when running CMake again you can use the current
folder '.' instead of the path to the LAMMPS cmake folder as the
required argument to the CMake command. Either way the existing
settings will be inherited unless the CMakeCache.txt file is removed.

If you later want to change a setting you can rerun cmake in the build
directory with different setting. Please note that some automatically
detected variables will not change their value when you rerun cmake.
In these cases it is usually better to first remove all the
files/directories in the build directory, or start with a fresh build
directory.

----------

**Curses version (terminal-style menu) of CMake**\ :

.. code-block:: bash

   ccmake ../cmake

You initiate the configuration and build environment generation steps
separately. For the first you have to type **c**\ , for the second you
have to type **g**\ . You may need to type **c** multiple times, and may be
required to edit some of the entries of CMake configuration variables
in between.  Please see the `ccmake manual <https://cmake.org/cmake/help/latest/manual/ccmake.1.html>`_ for
more information.

----------

**GUI version of CMake**\ :

.. code-block:: bash

   cmake-gui ../cmake

You initiate the configuration and build environment generation steps
separately. For the first you have to click on the **Configure** button,
for the second you have to click on the **Generate** button.  You may
need to click on **Configure** multiple times, and may be required to
edit some of the entries of CMake configuration variables in between.
Please see the `cmake-gui manual <https://cmake.org/cmake/help/latest/manual/cmake-gui.1.html>`_
for more information.


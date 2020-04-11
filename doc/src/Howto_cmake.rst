Using CMake with LAMMPS tutorial
================================

Thus a configuration can be quickly modified by directing CMake to the
location of this cache file and then using options that are supposed to
be altered.


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

You must have CMake version 3.10 or later on your system to build
LAMMPS.  Installation instructions for CMake are below.

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


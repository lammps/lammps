Build LAMMPS with CMake
=======================

This page describes how to use CMake in general to build LAMMPS.
Details for specific compile time settings and options for optional
add-on packages are discussed with those packages.  Links to those
pages on the :doc:`Build overview <Build>` page.

The following text assumes some familiarity with CMake and focuses on
using the command line tool ``cmake`` and what settings are supported
for building LAMMPS.  A more detailed tutorial on how to use ``cmake``
itself, the text mode or graphical user interface, change the generated
output files for different build tools and development environments is
on a :doc:`separate page <Howto_cmake>`.

LAMMPS currently requires that CMake version 3.10 or later is available.

Advantages of using CMake
^^^^^^^^^^^^^^^^^^^^^^^^^

CMake is an alternative to compiling LAMMPS in the traditional way
through :doc:`(manually customized) makefiles <Build_make>` and a rather
recent addition to LAMMPS thanks to the efforts of Christoph Junghans
(LANL) and Richard Berger (Temple U).  Using CMake has multiple
advantages that are specifically helpful for people with limited
experience in compiling software or for people that want to modify or
extend LAMMPS.

- CMake can detect available hardware, tools, features, and libraries
  and adapt the LAMMPS build configuration accordingly.
- CMake can output files for different build tools and also can generate
  support files for use with popular integrated development environments
  (IDEs).
- CMake will build all components in a single build operation.
- CMake supports out-of-source compilation, so multiple configurations
  and settings with different choices of LAMMPS packages can be
  configured and built concurrently from the same source tree.
- CMake simplifies packaging of LAMMPS for Linux distributions,
  environment modules, or automated build tools like `Homebrew
  <https://brew.sh/>`_.

.. _cmake_build:

Getting started
^^^^^^^^^^^^^^^

Building LAMMPS with CMake is a two-step process.  First you use CMake
to create a build environment in a new directory.  For that purpose you
can use either the command-line utility ``cmake`` (or ``cmake3``), the
text-mode UI utility ``ccmake`` (or ``ccmake3``) or the graphical
utility ``cmake-gui``, or use them interchangeably.

Here is a minimal example using the command line version of CMake to
build LAMMPS with no add-on packages enabled and no customization:

.. code-block:: bash

   cd lammps                # change to the LAMMPS distribution directory
   mkdir build; cd build    # create and use a build directory
   cmake ../cmake           # configuration reading CMake scripts from ../cmake
   cmake --build .          # compilation (or type "make")

This will create and change into a folder called ``build``, then run the
configuration step to generate build files for the default build command
and then launch that build command to compile LAMMPS.  During the
configuration step CMake will try to detect whether support for MPI,
OpenMP, FFTW, gzip, JPEG, PNG, and ffmpeg are available and enable the
corresponding configuration settings.  The progress of this
configuration can be followed on the screen and a summary of selected
options and settings will be printed at the end.  The ``cmake --build
.`` command will launch the compilation, which, if successful, will
ultimately produce a library ``liblammps.a`` and the LAMMPS executable
``lmp`` inside the ``build`` folder.

If your machine has multiple CPU cores (most do these days), you can
speed this up by compiling sources in parallel with ``make -j N`` (with
N being the maximum number of concurrently executed tasks).  Also
installation of the `ccache <https://ccache.dev/>`_ (= Compiler Cache)
software may speed up repeated compilation significantly, e.g. during code
development.

After compilation, you may optionally install the LAMMPS executable into
your system with:

.. code-block:: bash

   make install    # optional, copy compiled files into installation location

This will install the LAMMPS executable and library, some tools (if configured)
and additional files like LAMMPS API headers, manpages, potential and force field
files.  The location of the installation tree defaults to ``${HOME}/.local``.



.. _cmake_options:

CMake configuration and build options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMake commands have one mandatory argument: a folder containing a
file called ``CMakeLists.txt`` (for LAMMPS it is located in the
``cmake`` folder) or a folder containing a file called
``CMakeCache.txt``, which is generated at the end of the CMake
configuration step. The cache file contains all current CMake settings.
Thus a configuration can be quickly modified by directing CMake to the
location of this cache file and then using options that are supposed to
be altered.

To modify settings, enable or disable features, you need to set *variables*
with either the *-D* command line flag or edit them .  This can be used several times.
There are 3 major command line options used to change settings during
CMake configuration.

Generating files for alternate build tools (e.g. Ninja) and project files
for IDEs like Eclipse, CodeBlocks, or Kate can be selected using the *-G*
command line flag.  A list of available settings is given when running
``cmake --help``. Further details about features and settings for CMake
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

.. note::

   To perform a CMake build, no packages can be installed or a build
   been previously attempted in the LAMMPS src directory by using ``make``
   commands to :doc:`perform a conventional LAMMPS build <Build_make>`.
   CMake detects if this is the case and generates an error, telling you
   to type ``make no-all purge`` in the src directory to un-install all
   packages.  The purge removes all the \*.h files auto-generated by
   make.

You must have CMake version 3.10 or later on your system to build
LAMMPS.  Installation instructions for CMake are below.

After the initial build, if you edit LAMMPS source files, or add your
own new files to the source directory, you can just re-type make from
your build directory and it will re-compile only the files that have
changed.  If you want to change CMake options you can run cmake (or
ccmake or cmake-gui) again from the same build directory and alter
various options; see details below.  Or you can remove the entire build
folder, recreate the directory and start over.

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

----------

**Installing CMake**

Check if your machine already has CMake installed:

.. code-block:: bash

   which cmake             # do you have it?
   which cmake3            # version 3 may have this name
   cmake --version         # what specific version you have

On clusters or supercomputers which use environment modules to manage
software packages, do this:

.. code-block:: bash

   module list            # is a module for cmake already loaded?
   module avail           # is a module for cmake available?
   module load cmake3     # load cmake module with appropriate name

Most Linux distributions offer pre-compiled cmake packages through
their package management system. If you do not have CMake or a new
enough version, you can download the latest version at
`https://cmake.org/download/ <https://cmake.org/download/>`_.
Instructions on how to install it on various platforms can be found
`on this page <https://cmake.org/install/>`_.

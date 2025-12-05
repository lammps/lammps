Build LAMMPS with CMake
-----------------------

This page describes how to use `CMake <https://cmake.org>`_ in general
to build LAMMPS.  Details for specific compile time settings and options
to enable and configure add-on packages are discussed with those
packages.  Links to those pages on the :doc:`Build overview <Build>`
page.

The following text assumes some familiarity with CMake and focuses on
using the command-line tool ``cmake`` and what settings are supported
for building LAMMPS.  A more detailed tutorial on how to use CMake
itself, the text mode or graphical user interface, to change the
generated output files for different build tools and development
environments is on a :doc:`separate page <Howto_cmake>`.

.. note::

   LAMMPS currently requires CMake version 3.20 or later.

.. warning::

   You **must not** mix the :doc:`traditional make based <Build_make>`
   LAMMPS build procedure with using CMake.  No packages may be
   installed or a build been previously attempted in the LAMMPS source
   directory by using ``make <machine>``.  CMake will detect if this is
   the case and generate an error.  To remove conflicting files from the
   ``src`` you can use the command ``make no-all purge`` which will
   uninstall all packages and delete all auto-generated files.


Advantages of using CMake
^^^^^^^^^^^^^^^^^^^^^^^^^

CMake is the preferred way of compiling LAMMPS in contrast to the legacy
build system based on GNU make and through :doc:`(manually customized)
makefiles <Build_make>`.  Using CMake has multiple advantages that are
specifically helpful for people with limited experience in compiling
software or for people that want to modify or extend LAMMPS.

- CMake can detect available hardware, tools, features, and libraries
  and adapt the LAMMPS default build configuration accordingly.
- CMake can generate files for different build tools and integrated
  development environments (IDE).
- CMake supports customization of settings with a command-line, text
  mode, or graphical user interface.  No manual editing of files,
  knowledge of file formats or complex command-line syntax is required.
- All enabled components are compiled in a single build operation.
- Automated dependency tracking for all files and configuration options.
- Support for true out-of-source compilation.  Multiple configurations
  and settings with different choices of LAMMPS packages, settings, or
  compilers can be configured and built concurrently from the same
  source tree.
- Simplified packaging of LAMMPS for Linux distributions, environment
  modules, or automated build tools like `Spack <https://spack.io>`_
  or `Homebrew <https://brew.sh/>`_.
- Integration of automated unit and regression testing.

.. _cmake_build:

Getting started
^^^^^^^^^^^^^^^

Building LAMMPS with CMake is a two-step process.  In the first step,
you use CMake to generate a build environment in a new directory.  For
that purpose you can use either the command-line utility ``cmake`` (or
``cmake3``), the text-mode UI utility ``ccmake`` (or ``ccmake3``) or the
graphical utility ``cmake-gui``, or use them interchangeably.  The
second step is then the compilation and linking of all objects,
libraries, and executables using the selected build tool.  Here is a
minimal example using the command-line version of CMake to build LAMMPS
with no add-on packages enabled and no customization:

.. code-block:: bash

   cd lammps                # change to the LAMMPS source distribution directory
   cmake -S cmake -B build  # configure the "build" folder with CMake scripts from "cmake"
   cmake --build build      # compilation (or type "make -C build")

This will create a folder called ``build``, then run the configuration
step to generate build files for the default build command and then
launch that build command to compile LAMMPS.  During the configuration
step CMake will try to detect whether support for MPI, OpenMP, FFTW,
gzip, JPEG, PNG, and ffmpeg are available and enable the corresponding
configuration settings.  The progress of this configuration can be
followed on the screen and a summary of selected options and settings
will be printed at the end.  The ``cmake --build build`` command will
launch the compilation, which, if successful, will ultimately produce a
library ``liblammps.a`` and the LAMMPS executable ``lmp`` inside the
``build`` folder.

Compilation can take a long time, since LAMMPS is a large project with
many features. If your machine has multiple CPU cores (most do these
days), you can speed this up by compiling sources in parallel with
adding ``--parallel N`` to the ``cmake`` command line (with *N* being
the maximum number of concurrently executed tasks).  Installation of the
`ccache <https://ccache.dev/>`_ (= Compiler Cache) software may speed up
repeated compilation even more, e.g. during code development, especially
when repeatedly switching between branches.

After the initial build, whenever you edit LAMMPS source files, enable
or disable packages, change compiler flags or build options, you must
re-compile and relink the LAMMPS executable with ``cmake --build build``
(or ``make -C build``).  If the compilation fails for some reason, try
running ``cmake build`` and then compile again. The included dependency
tracking should make certain that only the necessary subset of files is
re-compiled.  You can also delete compiled objects, libraries, and
executables with ``cmake --build build --target clean`` (or ``make -C
build clean``).

After compilation, you may optionally install the LAMMPS executable into
your system with:

.. code-block:: bash

   cmake --install build    # optional, copy compiled files into installation location

This will install the LAMMPS executable and library, some tools (if
configured) and additional files like LAMMPS API headers, manpages,
potential and force field files.  The location of the installation tree
defaults to ``${HOME}/.local``.

.. note::

   If you have set `-D CMAKE_INSTALL_PREFIX` to install LAMMPS into a
   system location on a Linux machine, you also have to run (as root)
   the `ldconfig` program to update the cache file for fast lookup of
   system shared libraries.

.. admonition:: Using the installed library
   :class: Hint

   The CMake installation functionality is an experimental work in
   progress and thus not without problems, especially when writing your
   own program that is trying to use the LAMMPS C++ classes directly.

   While there is a well-defined :ref:`C-language interface
   <lammps_c_api>` with the ``library.h`` header file, there is no
   equivalent for the C++ interface yet.  When installing LAMMPS,
   only the core header files are copied into the installation folder
   and thus only high-level access to C++ features is available.

   The following is a minimal CMake example file for using the installed
   LAMMPS package which represents the current state of development:

   .. code-block:: cmake

      cmake_minimum_required(VERSION 3.20)
      project(simpleCC CXX)
      # set this to the LAMMPS installation location
      if(NOT CMAKE_PREFIX_PATH)
        set(CMAKE_PREFIX_PATH $ENV{HOME}/.local)
      endif()
      find_package(LAMMPS REQUIRED)
      add_executable(simpleCC simple.cpp)
      target_link_libraries(simpleCC PRIVATE LAMMPS::LAMMPS)

   The ``CMAKE_PREFIX_PATH`` setting tells CMake where to find the
   generated CMake configuration files for the `find_package()
   <https://cmake.org/cmake/help/latest/command/find_package.html>`_
   CMake command.  You can also specify a required minimal version or
   version range.  For that a numeric representation in the "YYYY.MM.DD"
   format has to be used: the 10 September 2025 release thus becomes
   version 2025.09.10.  The include statements in the ``simple.cpp``
   source file have to be prefixed with ``lammps/`` as follows:

   .. code-block:: C++

      #include <lammps/lammps.h>
      #include <lammps/input.h>
      #include <lammps/atom.h>
      #include <lammps/library.h>

   Using the ``LAMMPS::LAMMPS`` target imported from the installed
   LAMMPS CMake configuration files should set up the include and linker
   flags and folders automatically.  Below is the output for an example
   session (note how it checks for and includes MPI and OpenMP support
   since that specific LAMMPS library was set up this way):

   .. code-block:: console

      $ cmake -S . -B build -D CMAKE_PREFIX_PATH=$HOME/Downloads/test-install
      -- The CXX compiler identification is GNU 15.2.1
      -- Detecting CXX compiler ABI info
      -- Detecting CXX compiler ABI info - done
      -- Check for working CXX compiler: /usr/lib64/ccache/c++ - skipped
      -- Detecting CXX compile features
      -- Detecting CXX compile features - done
      -- Found MPI_CXX: /usr/lib64/mpich/lib/libmpicxx.so (found version "4.1")
      -- Found MPI: TRUE (found version "4.1") found components: CXX
      -- Found OpenMP_CXX: -fopenmp (found version "4.5")
      -- Found OpenMP: TRUE (found version "4.5") found components: CXX
      -- Found LAMMPS: /home/akohlmey/Downloads/test-install/lib64/liblammps.so.0
      -- Configuring done (0.9s)
      -- Generating done (0.0s)
      -- Build files have been written to: /home/akohlmey/Downloads/test-simple/build
      $ cmake --build build
      [ 50%] Building CXX object CMakeFiles/simpleCC.dir/simple.cpp.o
      [100%] Linking CXX executable simpleCC
      [100%] Built target simpleCC

.. _cmake_options:

Configuration and build options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMake commands have one mandatory argument: a folder containing a
file called ``CMakeLists.txt`` (for LAMMPS it is located in the
``cmake`` folder, in that case the current working directory becomes
the build folder) or a build folder containing a file called
``CMakeCache.txt``, which is generated at the end of the CMake
configuration step.  The cache file contains all current CMake settings.
This is a "legacy mode" or running CMake and thus often found
when searching the web.  We recommend to use the ``-S`` and ``-B``
folders to explicitly set the path to the folder containing the
``CMakeLists.txt`` file and the build folder, respectively.

To modify settings, enable or disable features, you need to set
*variables* with either the ``-D`` command-line flag (``-D
VARIABLE1_NAME=value``) or change them in the text mode of the graphical
user interface.  The ``-D`` flag can be used several times in one command.

For your convenience, we provide :ref:`CMake presets <cmake_presets>`
that combine multiple settings to enable optional LAMMPS packages or use
a different compiler tool chain.  Those are loaded with the ``-C`` flag
(``-C ../cmake/presets/basic.cmake``).  This step would only be needed
once, as the settings from the preset files are stored in the
``CMakeCache.txt`` file. It is also possible to customize the build
by adding one or more ``-D`` flags to the CMake command.

Generating files for alternate build tools (e.g. Ninja) and project files
for IDEs like Eclipse, CodeBlocks, or Kate can be selected using the ``-G``
command-line flag.  A list of available generator settings for your
specific CMake version is given when running ``cmake --help``.

.. _cmake_multiconfig:

Multi-configuration build systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Throughout this manual, it is mostly assumed that LAMMPS is being built
on a Unix-like operating system with "make" as the underlying "builder",
since this is the most common case.  In this case the build
"configuration" is chose using ``-D CMAKE_BUILD_TYPE=<configuration>``
with ``<configuration>`` being one of "Release", "Debug",
"RelWithDebInfo", or "MinSizeRel".  Some build tools, however, can also
use or even require having a so-called multi-configuration build system
setup.  For a multi-configuration build, the built type (or
configuration) is selected at compile time using the same build
files. E.g.  with:

.. code-block:: bash

   cmake --build build-multi --config Release

In that case the resulting binaries are not in the build folder directly
but in subdirectories corresponding to the build type (i.e. Release in
the example from above).  Similarly, for running unit tests the
configuration is selected with the ``-C`` flag:

.. code-block:: bash

   ctest -C Debug

The CMake scripts in LAMMPS have basic support for being compiled using
a multi-config build system, but not all of it has been ported.  This is
in particular applicable to compiling packages that require additional
libraries that would be downloaded and compiled by CMake.  The
``windows.cmake`` preset file tries to keep track of which packages can
be compiled natively with the MSVC compilers out-of-the box.  Not all of
the external libraries are portable to Windows, either.


Installing CMake
^^^^^^^^^^^^^^^^

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
   module load cmake      # load cmake module with appropriate name

Most Linux distributions offer pre-compiled cmake packages through their
package management system. If you do not have CMake or a recent enough
version (Note: for CentOS 7.x you need to enable the EPEL repository),
you can download the latest version from `https://cmake.org/download/
<https://cmake.org/download/>`_.  Links to more details on CMake can
be found `on this page <https://cmake.org/resources/>`_.

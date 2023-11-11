Build LAMMPS with CMake
-----------------------

This page describes how to use `CMake <https://cmake.org>`_ in general
to build LAMMPS.  Details for specific compile time settings and options
to enable and configure add-on packages are discussed with those
packages.  Links to those pages on the :doc:`Build overview <Build>`
page.

The following text assumes some familiarity with CMake and focuses on
using the command line tool ``cmake`` and what settings are supported
for building LAMMPS.  A more detailed tutorial on how to use CMake
itself, the text mode or graphical user interface, to change the
generated output files for different build tools and development
environments is on a :doc:`separate page <Howto_cmake>`.

.. note::

   LAMMPS currently requires that CMake version 3.16 or later is available.

.. warning::

   You must not mix the :doc:`traditional make based <Build_make>`
   LAMMPS build procedure with using CMake.  No packages may be
   installed or a build been previously attempted in the LAMMPS source
   directory by using ``make <machine>``.  CMake will detect if this is
   the case and generate an error.  To remove conflicting files from the
   ``src`` you can use the command ``make no-all purge`` which will
   uninstall all packages and delete all auto-generated files.


Advantages of using CMake
^^^^^^^^^^^^^^^^^^^^^^^^^

CMake is an alternative to compiling LAMMPS in the traditional way
through :doc:`(manually customized) makefiles <Build_make>`.  Using
CMake has multiple advantages that are specifically helpful for
people with limited experience in compiling software or for people
that want to modify or extend LAMMPS.

- CMake can detect available hardware, tools, features, and libraries
  and adapt the LAMMPS default build configuration accordingly.
- CMake can generate files for different build tools and integrated
  development environments (IDE).
- CMake supports customization of settings with a command line, text
  mode, or graphical user interface.  No manual editing of files,
  knowledge of file formats or complex command line syntax is required.
- All enabled components are compiled in a single build operation.
- Automated dependency tracking for all files and configuration options.
- Support for true out-of-source compilation. Multiple configurations
  and settings with different choices of LAMMPS packages, settings, or
  compilers can be configured and built concurrently from the same
  source tree.
- Simplified packaging of LAMMPS for Linux distributions, environment
  modules, or automated build tools like `Homebrew <https://brew.sh/>`_.
- Integration of automated unit and regression testing (the LAMMPS side
  of this is still under active development).

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
minimal example using the command line version of CMake to build LAMMPS
with no add-on packages enabled and no customization:

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

Compilation can take a long time, since LAMMPS is a large project with
many features. If your machine has multiple CPU cores (most do these
days), you can speed this up by compiling sources in parallel with
``make -j N`` (with N being the maximum number of concurrently executed
tasks).  Installation of the `ccache <https://ccache.dev/>`_ (= Compiler
Cache) software may speed up repeated compilation even more, e.g. during
code development, especially when repeatedly switching between branches.

After the initial build, whenever you edit LAMMPS source files, enable
or disable packages, change compiler flags or build options, you must
re-compile and relink the LAMMPS executable with ``cmake --build .`` (or
``make``).  If the compilation fails for some reason, try running
``cmake .`` and then compile again. The included dependency tracking
should make certain that only the necessary subset of files is
re-compiled.  You can also delete compiled objects, libraries, and
executables with ``cmake --build . --target clean`` (or ``make clean``).

After compilation, you may optionally install the LAMMPS executable into
your system with:

.. code-block:: bash

   make install    # optional, copy compiled files into installation location

This will install the LAMMPS executable and library, some tools (if
configured) and additional files like LAMMPS API headers, manpages,
potential and force field files.  The location of the installation tree
defaults to ``${HOME}/.local``.

.. _cmake_options:

Configuration and build options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMake commands have one mandatory argument: a folder containing a
file called ``CMakeLists.txt`` (for LAMMPS it is located in the
``cmake`` folder) or a build folder containing a file called
``CMakeCache.txt``, which is generated at the end of the CMake
configuration step.  The cache file contains all current CMake settings.

To modify settings, enable or disable features, you need to set
*variables* with either the *-D* command line flag (``-D
VARIABLE1_NAME=value``) or change them in the text mode of the graphical
user interface.  The *-D* flag can be used several times in one command.

For your convenience, we provide :ref:`CMake presets <cmake_presets>`
that combine multiple settings to enable optional LAMMPS packages or use
a different compiler tool chain.  Those are loaded with the *-C* flag
(``-C ../cmake/presets/basic.cmake``).  This step would only be needed
once, as the settings from the preset files are stored in the
``CMakeCache.txt`` file. It is also possible to customize the build
by adding one or more *-D* flags to the CMake command line.

Generating files for alternate build tools (e.g. Ninja) and project files
for IDEs like Eclipse, CodeBlocks, or Kate can be selected using the *-G*
command line flag.  A list of available generator settings for your
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
configuration is selected with the *-C* flag:

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
<https://cmake.org/download/>`_.  Instructions on how to install it on
various platforms can be found `on this page
<https://cmake.org/install/>`_.

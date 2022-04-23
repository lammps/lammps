Build LAMMPS with make
======================

Building LAMMPS with traditional makefiles requires that you have a
``Makefile.<machine>`` file appropriate for your system in either the
``src/MAKE``, ``src/MAKE/MACHINES``, ``src/MAKE/OPTIONS``, or
``src/MAKE/MINE`` directory (see below).  It can include various options
for customizing your LAMMPS build with a number of global compilation
options and features.

Requirements
^^^^^^^^^^^^

Those makefiles are written for and tested with GNU make and may not
be compatible with other make programs.  In most cases, if the "make"
program is not GNU make, then there will be a GNU make program
available under the name "gmake".  If GNU make or a compatible make is
not available, you may have to first install it or switch to building
with :doc:`CMake <Build_cmake>`.  The makefiles of the traditional
make based build process and the scripts they are calling expect a few
additional tools to be available and functioning.

  * a working C/C++ compiler toolchain supporting the C++11 standard; on
    Linux these are often the GNU compilers. Some older compilers
    require adding flags like ``-std=c++11`` to enable the C++11 mode.
  * a Bourne shell compatible "Unix" shell program (often this is ``bash``)
  * a few shell utilities: ``ls``, ``mv``, ``ln``, ``rm``, ``grep``, ``sed``, ``tr``, ``cat``, ``touch``, ``diff``, ``dirname``
  * python (optional, required for ``make lib-<pkg>`` in the src folder).
    python scripts are currently tested with python 2.7 and 3.6. The procedure
    for :doc:`building the documentation <Build_manual>` requires python 3.5 or later.

Getting started
^^^^^^^^^^^^^^^

To include LAMMPS packages (i.e. optional commands and styles) you must
enable (or "install") them first, as discussed on the :doc:`Build
package <Build_package>` page.  If a packages requires (provided or
external) libraries, you must configure and build those libraries
**before** building LAMMPS itself and especially **before** enabling
such a package with ``make yes-<package>``.  :doc:`Building LAMMPS with
CMake <Build_cmake>` can automate much of this for many types of
machines, especially workstations, desktops, and laptops, so we suggest
you try it first when building LAMMPS in those cases.

The commands below perform a default LAMMPS build, producing the LAMMPS
executable ``lmp_serial`` and ``lmp_mpi`` in ``lammps/src``:

.. code-block:: bash

   cd lammps/src   # change to main LAMMPS source folder
   make serial     # build a serial LAMMPS executable using GNU g++
   make mpi        # build a parallel LAMMPS executable with MPI
   make            # see a variety of make options

Compilation can take a long time, since LAMMPS is a large project with
many features. If your machine has multiple CPU cores (most do these
days), you can speed this up by compiling sources in parallel with
``make -j N`` (with N being the maximum number of concurrently executed
tasks).  Also installation of the `ccache <https://ccache.dev/>`_ (=
Compiler Cache) software may speed up repeated compilation even more,
e.g. during code development.

After the initial build, whenever you edit LAMMPS source files, or add
or remove new files to the source directory (e.g. by installing or
uninstalling packages), you must re-compile and relink the LAMMPS
executable with the same ``make <machine>`` command.  The makefile's
dependency tracking should insure that only the necessary subset of
files are re-compiled.  If you change settings in the makefile, you have
to recompile *everything*.  To delete all objects you can use ``make
clean-<machine>``.

.. note::

   Before the actual compilation starts, LAMMPS will perform several
   steps to collect information from the configuration and setup that
   is then embedded into the executable.  When you build LAMMPS for
   the first time, it will also compile a tool to quickly assemble
   a list of dependencies, that are required for the make program to
   correctly detect which parts need to be recompiled after changes
   were made to the sources.

Customized builds and alternate makefiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``src/MAKE`` directory tree contains the ``Makefile.<machine>``
files included in the LAMMPS distribution.  Typing ``make example`` uses
``Makefile.example`` from one of those folders, if available.  Thus the
``make serial`` and ``make mpi`` lines above use
``src/MAKE/Makefile.serial`` and ``src/MAKE/Makefile.mpi``,
respectively.  Other makefiles are in these directories:

.. code-block:: bash

   OPTIONS      # Makefiles which enable specific options
   MACHINES     # Makefiles for specific machines
   MINE         # customized Makefiles you create (you may need to create this folder)

Simply typing ``make`` lists all the available ``Makefile.<machine>``
files with a single line description toward the end of the output.  A
file with the same name can appear in multiple folders (not a good
idea).  The order the directories are searched is as follows:
``src/MAKE/MINE``, ``src/MAKE``, ``src/MAKE/OPTIONS``,
``src/MAKE/MACHINES``.  This gives preference to a customized file you
put in ``src/MAKE/MINE``.  If you create your own custom makefile under
a new name, please edit the first line with the description and machine
name, so you will not confuse yourself, when looking at the machine
summary.

Makefiles you may wish to try include these (some require a package
first be installed).  Many of these include specific compiler flags
for optimized performance.  Please note, however, that some of these
customized machine Makefile are contributed by users.  Since both
compilers, OS configurations, and LAMMPS itself keep changing, their
settings may become outdated:

.. code-block:: bash

   make mac             # build serial LAMMPS on a Mac
   make mac_mpi         # build parallel LAMMPS on a Mac
   make intel_cpu       # build with the INTEL package optimized for CPUs
   make knl             # build with the INTEL package optimized for KNLs
   make opt             # build with the OPT package optimized for CPUs
   make omp             # build with the OPENMP package optimized for OpenMP
   make kokkos_omp      # build with the KOKKOS package for OpenMP
   make kokkos_cuda_mpi # build with the KOKKOS package for GPUs
   make kokkos_phi      # build with the KOKKOS package for KNLs

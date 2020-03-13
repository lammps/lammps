Build LAMMPS with make
======================

Building LAMMPS with traditional makefiles requires that you have a
Makefile."machine" file appropriate for your system in the src/MAKE,
src/MAKE/MACHINES, src/MAKE/OPTIONS, or src/MAKE/MINE directory (see
below).  It can include various options for customizing your LAMMPS
build with a number of global compilation options and features.

Those makefiles are written for and tested with GNU make and may not
be compatible with other make programs.  In most cases, if the "make"
program is not GNU make, then there will be a GNU make program
available under the name "gmake".  If GNU make or a compatible make is
not available, you may have to first install it or switch to building
with :doc:`CMake <Build_cmake>`.  The makefiles of the traditional
make based build process and the scripts they are calling expect a few
additional tools to be available and functioning.

  * a Bourne shell compatible "Unix" shell program (often this is bash)
  * a few shell utilities: ls, mv, ln, rm, grep, sed, tr, cat, touch, diff, dirname
  * python (optional, required for "make lib-XXX" in the src folder)

To include LAMMPS packages (i.e. optional commands and styles) you
must enable them first, as discussed on the :doc:`Build package
<Build_package>` doc page.  If a packages requires (provided or
external) libraries, you must configure and build those libraries
**before** building LAMMPS itself and especially **before** enabling
such a package with "make yes-<package>".  Building :doc:`LAMMPS
with CMake <Build_cmake>` can automate much of this for many types of
machines, especially workstations, desktops, and laptops, so we suggest
you try it first when building LAMMPS in those cases.

The commands below perform a default LAMMPS build, producing the LAMMPS
executable lmp_serial and lmp_mpi in lammps/src:

.. code-block:: bash

   cd lammps/src
   make serial     # build a serial LAMMPS executable
   make mpi        # build a parallel LAMMPS executable with MPI
   make            # see a variety of make options

This initial compilation can take a long time, since LAMMPS is a large
project with many features. If your machine has multiple CPU cores
(most do these days), using a command like "make -jN mpi" (with N =
the number of available CPU cores) can be much faster.  If you plan to
do development on LAMMPS or need to re-compile LAMMPS repeatedly, the
installation of the ccache (= Compiler Cache) software may speed up
compilation even more.

After the initial build, whenever you edit LAMMPS source files, or add
or remove new files to the source directory (e.g. by installing or
uninstalling packages), you must re-compile and relink the LAMMPS
executable with the same "make" command.  This makefiles dependencies
should insure that only the subset of files that need to be are
re-compiled.

.. note::

   Before the actual compilation starts, LAMMPS will perform several
   steps to collect information from the configuration and setup that
   is then embedded into the executable.  When you build LAMMPS for
   the first time, it will also compile a tool to quickly assemble
   a list of dependencies, that are required for the make program to
   correctly detect which parts need to be recompiled after changes
   were made to the sources.

----------

The lammps/src/MAKE tree contains the Makefile.machine files included
in the LAMMPS distribution.  Typing "make machine" uses
*Makefile.machine*\ .  Thus the "make serial" or "make mpi" lines above
use Makefile.serial and Makefile.mpi, respectively.  Other makefiles
are in these directories:

.. parsed-literal::

   OPTIONS      # Makefiles which enable specific options
   MACHINES     # Makefiles for specific machines
   MINE         # customized Makefiles you create (you may need to create this folder)

Typing "make" lists all the available Makefile.machine files.  A file
with the same name can appear in multiple folders (not a good idea).
The order the directories are searched is as follows: src/MAKE/MINE,
src/MAKE, src/MAKE/OPTIONS, src/MAKE/MACHINES.  This gives preference
to a customized file you put in src/MAKE/MINE.

Makefiles you may wish to try include these (some require a package
first be installed).  Many of these include specific compiler flags
for optimized performance.  Please note, however, that some of these
customized machine Makefile are contributed by users.  Since both
compilers, OS configurations, and LAMMPS itself keep changing, their
settings may become outdated:

.. code-block:: bash

   make mac             # build serial LAMMPS on a Mac
   make mac_mpi         # build parallel LAMMPS on a Mac
   make intel_cpu       # build with the USER-INTEL package optimized for CPUs
   make knl             # build with the USER-INTEL package optimized for KNLs
   make opt             # build with the OPT package optimized for CPUs
   make omp             # build with the USER-OMP package optimized for OpenMP
   make kokkos_omp      # build with the KOKKOS package for OpenMP
   make kokkos_cuda_mpi # build with the KOKKOS package for GPUs
   make kokkos_phi      # build with the KOKKOS package for KNLs

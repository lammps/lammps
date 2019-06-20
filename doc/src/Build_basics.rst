Basic build options
===================

The following topics are covered on this page, for building both with
CMake and make:

* :ref:`Serial vs parallel build <serial>`
* :ref:`Choice of compiler and compile/link options <compile>`
* :ref:`Build LAMMPS as an executable or a library <exe>`
* :ref:`Build the LAMMPS documentation <doc>`
* :ref:`Install LAMMPS after a build <install>`


----------


.. _serial:

Serial vs parallel build
-------------------------------------

LAMMPS can be built to run in parallel using the ubiquitous `MPI (message-passing interface) <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_
library.  Or it can built to run on a single processor (serial)
without MPI.  It can also be built with support for OpenMP threading
(see more discussion below).

**CMake variables**\ :


.. parsed-literal::

   -D BUILD_MPI=value        # yes or no, default is yes if CMake finds MPI, else no
   -D BUILD_OMP=value        # yes or no (default)
   -D LAMMPS_MACHINE=name    # name = mpi, serial, mybox, titan, laptop, etc
                             # no default value

The executable created by CMake (after running make) is lmp\_name.  If
the LAMMPS\_MACHINE variable is not specified, the executable is just
lmp.  Using BUILD\_MPI=no will produce a serial executable.

**Traditional make**\ :


.. parsed-literal::

   cd lammps/src
   make mpi                # parallel build, produces lmp_mpi using Makefile.mpi
   make serial             # serial build, produces lmp_serial using Makefile/serial
   make mybox          # uses Makefile.mybox to produce lmp_mybox

Serial build (see src/MAKE/Makefile.serial):


.. parsed-literal::

   MPI_INC =       -I../STUBS
   MPI_PATH =      -L../STUBS
   MPI_LIB =	-lmpi_stubs

For a parallel build, if MPI is installed on your system in the usual
place (e.g. under /usr/local), you do not need to specify the 3
variables MPI\_INC, MPI\_PATH, MPI\_LIB.  The MPI wrapper on the compiler
(e.g. mpicxx, mpiCC) knows where to find the needed include and
library files.  Failing this, these 3 variables can be used to specify
where the mpi.h file (MPI\_INC), and the MPI library files (MPI\_PATH)
are found, and the name of the library files (MPI\_LIB).

For a serial build, you need to specify the 3 variables, as shown
above.

For a serial LAMMPS build, use the dummy MPI library provided in
src/STUBS.  You also need to build the STUBS library for your platform
before making LAMMPS itself.  A "make serial" build does this for.
Otherwise, type "make mpi-stubs" from the src directory, or "make"
from the src/STUBS dir.  If the build fails, you will need to edit the
STUBS/Makefile for your platform.

The file STUBS/mpi.c provides a CPU timer function called MPI\_Wtime()
that calls gettimeofday() .  If your system doesn't support
gettimeofday() , you'll need to insert code to call another timer.
Note that the ANSI-standard function clock() rolls over after an hour
or so, and is therefore insufficient for timing long LAMMPS
simulations.

**CMake and make info**\ :

If you are installing MPI yourself, we recommend MPICH2 from Argonne
National Laboratory or OpenMPI.  MPICH can be downloaded from the
`Argonne MPI site <http://www.mcs.anl.gov/research/projects/mpich2/>`_.
OpenMPI can be downloaded from the `OpenMPI site <http://www.open-mpi.org>`_.  Other MPI packages should also work.
If you are running on a large parallel machine, your system admins or
the vendor should have already installed a version of MPI, which is
likely to be faster than a self-installed MPICH or OpenMPI, so find
out how to build and link with it.

The majority of OpenMP (threading) support in LAMMPS is provided by
the USER-OMP package; see the :doc:`Speed omp <Speed_omp>` doc page for
details. The USER-INTEL package also provides OpenMP support (it is
compatible with USER-OMP) and adds vectorization support when compiled
with the Intel compilers on top of that. Also, the KOKKOS package can
be compiled for using OpenMP threading.

However, there are a few commands in LAMMPS that have native OpenMP
support.  These are commands in the MPIIO, SNAP, USER-DIFFRACTION, and
USER-DPD packages.  In addition some packages support OpenMP threading
indirectly through the libraries they interface to: e.g. LATTE and
USER-COLVARS.  See the :doc:`Packages details <Packages_details>` doc
page for more info on these packages and the doc pages for their
respective commands for OpenMP threading info.

For CMake, if you use BUILD\_OMP=yes, you can use these packages and
turn on their native OpenMP support and turn on their native OpenMP
support at run time, by setting the OMP\_NUM\_THREADS environment
variable before you launch LAMMPS.

For building via conventional make, the CCFLAGS and LINKFLAGS
variables in Makefile.machine need to include the compiler flag that
enables OpenMP. For GNU compilers it is -fopenmp.  For (recent) Intel
compilers it is -qopenmp.  If you are using a different compiler,
please refer to its documentation.

.. _default-none-issues:

**OpenMP Compiler compatibility info**\ : 

Some compilers do not fully support the 'default(none)' directive
and others (e.g. GCC version 9 and beyond) may implement OpenMP 4.0
semantics, which are incompatible with the OpenMP 3.1 directives used
in LAMMPS (for maximal compatibility with compiler versions in use).
In those case, all 'default(none)' directives (which aid in detecting
incorrect and unwanted sharing) can be replaced with 'default(shared)'
while dropping all 'shared()' directives. The script
'src/USER-OMP/hack\_openmp\_for\_pgi\_gcc9.sh' can be used to automate
this conversion.


----------


.. _compile:

Choice of compiler and compile/link options
---------------------------------------------------------

The choice of compiler and compiler flags can be important for
performance.  Vendor compilers can produce faster code than
open-source compilers like GNU.  On boxes with Intel CPUs, we suggest
trying the `Intel C++ compiler <intel_>`_.

.. _intel: https://software.intel.com/en-us/intel-compilers



On parallel clusters or supercomputers which use "modules" for their
compile/link environments, you can often access different compilers by
simply loading the appropriate module before building LAMMPS.

**CMake variables**\ :


.. parsed-literal::

   -D CMAKE_CXX_COMPILER=name            # name of C++ compiler
   -D CMAKE_C_COMPILER=name              # name of C compiler
   -D CMAKE_Fortran_COMPILER=name        # name of Fortran compiler

   -D CMAKE_CXX_FLAGS=string             # flags to use with C++ compiler
   -D CMAKE_C_FLAGS=string               # flags to use with C compiler
   -D CMAKE_Fortran_FLAGS=string         # flags to use with Fortran compiler

By default CMake will use a compiler it finds and it will add
optimization flags appropriate to that compiler and any :doc:`accelerator packages <Speed_packages>` you have included in the build.

You can tell CMake to look for a specific compiler with these variable
settings.  Likewise you can specify the FLAGS variables if you want to
experiment with alternate optimization flags.  You should specify all
3 compilers, so that the small number of LAMMPS source files written
in C or Fortran are built with a compiler consistent with the one used
for all the C++ files:


.. parsed-literal::

   Building with GNU Compilers:
   cmake ../cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
   Building with Intel Compilers:
   cmake ../cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort
   Building with LLVM/Clang Compilers:
   cmake ../cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_Fortran_COMPILER=flang

.. note::

   When the cmake command completes, it prints info to the screen
   as to which compilers it is using, and what flags will be used in the
   compilation.  Note that if the top-level compiler is mpicxx, it is
   simply a wrapper on a real compiler.  The underlying compiler info is
   what will be listed in the CMake output.  You should check to insure
   you are using the compiler and optimization flags are the ones you
   want.

**Makefile.machine settings**\ :

Parallel build (see src/MAKE/Makefile.mpi):


.. parsed-literal::

   CC =		mpicxx
   CCFLAGS =	-g -O3
   LINK =		mpicxx
   LINKFLAGS =	-g -O

Serial build (see src/MAKE/Makefile.serial):


.. parsed-literal::

   CC =		g++
   CCFLAGS =	-g -O3
   LINK =		g++
   LINKFLAGS =	-g -O

The "compiler/linker settings" section of a Makefile.machine lists
compiler and linker settings for your C++ compiler, including
optimization flags.  You should always use mpicxx or mpiCC for
a parallel build, since these compiler wrappers will include
a variety of settings appropriate for your MPI installation.

.. note::

   If you build LAMMPS with any :doc:`accelerator packages <Speed_packages>` included, they have specific
   optimization flags that are either required or recommended for optimal
   performance.  You need to include these in the CCFLAGS and LINKFLAGS
   settings above.  For details, see the individual package doc pages
   listed on the :doc:`Speed packages <Speed_packages>` doc page.  Or
   examine these files in the src/MAKE/OPTIONS directory.  They
   correspond to each of the 5 accelerator packages and their hardware
   variants:


.. parsed-literal::

   Makefile.opt                   # OPT package
   Makefile.omp                   # USER-OMP package
   Makefile.intel_cpu             # USER-INTEL package for CPUs
   Makefile.intel_coprocessor     # USER-INTEL package for KNLs
   Makefile.gpu                   # GPU package
   Makefile.kokkos_cuda_mpi       # KOKKOS package for GPUs
   Makefile.kokkos_omp            # KOKKOS package for CPUs (OpenMP)
   Makefile.kokkos_phi            # KOKKOS package for KNLs (OpenMP)


----------


.. _exe:

Build LAMMPS as an executable or a library
----------------------------------------------------

LAMMPS can be built as either an executable or as a static or shared
library.  The LAMMPS library can be called from another application or
a scripting language.  See the :doc:`Howto couple <Howto_couple>` doc
page for more info on coupling LAMMPS to other codes.  See the
:doc:`Python <Python_head>` doc page for more info on wrapping and
running LAMMPS from Python via its library interface.

**CMake variables**\ :


.. parsed-literal::

   -D BUILD_EXE=value           # yes (default) or no
   -D BUILD_LIB=value           # yes or no (default)
   -D BUILD_SHARED_LIBS=value   # yes or no (default)

Setting BUILD\_EXE=no will not produce an executable.  Setting
BUILD\_LIB=yes will produce a static library named liblammps.a.
Setting both BUILD\_LIB=yes and BUILD\_SHARED\_LIBS=yes will produce a
shared library named liblammps.so.

**Traditional make**\ :


.. parsed-literal::

   cd lammps/src
   make machine               # build LAMMPS executable lmp_machine
   make mode=lib machine      # build LAMMPS static lib liblammps_machine.a
   make mode=shlib machine    # build LAMMPS shared lib liblammps_machine.so

The two library builds also create generic soft links, named
liblammps.a and liblammps.so, which point to the liblammps\_machine
files.

**CMake and make info**\ :

Note that for a shared library to be usable by a calling program, all
the auxiliary libraries it depends on must also exist as shared
libraries.  This will be the case for libraries included with LAMMPS,
such as the dummy MPI library in src/STUBS or any package libraries in
the lib/packages directory, since they are always built as shared
libraries using the -fPIC switch.  However, if a library like MPI or
FFTW does not exist as a shared library, the shared library build will
generate an error.  This means you will need to install a shared
library version of the auxiliary library.  The build instructions for
the library should tell you how to do this.

As an example, here is how to build and install the `MPICH library <mpich_>`_, a popular open-source version of MPI, distributed by
Argonne National Lab, as a shared library in the default
/usr/local/lib location:

.. _mpich: http://www-unix.mcs.anl.gov/mpi




.. parsed-literal::

   ./configure --enable-shared
   make
   make install

You may need to use "sudo make install" in place of the last line if
you do not have write privileges for /usr/local/lib.  The end result
should be the file /usr/local/lib/libmpich.so.


----------


.. _doc:

Build the LAMMPS documentation
----------------------------------------

**CMake variable**\ :


.. parsed-literal::

   -D BUILD_DOC=value       # yes or no (default)

This will create the HTML doc pages within the CMake build directory.
The reason to do this is if you want to "install" LAMMPS on a system
after the CMake build via "make install", and include the doc pages in
the install.

**Traditional make**\ :


.. parsed-literal::

   cd lammps/doc
   make html       # html doc pages
   make pdf        # single Manual.pdf file

This will create a lammps/doc/html dir with the HTML doc pages so that
you can browse them locally on your system.  Type "make" from the
lammps/doc dir to see other options.

.. note::

   You can also download a tarball of the documentation for the
   current LAMMPS version (HTML and PDF files), from the website
   `download page <http://lammps.sandia.gov/download.html>`_.


----------


.. _install:

Install LAMMPS after a build
------------------------------------------

After building LAMMPS, you may wish to copy the LAMMPS executable of
library, along with other LAMMPS files (library header, doc files) to
a globally visible place on your system, for others to access.  Note
that you may need super-user privileges (e.g. sudo) if the directory
you want to copy files to is protected.

**CMake variable**\ :


.. parsed-literal::

   cmake -D CMAKE_INSTALL_PREFIX=path [options ...] ../cmake
   make                        # perform make after CMake command
   make install                # perform the installation into prefix

**Traditional make**\ :

There is no "install" option in the src/Makefile for LAMMPS.  If you
wish to do this you will need to first build LAMMPS, then manually
copy the desired LAMMPS files to the appropriate system directories.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

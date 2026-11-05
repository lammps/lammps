Basic build options
===================

The following topics are covered on this page:

.. contents::
   :local:
   :depth: 2

----------

.. _build_serial:

Serial vs parallel build
------------------------

LAMMPS is written to use the ubiquitous `MPI (Message Passing Interface)
<https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ library API
for distributed memory parallel computation.  You need to have such a
library installed for building and running LAMMPS in parallel using a
domain decomposition parallelization.  It is compatible with the MPI
standard version 2.x and later.  LAMMPS can also be built into a
"serial" executable for use with a single processor using the bundled
``MPI STUBS`` library in the ``src/STUBS`` folder.

Independent of the distributed memory MPI parallelization, parts of
LAMMPS are also written with support for shared memory parallelization
using the `OpenMP <https://en.wikipedia.org/wiki/OpenMP>`_ threading
standard.  A more detailed discussion of that is below.

.. code-block:: bash

   -D BUILD_MPI=value        # on or off, default is on if CMake finds MPI
   -D BUILD_OMP=value        # on or off, default is on if a compatible
                             # compiler is detected
   -D LAMMPS_MACHINE=name    # name = mpi, serial, mybox, titan, laptop, etc
                             # no default value

The executable created by CMake (after compilation) is named ``lmp``
unless the ``LAMMPS_MACHINE`` option is set.  When setting
``LAMMPS_MACHINE=name``, the executable will be called ``lmp_name``.
Using ``BUILD_MPI=off`` will enforce building a serial executable using
the MPI STUBS library.

MPI and OpenMP support in LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are installing MPI yourself to build a parallel LAMMPS
executable, we recommend either MPICH or OpenMPI, which are regularly
used and tested with LAMMPS by the LAMMPS developers.  MPICH can be
downloaded from the `MPICH home page <https://www.mpich.org>`_, and
OpenMPI can be downloaded correspondingly from the `OpenMPI home page
<https://www.open-mpi.org>`_.  Other MPI packages should also work.  No
specific vendor provided and standard compliant MPI library is currently
known to be incompatible with LAMMPS.  If you are running on a large
parallel machine, your system admins or the vendor should have already
installed a version of MPI, which is likely to be faster than a
self-installed MPICH or OpenMPI, so you should study the provided
documentation to find out how to build and link with it.

.. note::

   Some Linux distributions (for example Ubuntu) build their MPI
   packages so that the MPI compiler wrappers add link-time
   optimization flags (``-flto``).  When building with CMake, these
   flags are then also used to compile LAMMPS, which can make the build
   fail at the link step, for example with the KOKKOS package and CUDA
   (errors mentioning ``fatbinData`` or ``lto-wrapper``).  To remove
   these flags, run CMake a second time in the same build folder with
   ``-D MPI_CXX_COMPILE_OPTIONS=""`` and then build again.  This
   setting has no effect in the first CMake run in a new build folder.

The majority of OpenMP (threading) support in LAMMPS is provided by the
``OPENMP`` package; see the :doc:`Speed_omp` page for details. The
``INTEL`` package also includes OpenMP threading (it is compatible with
``OPENMP`` and will usually fall back on styles from that package, if a
``INTEL`` does not exist) and adds vectorization support when compiled
with compatible compilers, in particular the Intel compilers on top of
OpenMP. Also, the ``KOKKOS`` package can be compiled to include OpenMP
threading.

In addition, there are a few commands in LAMMPS that have native OpenMP
support included as well.  These are commands in the ``ML-SNAP``,
``DIFFRACTION``, and ``DPD-REACT`` packages.  Furthermore, some packages
support OpenMP threading indirectly through the libraries they interface
to: e.g. ``KSPACE``, and ``COLVARS``.  See the :doc:`Packages details
<Packages_details>` page for more info on these packages, and the pages
for their respective commands for OpenMP threading info.

If you use ``-D BUILD_OMP=on`` when configuring a LAMMPS build, you can
use these packages **and** turn on their native OpenMP support at run time
by setting the ``OMP_NUM_THREADS`` environment variable before you
launch LAMMPS.

----------

.. _build_compiler:

Choice of compiler and compile/link options
---------------------------------------------------------

The choice of compiler and compiler flags can be important for maximum
performance.  Vendor provided compilers for a specific hardware can
produce faster code than open-source compilers like the GNU compilers.
On the most common x86 hardware, the most popular C++ compilers are
quite similar in their ability to optimize regular C/C++ source code at
high optimization levels.  When using the ``INTEL`` package, there is a
distinct advantage in using the `Intel C++ compiler <intel_>`_ due to
much improved vectorization through SSE and AVX instructions on
compatible hardware.  The source code in that package conditionally
includes compiler specific directives to enable these high degrees of
vectorization.  This may change over time as equivalent vectorization
directives are included into the OpenMP standard and other compilers
adopt them.

.. _intel: https://software.intel.com/en-us/intel-compilers

.. note::

   The Intel oneAPI compilers can be installed free of charge from
   Intel's package repositories, and we recommend using the latest
   release.  LAMMPS has been tested with the versions 2023.2.1,
   2024.2.1, 2025.3.3, and 2026.2.0 of the Intel oneAPI C++ compiler.
   The Intel compilers use the C++ library header files of the GCC
   version installed on the system, and older versions of the Intel
   compilers may not be able to process the header files of a newer GCC
   version.  For example, version 2023.2.1 cannot compile LAMMPS with
   the header files of GCC 16.  With CMake, a different GCC installation
   can be selected with ``-D
   CMAKE_CXX_FLAGS="--gcc-install-dir=/usr/lib/gcc/x86_64-redhat-linux/15"``
   (the path depends on the Linux distribution).

   The ``setvars.sh`` script of the Intel oneAPI installation sets up
   the environment for all installed oneAPI components.  If the Intel
   MPI library is installed, this includes the folder with its header
   files, which are then also used when LAMMPS is compiled with a
   different MPI library.  This mix can lead to errors when compiling,
   linking, or running LAMMPS.  So when using ``setvars.sh``, LAMMPS
   should be compiled with Intel MPI, or the Intel MPI package should
   not be installed.

   With the Intel compilers, some of the :ref:`unit tests <testing>`
   with the "unstable" label fail.  Those tests should be skipped with
   ``ctest -LE unstable``, or their failures can be ignored.

On parallel clusters or supercomputers which use "environment modules"
for their compile/link environments, you can often access different
compilers by simply loading the appropriate module before building
LAMMPS.

By default CMake will use the compiler it finds according to its
internal preferences, and it will add optimization flags
appropriate to that compiler and any :doc:`accelerator packages
<Speed_packages>` you have included in the build.  CMake will
check if the detected or selected compiler is compatible with the
C++ support requirements of LAMMPS and stop with an error, if this
is not the case.  A C++17 compatible compiler is required.

You can tell CMake to look for a specific compiler with setting
CMake variables (listed below) during configuration.  For a few
common choices, there are also presets in the ``cmake/presets``
folder.  You may also specify the corresponding ``CMAKE_*_FLAGS``
variables individually, if you want to experiment with alternate
optimization flags.  You should specify all 3 compilers, so that
the (few) LAMMPS source files written in C or Fortran are built
with a compiler consistent with the one used for the C++ files:

.. code-block:: bash

   -D CMAKE_CXX_COMPILER=name            # name of C++ compiler
   -D CMAKE_C_COMPILER=name              # name of C compiler
   -D CMAKE_Fortran_COMPILER=name        # name of Fortran compiler

   -D CMAKE_CXX_STANDARD=17              # put compiler in C++17 mode
   -D CMAKE_CXX_FLAGS=string             # flags to use with C++ compiler
   -D CMAKE_C_FLAGS=string               # flags to use with C compiler
   -D CMAKE_Fortran_FLAGS=string         # flags to use with Fortran compiler

A few example command lines are:

.. code-block:: bash

   # Building with GNU Compilers:
   cmake -S cmake -B build -D CMAKE_C_COMPILER=gcc \
         -D CMAKE_CXX_COMPILER=g++ -D CMAKE_Fortran_COMPILER=gfortran
   # Building with Intel oneAPI Compilers:
   cmake -S cmake -B build -D CMAKE_C_COMPILER=icx \
         -D CMAKE_CXX_COMPILER=icpx -D CMAKE_Fortran_COMPILER=ifx
   # Building with LLVM/Clang Compilers:
   cmake -S cmake -B build -D CMAKE_C_COMPILER=clang \
         -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_Fortran_COMPILER=flang
   # Building with PGI/Nvidia Compilers:
   cmake -S cmake -B build -D CMAKE_C_COMPILER=pgcc \
         -D CMAKE_CXX_COMPILER=pgc++ -D CMAKE_Fortran_COMPILER=pgfortran
   # Building with the NVHPC Compilers:
   cmake -S cmake -B build -D CMAKE_C_COMPILER=nvc \
         -D CMAKE_CXX_COMPILER=nvc++ -D CMAKE_Fortran_COMPILER=nvfortran

For compiling with the Clang/LLVM compilers a CMake preset is provided
that can be loaded with ``-C cmake/presets/clang.cmake``.  Similarly,
``-C cmake/presets/oneapi.cmake`` will switch to the LLVM based oneAPI
Intel compilers, ``-C cmake/presets/pgi.cmake`` will switch the compiler
to the PGI compilers, and ``-C cmake/presets/nvhpc.cmake`` will switch
to the NVHPC compilers.

.. note::

   When the ``cmake`` command completes, it prints a summary to the
   screen which compilers it is using and what flags and settings will
   be used for the compilation.  Note that if the top-level compiler is
   ``mpicxx``, it is simply a wrapper on a real compiler.  The
   underlying compiler info is what CMake will try to determine and
   report.  You should check to confirm you are using the compiler and
   optimization flags you want.

----------

.. _build_exe:

Build the LAMMPS executable and library
---------------------------------------

LAMMPS is always built as a library of C++ classes plus an executable.
The executable is a simple ``main()`` function that sets up MPI and then
creates a LAMMPS class instance from the LAMMPS library, which will then
process commands provided via a file or from the console input.  The
LAMMPS library can also be called from another application or a
scripting language.  See the :doc:`Howto couple <Howto_couple>` doc page
for more info on coupling LAMMPS to other codes.  See the :doc:`Python
<Python_head>` page for more info on wrapping and running LAMMPS from
Python via its library interface.

You can select through setting CMake variables between building a shared
or a static LAMMPS library and what kind of suffix is added to them (in
case you want to concurrently install multiple variants of binaries with
different settings). If none are set, defaults are applied.

.. code-block:: bash

   -D BUILD_SHARED_LIBS=value   # on or off (default)
   -D LAMMPS_MACHINE=name       # name = mpi, serial, mybox, titan, laptop, etc
                                # no default value

The compilation will always produce a LAMMPS library and an executable
linked to it.  By default, this will be a static library named
``liblammps.a`` and an executable named ``lmp``.  Setting
``BUILD_SHARED_LIBS=on`` will instead produce a shared library called
``liblammps.so`` (or ``liblammps.dylib`` or ``liblammps.dll`` depending
on the platform).  If ``LAMMPS_MACHINE=name`` is set in addition, the name
of the generated libraries will be changed to either
``liblammps_name.a`` or ``liblammps_name.so``\ , respectively and the
executable will be called ``lmp_name``.

Additional information
^^^^^^^^^^^^^^^^^^^^^^

Note that for creating a shared library, all the libraries it depends on
must be compiled to be compatible with shared libraries.  This should be
the case for code included with LAMMPS, such as the MPI stubs in
``src/STUBS`` or the libraries in the ``lib`` directory, since CMake
compiles them as position independent code by default (e.g. using the
``-fPIC`` compiler switch with GCC on Linux).
However, if an auxiliary library (like MPI or FFTW) does not exist as a
compatible format, the shared library linking step may generate an
error.  This means you will need to install a compatible version of the
auxiliary library.  The build instructions for that library should tell
you how to do this.

----------

.. _build_debug:

Including or removing debug support
-----------------------------------

By default the compilation settings will include the ``-g`` flag which
instructs the compiler to include debug information (e.g. which line of
source code a particular instruction correspond to).  This can be
extremely useful in case LAMMPS crashes and can help to provide crucial
information in :doc:`tracking down the origin of a crash <Errors_debug>`
and help the LAMMPS developers fix bugs in the source code.  However,
this increases the storage requirements for object files, libraries, and
the executable 3-5 fold.

If this is a concern, you can change the compilation settings to use
``-D CMAKE_BUILD_TYPE=Release`` or explicitly set the applicable
compiler flags or remove the debug information from the LAMMPS
executable.

To remove debug info from existing executables (e.g. after testing or
debugging), you can use the ``strip`` tool (e.g. ``strip -g lmp
liblammps.so.0``) to remove the debug information from the executable
files.

----------

.. _build_tools:

Build LAMMPS tools
------------------

Some tools described in :doc:`Auxiliary tools <Tools>` can be built directly.

.. code-block:: bash

   -D BUILD_TOOLS=value         # on or off (default). Build binary2txt,
                                # find-symmetry, msi2lmp, phana, reformat-json,
                                # stl_bin2txt, and (with Fortran)
                                # chain.x and micelle2d.x
   -D BUILD_LAMMPS_GUI=value    # on or off (default). Build LAMMPS-GUI
   -D BUILD_WHAM=value          # on (default). Download and build WHAM;
                                # only available for BUILD_LAMMPS_GUI=on

The generated binaries will also become part of the LAMMPS installation
(see below).

----------

Install LAMMPS after a build
------------------------------------------

After building LAMMPS, you may wish to copy the LAMMPS executable or
library, along with other LAMMPS files (library header, doc files), to a
globally visible place on your system, for others to access.  Note that
you may need super-user privileges (e.g. sudo) if the directory you want
to copy files to is protected.

.. code-block:: bash

   cmake -S cmake -B build -D CMAKE_INSTALL_PREFIX=path [options ...]
   cmake --build build           # compile LAMMPS
   cmake --install build         # perform the installation into prefix

During the installation process CMake will by default remove any runtime
path settings for loading shared libraries.  Because of this you may
have to set or modify the ``LD_LIBRARY_PATH`` (or ``DYLD_LIBRARY_PATH``)
environment variable, if you are installing LAMMPS into a non-system
location and/or are linking to libraries in a non-system location that
depend on such runtime path settings.
As an alternative, you may set the CMake variable ``LAMMPS_INSTALL_RPATH``
to ``on`` and then the runtime paths for any linked shared libraries
and the library installation folder for the LAMMPS library will be
embedded and thus the requirement to set environment variables is avoided.
The ``off`` setting is usually preferred for packaged binaries or when
setting up environment modules, the ``on`` setting is more convenient
for installing software into a non-system or personal folder.

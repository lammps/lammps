Basic build options
===================

The following topics are covered on this page, for building both with
CMake and make:

* :ref:`Serial vs parallel build <serial>`
* :ref:`Choice of compiler and compile/link options <compile>`
* :ref:`Build the LAMMPS executable and library <exe>`
* :ref:`Debug support <debug>`
* :ref:`Build the LAMMPS documentation <doc>`
* :ref:`Install LAMMPS after a build <install>`

----------

.. _serial:

Serial vs parallel build
------------------------

LAMMPS is written to use the ubiquitous `MPI (Message Passing Interface)
<https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ library API
for distributed memory parallel computation.  You need to have such a
library installed for building and running LAMMPS in parallel using a
domain decomposition parallelization.  It is compatible with the MPI
standard version 2.x and later.  LAMMPS can also be built into a
"serial" executable for use with a single processor using the bundled
MPI STUBS library.

Independent of the distributed memory MPI parallelization, parts of
LAMMPS are also written with support for shared memory parallelization
using the `OpenMP <https://en.wikipedia.org/wiki/OpenMP>`_ threading
standard. A more detailed discussion of that is below.

**CMake build**\ :

.. code-block:: bash

   -D BUILD_MPI=value        # yes or no, default is yes if CMake finds MPI, else no
   -D BUILD_OMP=value        # yes or no, default is yes if a compatible compiler is detected
   -D LAMMPS_MACHINE=name    # name = mpi, serial, mybox, titan, laptop, etc
                             # no default value

The executable created by CMake (after running make) is named ``lmp`` unless
the ``LAMMPS_MACHINE`` option is set.  When setting ``LAMMPS_MACHINE=name``
the executable will be called ``lmp_name``.  Using ``BUILD_MPI=no`` will
enforce building a serial executable using the MPI STUBS library.

**Traditional make**\ :

The build with traditional makefiles has to be done inside the source folder ``src``.

.. code-block:: bash

   make mpi                # parallel build, produces lmp_mpi using Makefile.mpi
   make serial             # serial build, produces lmp_serial using Makefile/serial
   make mybox              # uses Makefile.mybox to produce lmp_mybox

Any ``make machine`` command will look up the make settings from a file
``Makefile.machine`` in the folder ``src/MAKE`` or one of its
sub-directories ``MINE``, ``MACHINES``, or ``OPTIONS``, create a folder
``Obj_machine`` with all objects and generated files and an executable
called ``lmp_machine``\ .  The standard parallel build with ``make mpi``
assumes a standard MPI installation with MPI compiler wrappers where all
necessary compiler and linker flags to get access and link with the
suitable MPI headers and libraries are set by the wrapper programs.  For
other cases or the serial build, you have to adjust the make file
variables ``MPI_INC``, ``MPI_PATH``, ``MPI_LIB`` as well as ``CC`` and
``LINK``\ .  To enable OpenMP threading usually a compiler specific flag
needs to be added to the compile and link commands.  For the GNU
compilers, this is ``-fopenmp``\ , which can be added to the ``CC`` and
``LINK`` makefile variables.

For the serial build the following make variables are set (see src/MAKE/Makefile.serial):

.. code-block:: make

   CC =            g++
   LINK =          g++
   MPI_INC =       -I../STUBS
   MPI_PATH =      -L../STUBS
   MPI_LIB =       -lmpi_stubs

You also need to build the STUBS library for your platform before making
LAMMPS itself.  A ``make serial`` build does this for you automatically,
otherwise, type ``make mpi-stubs`` from the src directory, or ``make``
from the ``src/STUBS`` dir.  If the build fails, you may need to edit
the ``STUBS/Makefile`` for your platform.  The stubs library does not
provide MPI/IO functions required by some LAMMPS packages,
e.g. ``MPIIO`` or ``USER-LB``, and thus is not compatible with those
packages.

.. note::

   The file ``src/STUBS/mpi.c`` provides a CPU timer function called
   ``MPI_Wtime()`` that calls ``gettimeofday()``.  If your operating system
   does not support ``gettimeofday()``, you will need to insert code to
   call another timer.  Note that the ANSI-standard function ``clock()``
   rolls over after an hour or so, and is therefore insufficient for
   timing long LAMMPS simulations.

**MPI and OpenMP support info**\ :

If you are installing MPI yourself to build a parallel LAMMPS
executable, we recommend either MPICH or OpenMPI which are regularly
used and tested with LAMMPS by the LAMMPS developers.  MPICH can be
downloaded from the `MPICH home page <https://www.mpich.org>`_ and
OpenMPI can be downloaded correspondingly from the `OpenMPI home page
<https://www.open-mpi.org>`_.  Other MPI packages should also work.  No
specific vendor provided and standard compliant MPI library is currently
known to be incompatible with LAMMPS.  If you are running on a large
parallel machine, your system admins or the vendor should have already
installed a version of MPI, which is likely to be faster than a
self-installed MPICH or OpenMPI, so you should study the provided
documentation to find out how to build and link with it.

The majority of OpenMP (threading) support in LAMMPS is provided by the
``USER-OMP`` package; see the :doc:`Speed omp <Speed_omp>` doc page for
details. The ``USER-INTEL`` package also includes OpenMP threading (it
is compatible with ``USER-OMP`` and will usually fall back on styles
from that package, if a ``USER-INTEL`` does not exist) and adds
vectorization support when compiled with compatible compilers, in
particular the Intel compilers on top of OpenMP. Also, the ``KOKKOS``
package can be compiled to include OpenMP threading.

In addition, there are a few commands in LAMMPS that have native OpenMP
support included as well.  These are commands in the ``MPIIO``,
``SNAP``, ``USER-DIFFRACTION``, and ``USER-DPD`` packages.  In addition
some packages support OpenMP threading indirectly through the libraries
they interface to: e.g. ``LATTE``, ``KSPACE``, and ``USER-COLVARS``.
See the :doc:`Packages details <Packages_details>` doc page for more
info on these packages and the doc pages for their respective commands
for OpenMP threading info.

For CMake, if you use ``BUILD_OMP=yes``, you can use these packages
and turn on their native OpenMP support and turn on their native OpenMP
support at run time, by setting the ``OMP_NUM_THREADS`` environment
variable before you launch LAMMPS.

For building via conventional make, the ``CCFLAGS`` and ``LINKFLAGS``
variables in Makefile.machine need to include the compiler flag that
enables OpenMP. For GNU compilers it is ``-fopenmp``\ .  For (recent) Intel
compilers it is ``-qopenmp``\ .  If you are using a different compiler,
please refer to its documentation.

.. _default-none-issues:

**OpenMP Compiler compatibility info**\ :

Some compilers do not fully support the ``default(none)`` directive
and others (e.g. GCC version 9 and beyond) may implement OpenMP 4.0
semantics, which are incompatible with the OpenMP 3.1 semantics used
in LAMMPS (for maximal compatibility with compiler versions in use).
In those case, all ``default(none)`` directives (which aid in detecting
incorrect and unwanted sharing) can be replaced with ``default(shared)``
while dropping all ``shared()`` directives. The script
'src/USER-OMP/hack_openmp_for_pgi_gcc9.sh' can be used to automate
this conversion.

----------

.. _compile:

Choice of compiler and compile/link options
---------------------------------------------------------

The choice of compiler and compiler flags can be important for maximum
performance.  Vendor provided compilers for a specific hardware can
produce faster code than open-source compilers like the GNU compilers.
On the most common x86 hardware most popular C++ compilers are quite
similar in performance of C/C++ code at high optimization levels.  When
using the ``USER-INTEL`` package, there is a distinct advantage in using
the `Intel C++ compiler <intel_>`_ due to much improved vectorization
through SSE and AVX instructions on compatible hardware as the source
code includes changes and Intel compiler specific directives to enable
high degrees of vectorization.  This may change over time as equivalent
vectorization directives are included into OpenMP standard revisions and
other compilers adopt them.

.. _intel: https://software.intel.com/en-us/intel-compilers

On parallel clusters or supercomputers which use "environment modules"
for their compile/link environments, you can often access different
compilers by simply loading the appropriate module before building
LAMMPS.

**CMake build**\ :

By default CMake will use a compiler it finds according to internal
preferences and it will add optimization flags appropriate to that
compiler and any :doc:`accelerator packages <Speed_packages>` you have
included in the build.

You can tell CMake to look for a specific compiler with setting CMake
variable during configuration.  For a few common choices, there are also
presets in the ``cmake/presets`` folder.  For convenience, there is a
``CMAKE_TUNE_FLAGS`` variable that can be set to apply global compiler
options.  More on that below, but you can also specify the corresponding
``CMAKE_*_FLAGS`` variables individually if you want to experiment with
alternate optimization flags.  You should specify all 3 compilers, so
that the (few) LAMMPS source files written in C or Fortran are built
with a compiler consistent with the one used for the C++ files:

.. code-block:: bash

   -D CMAKE_CXX_COMPILER=name            # name of C++ compiler
   -D CMAKE_C_COMPILER=name              # name of C compiler
   -D CMAKE_Fortran_COMPILER=name        # name of Fortran compiler

   -D CMAKE_CXX_FLAGS=string             # flags to use with C++ compiler
   -D CMAKE_C_FLAGS=string               # flags to use with C compiler
   -D CMAKE_Fortran_FLAGS=string         # flags to use with Fortran compiler

A few example command lines are:

.. code-block:: bash

   # Building with GNU Compilers:
   cmake ../cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
   # Building with Intel Compilers:
   cmake ../cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort
   # Building with LLVM/Clang Compilers:
   cmake ../cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_Fortran_COMPILER=flang

For compiling with the Clang/LLVM compilers a CMake preset is provided that
can be loaded with `-C ../cmake/presets/clang.cmake`.  Similarly,
`-C ../cmake/presets/intel.cmake` should switch the 

In addition you can set ``CMAKE_TUNE_FLAGS`` to specifically add
compiler flags to tune for optimal performance on given hosts. By
default these are initialized to some compiler specific flags, to
optimize the LAMMPS executable with optimizations and instructions
available on the host where LAMMPS is compiled. For example, for Intel
compilers this would be ``-xHost`` and for GNU compilers this would be
``-march=native``. To turn these flags off, do ``-D CMAKE_TUNE_FLAGS=``.

.. note::

   When the cmake command completes, it prints a summary to the screen
   which compilers it is using and what flags and settings will be used
   for the  compilation.  Note that if the top-level compiler is mpicxx,
   it is simply a wrapper on a real compiler.  The underlying compiler
   info is what CMake will try to determine and report.  You should check
   to confirm you are using the compiler and optimization flags you want.

**Makefile.machine settings for traditional make**\ :

The "compiler/linker settings" section of a Makefile.machine lists
compiler and linker settings for your C++ compiler, including
optimization flags.  For a parallel build it is recommended to use
``mpicxx`` or ``mpiCC``, since these compiler wrappers will include a
variety of settings appropriate for your MPI installation and thus
avoiding the guesswork of finding the right flags.

Parallel build (see ``src/MAKE/Makefile.mpi``):

.. code-block:: bash

   CC =            mpicxx
   CCFLAGS =       -g -O3
   LINK =          mpicxx
   LINKFLAGS =     -g -O

Serial build with GNU gcc (see ``src/MAKE/Makefile.serial``):

.. code-block:: make

   CC =            g++
   CCFLAGS =       -g -O3
   LINK =          g++
   LINKFLAGS =     -g -O

.. note::

   If you build LAMMPS with any :doc:`accelerator packages <Speed_packages>`
   included, there may be specific optimization flags that are either
   required or recommended to enable required features and to achieve
   optimal performance.  You need to include these in the CCFLAGS and
   LINKFLAGS settings above.  For details, see the individual package
   doc pages listed on the :doc:`Speed packages <Speed_packages>` doc
   page.  Or examine these files in the src/MAKE/OPTIONS directory.
   They correspond to each of the 5 accelerator packages and their
   hardware variants:

.. code-block:: bash

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

Build the LAMMPS executable and library
---------------------------------------

LAMMPS is always built as a library of C++ classes plus an executable.
The executable is a simple ``main()`` function that sets up MPI and then
creates a LAMMPS class instance from the LAMMPS library, which
will then process commands provided via a file or from the console
input.  The LAMMPS library can also be called from another application
or a scripting language.  See the :doc:`Howto couple <Howto_couple>` doc
page for more info on coupling LAMMPS to other codes.  See the
:doc:`Python <Python_head>` doc page for more info on wrapping and
running LAMMPS from Python via its library interface.

**CMake build**\ :

For CMake builds, you can select through setting CMake variables between
building a shared or a static LAMMPS library and what kind of suffix is
added to them (in case you want to concurrently install multiple variants
of binaries with different settings). If none are set, defaults are applied.

.. code-block:: bash

   -D BUILD_SHARED_LIBS=value   # yes or no (default)
   -D LAMMPS_MACHINE=name       # name = mpi, serial, mybox, titan, laptop, etc
                                # no default value

The compilation will always produce a LAMMPS library and an executable
linked to it.  By default this will be a static library named
``liblammps.a`` and an executable named ``lmp`` Setting
``BUILD_SHARED_LIBS=yes`` will instead produce a shared library called
``liblammps.so`` (or ``liblammps.dylib`` or ``liblammps.dll`` depending
on the platform) If ``LAMMPS_MACHINE=name`` is set in addition, the name
of the generated libraries will be changed to either
``liblammps_name.a`` or ``liblammps_name.so``\ , respectively and the
executable will be called ``lmp_name``.

**Traditional make**\ :

With the traditional makefile based build process, the choice of
the generated executable or library depends on the "mode" setting.
Several options are available and ``mode=static`` is the default.

.. code-block:: bash

   make machine               # build LAMMPS executable lmp_machine
   make mode=static machine   # same as "make machine"
   make mode=shared machine   # build LAMMPS shared lib liblammps_machine.so instead

The "static" build will generate a static library called
``liblammps_machine.a`` and an executable named ``lmp_machine``\ , while
the "shared" build will generate a shared library
``liblammps_machine.so`` instead and ``lmp_machine`` will be linked to
it.  The build step will also create generic soft links, named
``liblammps.a`` and ``liblammps.so``\ , which point to the specific
``liblammps_machine.a/so`` files.

**CMake and make info**\ :

Note that for creating a shared library, all the libraries it depends on
must be compiled to be compatible with shared libraries.  This should be
the case for libraries included with LAMMPS, such as the dummy MPI
library in ``src/STUBS`` or any package libraries in the ``lib``
directory, since they are always built in a shared library compatible
way using the ``-fPIC`` compiler switch.  However, if an auxiliary
library (like MPI or FFTW) does not exist as a compatible format, the
shared library linking step may generate an error.  This means you will
need to install a compatible version of the auxiliary library.  The
build instructions for that library should tell you how to do this.

As an example, here is how to build and install the `MPICH library
<mpich_>`_, a popular open-source version of MPI, as a shared library
in the default /usr/local/lib location:

.. _mpich: https://www.mpich.org

.. code-block:: bash

   ./configure --enable-shared
   make
   make install

You may need to use ``sudo make install`` in place of the last line if
you do not have write privileges for ``/usr/local/lib`` or use the
``--prefix`` configuration option to select an installation folder,
where you do have write access.  The end result should be the file
``/usr/local/lib/libmpich.so``.  On many Linux installations the folder
``${HOME}/.local`` is an alternative to using ``/usr/local`` and does
not require superuser or sudo access.  In that case the configuration
step becomes:

.. code-block:: bash

  ./configure --enable-shared --prefix=${HOME}/.local

Avoiding to use "sudo" for custom software installation (i.e. from source
and not through a package manager tool provided by the OS) is generally
recommended to ensure the integrity of the system software installation.

----------

.. _debug:

Debug support
-------------

By default the compilation settings will include the *-g* flag which
instructs the compiler to include debug information (e.g. which line of
source code particular instructions correspond to).  This can be
extremely useful in case LAMMPS crashes and can help to provide crucial
information in :doc:`tracking down the origin of a crash <Errors_debug>`
and possibly help fix a bug in the source code.  However, this increases
the storage requirements for object files, libraries, and the executable
3-5 fold.  If this is a concern, you can change the compilation settings
(either by editing the machine makefile or setting the compiler flags or
build time when using CMake).  If you are only concerned about the
executable being too large, you can use the ``strip`` tool (e.g. ``strip
lmp_serial``) to remove the debug information from the file.

----------

.. _doc:

Build the LAMMPS documentation
----------------------------------------

The LAMMPS manual is written in `reStructuredText <rst_>`_ format which
can be translated to different output format using the `Sphinx <sphinx_>`_
document generator tool.  Currently the translation to HTML and PDF (via
LaTeX) are supported.  For that to work a Python 3 interpreter and
internet access is required.  For the documentation build a python
based virtual environment is set up in the folder doc/docenv and various
python packages are installed into that virtual environment via the pip
tool.  The actual translation is then done via make commands.

.. _rst: https://docutils.readthedocs.io/en/sphinx-docs/user/rst/quickstart.html
.. _sphinx: https://sphinx-doc.org

**Documentation make option**\ :

The following make commands can be issued in the doc folder of the
LAMMPS source distribution.

.. code-block:: bash

  make html          # create HTML doc pages in html directory
  make pdf           # create Developer.pdf and Manual.pdf in this directory
  make fetch         # fetch HTML and PDF files from LAMMPS web site
  make clean         # remove all intermediate files
  make clean-all     # reset the entire doc build environment
  make anchor_check  # scan for duplicate anchor labels
  make style_check   # check for complete and consistent style lists
  make package_check # check for complete and consistent package lists
  make spelling      # spell-check the manual

Thus "make html" will create a "doc/html" directory with the HTML format
manual pages so that you can browse them with a web browser locally on
your system.

.. note::

   You can also download a tarball of the documentation for the
   current LAMMPS version (HTML and PDF files), from the website
   `download page <https://lammps.sandia.gov/download.html>`_.

**CMake build option**\ :

It is also possible to create the HTML version of the manual within
the :doc:`CMake build directory <Build_cmake>`.  The reason for this
option is to include the installation of the HTML manual pages into
the "install" step when installing LAMMPS after the CMake build via
``make install``.

.. code-block:: bash

   -D BUILD_DOC=value       # yes or no (default)

----------

.. _tools:

Build LAMMPS tools
------------------------------

Some tools described in :doc:`Auxiliary tools <Tools>` can be built directly
using CMake or Make.

**CMake build3**\ :

.. code-block:: bash

   -D BUILD_TOOLS=value       # yes or no (default)

The generated binaries will also become part of the LAMMPS installation
(see below).

**Traditional make**\ :

.. code-block:: bash

   cd lammps/tools
   make all              # build all binaries of tools
   make binary2txt       # build only binary2txt tool
   make chain            # build only chain tool
   make micelle2d        # build only micelle2d tool
   make thermo_extract   # build only thermo_extract tool

----------

.. _install:

Install LAMMPS after a build
------------------------------------------

After building LAMMPS, you may wish to copy the LAMMPS executable of
library, along with other LAMMPS files (library header, doc files) to
a globally visible place on your system, for others to access.  Note
that you may need super-user privileges (e.g. sudo) if the directory
you want to copy files to is protected.

**CMake build**\ :

.. code-block:: bash

   cmake -D CMAKE_INSTALL_PREFIX=path [options ...] ../cmake
   make                        # perform make after CMake command
   make install                # perform the installation into prefix

**Traditional make**\ :

There is no "install" option in the ``src/Makefile`` for LAMMPS.  If
you wish to do this you will need to first build LAMMPS, then manually
copy the desired LAMMPS files to the appropriate system directories.

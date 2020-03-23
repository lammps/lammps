Link LAMMPS as a library to another code
========================================

LAMMPS can be used as a library by another application, including
Python scripts.  The files ``src/library.cpp`` and ``library.h``
define the C-style API for using LAMMPS as a library.  See the
:doc:`Howto library <Howto_library>` page for a description of the
interface and how to use it for your needs.

The :doc:`Build basics <Build_basics>` doc page explains how to build
LAMMPS as either a shared or static library.  This results in one of
these 2 files:

.. code-block:: bash

   liblammps.so      # shared library
   liblammps.a       # static library

.. note::

   Care should be taken to use the same MPI library for the calling
   code and the LAMMPS library.  The library.h file includes mpi.h and
   uses definitions from it so those need to be available and
   consistent.  When LAMMPS is compiled with the MPI STUBS library,
   then its mpi.h file needs to be included.  While it is technically
   possible to use a full MPI library in the calling code and link to
   a serial LAMMPS library compiled with MPI STUBS, it is recommended
   to use the *same* MPI library for both, and then use MPI_Comm_split()
   in the calling code to pass a suitable communicator with a subset
   of MPI ranks to the function creating the LAMMPS instance.

----------

**Link with LAMMPS as a static library**\ :

The calling application can link to LAMMPS as a static library with
compilation and link commands as in the examples shown below.  These
are examples for a code written in C in the file ``caller.c``.
The benefit of linking to a static library is, that the resulting
executable is independent of that library since all required
executable code from the library is copied into the calling executable.

*CMake build*\ :

This assumes that LAMMPS has been configured with ``-D BUILD_LIB=yes``
and installed with ``make install`` and the ``PKG_CONFIG_PATH`` environment
variable updated to include the ``liblammps.pc`` file installed into the
configured destination folder, if needed.  The commands to compile and
link the coupled executable are then:

.. code-block:: bash

   mpicc -c -O $(pkgconf liblammps --cflags) caller.c
   mpicxx -o caller caller.o -$(pkgconf liblammps --libs)

*Traditional make*\ :

This assumes that LAMMPS has been compiled in the folder
``${HOME}/lammps/src`` with ``make mode=lib mpi``. The commands to
compile and link the coupled executable are then:

.. code-block:: bash

   mpicc -c -O -I${HOME}/lammps/src caller.c
   mpicxx -o caller caller.o -L${HOME}/lammps/src -llammps

The ``-I`` argument is the path to the location of the ``library.h``
header file containing the interface to the LAMMPS C-style library
interface.  The ``-L`` argument is the path to where the ``liblammps.a``
file is located.  The ``-llammps`` argument is shorthand for telling the
compiler to link the file ``liblammps.a``.

However, it is only as simple as shown above for the case of a plain
LAMMPS library without any optional packages that depend on libraries
(bundled or external).  Otherwise, you need to include all flags,
libraries, and paths for the coupled executable, that are also
required to link the LAMMPS executable.

*CMake build*\ :

When using CMake, additional libraries with sources in the lib folder
are built, but not included in ``liblammps.a`` and (currently) not
installed with ``make install`` and not included in the ``pkgconfig``
configuration file.  They can be found in the top level build folder,
but you have to determine the necessary link flags manually.  It is
therefore recommended to either use the traditional make procedure to
build and link with a static library or build and link with a shared
library instead.

*Traditional make*\ :

After you have compiled a static LAMMPS library using the conventional
build system for example with ``make mode=lib serial``. And you also
have installed the POEMS package after building its bundled library in
lib/poems. Then the commands to build and link the coupled executable
change to:

.. code-block:: bash

   gcc -c -O -I${HOME}/lammps/src/STUBS -I${HOME}/lammps/src -caller.c
   g++ -o caller caller.o -L${HOME}/lammps/lib/poems \
     -L${HOME}/lammps/src/STUBS -L${HOME}/lammps/src -llammps -lpoems -lmpi_stubs

Note, that you need to link with ``g++`` instead of ``gcc``, since LAMMPS
is C++ code.  You can display the currently applied settings for building
LAMMPS for the "serial" machine target by using the command:

.. code-block:: bash

   make mode=print serial

Which should output something like:

.. code-block:: bash

   # Compiler:
   CXX=g++
   # Linker:
   LD=g++
   # Compilation:
   CXXFLAGS=-g -O3 -DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64 -I${HOME}/lammps/lib/poems -I${HOME}/lammps/src/STUBS
   # Linking:
   LDFLAGS=-g -O
   # Libraries:
   LDLIBS=-L${HOME}/lammps/lib/poems -L${HOME}/lammps/src/STUBS -lpoems -lmpi_stubs

From this you can gather the necessary paths and flags.  With
makefiles for other *machine* configurations you need to do the
equivalent and replace *serial* with the corresponding *machine* name
of the makefile.

----------

**Link with LAMMPS as a shared library**\ :

When linking to LAMMPS built as a shared library, the situation
becomes much simpler, as all dependent libraries and objects are
included in the shared library, which is - technically speaking -
effectively a regular LAMMPS executable that is missing the `main()`
function.  Thus those libraries need not to be specified when linking
the calling executable.  Only the ``-I`` flags are needed.  So the
example case from above of the serial version static LAMMPS library
with the POEMS package installed becomes:

*CMake build*\ :

The commands with a shared LAMMPS library compiled with the CMake
build process are the same as for the static library.

.. code-block:: bash

   mpicc -c -O $(pkgconf liblammps --cflags) caller.c
   mpicxx -o caller caller.o -$(pkgconf --libs)

*Traditional make*\ :

The commands with a shared LAMMPS library compiled with the
traditional make build using ``make mode=shlib serial`` becomes:

.. code-block:: bash

   gcc -c -O -I${HOME}/lammps/src/STUBS -I${HOME}/lammps/src -caller.c
   g++ -o caller caller.o -L${HOME}/lammps/src -llammps

*Locating liblammps.so at runtime*\ :

However, now the ``liblammps.so`` file is required at runtime and needs
to be in a folder, where the shared linker program of the operating
system can find it.  This would be either a folder like ``/usr/local/lib64``
or ``${HOME}/.local/lib64`` or a folder pointed to by the ``LD_LIBRARY_PATH``
environment variable. You can type

.. code-block:: bash

   printenv LD_LIBRARY_PATH

to see what directories are in that list.

Or you can add the LAMMPS src directory (or the directory you performed
a CMake style build in) to your ``LD_LIBRARY_PATH``, so that the current
version of the shared library is always available to programs that use it.

For the Bourne or Korn shells (/bin/sh, /bin/ksh, /bin/bash etc.), you
would add something like this to your ``${HOME}/.profile`` file:

.. code-block:: bash

   LD_LIBRARY_PATH ${LD_LIBRARY_PATH-/usr/lib64}:${HOME}/lammps/src
   export LD_LIBRARY_PATH

For the csh or tcsh shells, you would equivalently add something like this
to your ``${HOME}/.cshrc`` file:

.. code-block:: csh

   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HOME}/lammps/src

You can verify whether all required shared libraries are found with the
`ldd` tool.  Example:

.. code-block:: bash

   $ LD_LIBRARY_PATH=/home/user/lammps/src ldd caller
        linux-vdso.so.1 (0x00007ffe729e0000)
        liblammps.so => /home/user/lammps/src/liblammps.so (0x00007fc91bb9e000)
        libstdc++.so.6 => /lib64/libstdc++.so.6 (0x00007fc91b984000)
        libm.so.6 => /lib64/libm.so.6 (0x00007fc91b83e000)
        libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00007fc91b824000)
        libc.so.6 => /lib64/libc.so.6 (0x00007fc91b65b000)
        /lib64/ld-linux-x86-64.so.2 (0x00007fc91c094000)

If a required library is missing, you would get a 'not found' entry:

.. code-block:: bash

   $  ldd caller
        linux-vdso.so.1 (0x00007ffd672fe000)
        liblammps.so => not found
        libstdc++.so.6 => /usr/lib64/libstdc++.so.6 (0x00007fb7c7e86000)
        libm.so.6 => /usr/lib64/libm.so.6 (0x00007fb7c7d40000)
        libgcc_s.so.1 => /usr/lib64/libgcc_s.so.1 (0x00007fb7c7d26000)
        libc.so.6 => /usr/lib64/libc.so.6 (0x00007fb7c7b5d000)
        /lib64/ld-linux-x86-64.so.2 (0x00007fb7c80a2000)


Link LAMMPS as a library to another code
========================================

LAMMPS can be used as a library by another application, including
Python scripts.  The files src/library.cpp and library.h define the
C-style API for using LAMMPS as a library.  See the :doc:`Howto library <Howto_library>` doc page for a description of the
interface and how to extend it for your needs.

The :doc:`Build basics <Build_basics>` doc page explains how to build
LAMMPS as either a shared or static library.  This results in one of
these 2 files:

.. parsed-literal::

   liblammps.so      # shared library
   liblammps.a       # static library


----------


**Link with LAMMPS as a static library**\ :

The calling application can link to LAMMPS as a static library with
a compilation and link command like this (assuming a code written in
C in the file *caller.c*):

.. code-block:: bash

   mpicc -c -O -I${HOME}/lammps/src caller.c
   mpicxx -o caller caller.o -L${HOME}/lammps/src -llammps

The *-I* argument is the path to the location of the *library.h*
header file containing the interface to the LAMMPS C-style library
interface.  The *-L* argument is the path to where the *liblammps.a*
file is located.  The *-llammps* argument is shorthand for telling the
compiler to link the file *liblammps.a*\ .

The benefit of linking ot a static library is, that the resulting
executable is independent of that library since all used executable
code is copied into the calling executable.  However, it is only as
simple as shown for the case of a plain LAMMPS library without any
optional packages and libraries.  Otherwise, you need to include all
flags, libraries, and paths that are required to link the LAMMPS
executable.  Assuming you have compiled LAMMPS using the conventional
build system with "make serial" and also have the POEMS package
installed, the command changes to:

.. code-block:: bash

   gcc -c -O -I${HOME}/lammps/src/STUBS -I${HOME}/lammps/src -caller.c
   mpicxx -o caller caller.o -L${HOME}/lammps/src/../lib/poems \
     -L${HOME}/lammps/src/STUBS -L${HOME}/lammps/src -llammps -lpoems -lmpi_stubs 

You can display the currently applied settings for the "serial" machine
target by using the command:

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

----------

**Link with LAMMPS as a shared library**\ :

When linking to a shared library, the situation becomes much simpler,
as all dependent libraries and objects are included in the shared
library, which is - technically speaking - very similar to a regular
LAMMPS executable that is missing the `main()` function. Thus those
libraries need not to be specified when linking the calling
executable.  So the example case from above of the serial version
library with the POEMS package installed becomes:

.. code-block:: bash

   gcc -c -O -I${HOME}/lammps/src/STUBS -I${HOME}/lammps/src -caller.c
   mpicxx -o caller caller.o -L${HOME}/lammps/src -llammps

However, now the `liblammps.so` file is required at runtime and needs
to be in a folder, where the shared linker program of the operating
system can find it.  This would be either a folder like "/usr/local/lib64"
or "${HOME}/.local/lib64" or a folder pointed to by the LD\_LIBRARY\_PATH
environment variable. You can type

.. code-block:: bash

   printenv LD_LIBRARY_PATH

to see what directories are in that list.

Or you can add the LAMMPS src directory (or the directory you performed
a CMake style build in) to your LD\_LIBRARY\_PATH, so that the current
version of the shared library is always available to programs that use it.

For the Bourne or Korn shells (/bin/sh, /bin/ksh, /bin/bash etc.), you
would add something like this to your ~/.profile file:

.. code-block:: bash

   LD_LIBRARY_PATH ${LD_LIBRARY_PATH-/usr/lib64}:${HOME}/lammps/src
   export LD_LIBRARY_PATH

For the csh or tcsh shells, you would equivalently add something like this
to your ~/.cshrc file:


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


----------


**Calling the LAMMPS library**\ :

Either flavor of library (static or shared) allows one or more LAMMPS
objects to be instantiated from the calling program. When used from a
C++ program, most of the symbols and functions in LAMMPS are wrapped
in a LAMMPS\_NS namespace; you can safely use any of its classes and
methods from within the calling code, as needed, and you will not incur
conflicts with functions and variables in your code that share the name.
This, however, does not extend to all additional libraries bundled with
LAMMPS in the lib folder and some of the low-level code of some packages.

To be compatible with C, Fortran, Python programs, the library has a simple
C-style interface, provided in src/library.cpp and src/library.h.

See the :doc:`Python library <Python_library>` doc page for a
description of the Python interface to LAMMPS, which wraps the C-style
interface from a shared library through the ctypes python module.

See the sample codes in examples/COUPLE/simple for examples of C++ and
C and Fortran codes that invoke LAMMPS through its library interface.
Other examples in the COUPLE directory use coupling ideas discussed on
the :doc:`Howto couple <Howto_couple>` doc page.

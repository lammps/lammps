Link LAMMPS as a library to another code
========================================

LAMMPS can be used as a library by another application, including
Python scripts.  The files src/library.cpp and library.h define the
C-style API for using LAMMPS as a library.  See the :doc:`Howto library <Howto_library>` doc page for a description of the
interface and how to extend it for your needs.

The :doc:`Build basics <Build_basics>` doc page explains how to build
LAMMPS as either a shared or static library.  This results in one of
these 2 files:

liblammps.so      # shared library
liblammps.a       # static library


----------


**Link with LAMMPS as a static library**\ :

The calling application can link to LAMMPS as a static library with a
link command like this:

g++ caller.o -L/home/sjplimp/lammps/src -llammps -o caller

The -L argument is the path to where the liblammps.a file is.  The
-llammps argument is shorthand for the file liblammps.a.


----------


**Link with LAMMPS as a shared library**\ :

If you wish to link to liblammps.so, the operating system finds shared
libraries to load at run-time using the environment variable
LD\_LIBRARY\_PATH.  To enable this you can do one of two things:

(1) Copy the liblammps.so file to a location the system can find it,
such as /usr/local/lib.  I.e. a directory already listed in your
LD\_LIBRARY\_PATH variable.  You can type


.. parsed-literal::

   printenv LD_LIBRARY_PATH

to see what directories are in that list.

(2) Add the LAMMPS src directory (or the directory you perform CMake
build in) to your LD\_LIBRARY\_PATH, so that the current version of the
shared library is always available to programs that use it.

For the csh or tcsh shells, you would add something like this to your
~/.cshrc file:


.. parsed-literal::

   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home/sjplimp/lammps/src


----------


**Calling the LAMMPS library**\ :

Either flavor of library (static or shared) allows one or more LAMMPS
objects to be instantiated from the calling program.

When used from a C++ program, all of LAMMPS is wrapped in a LAMMPS\_NS
namespace; you can safely use any of its classes and methods from
within the calling code, as needed.

When used from a C or Fortran program, the library has a simple
C-style interface, provided in src/library.cpp and src/library.h.

See the :doc:`Python library <Python_library>` doc page for a
description of the Python interface to LAMMPS, which wraps the C-style
interface.

See the sample codes in examples/COUPLE/simple for examples of C++ and
C and Fortran codes that invoke LAMMPS through its library interface.
Other examples in the COUPLE directory use coupling ideas discussed on
the :doc:`Howto couple <Howto_couple>` doc page.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html

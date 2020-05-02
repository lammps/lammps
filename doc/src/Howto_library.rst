Library interface to LAMMPS
===========================

As described on the :doc:`Build basics <Build_basics>` doc page, LAMMPS
can be built as a static or shared library, so that it can be called by
another code, used in a :doc:`coupled manner <Howto_couple>` with other
codes, or driven through a :doc:`Python interface <Python_head>`.

At the core of LAMMPS is the ``LAMMPS`` class which encapsulates the
state of the simulation program through the state of the various class
instances that it is composed of.  So a calculation using LAMMPS
requires to create an instance of the ``LAMMPS`` class and then send it
(text) commands, either individually or from a file, or perform other
operations that modify the state stored inside that instance or drive
simulations.  This is essentially what the ``src/main.cpp`` file does
as well for the standalone LAMMPS executable with reading commands
either from an input file or stdin.

Creating a LAMMPS instance can be done by using C++ code directly or
through a C-style interface library to LAMMPS that is provided in the
files ``src/library.cpp`` and ``src/library.h``.  This
:ref:`C language API <lammps_c_api>`, can be used from C and C++,
and is also the basis for the :doc:`Python <pg_python>` and
:doc:`Fortran <pg_fortran>` interfaces or wrappers included in the
LAMMPS source code.

The ``examples/COUPLE`` and ``python/examples`` directories contain some
example programs written in C++, C, Fortran, and Python, which show how
a driver code can link to LAMMPS as a library, run LAMMPS on a subset of
processors (so the others are available to run some other code
concurrently), grab data from LAMMPS, change it, and send it back into
LAMMPS.

A detailed documentation of the available APIs and examples of how to
use them can be found in the :doc:`Programmer Documentation
<pg_library>` section of this manual.


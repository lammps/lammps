LAMMPS library C API
********************

Overview
========

As described on the :doc:`Build basics <Build_basics>` doc page, LAMMPS
can be built as a library, so that it can be called by another code,
used in a :doc:`coupled manner <Howto_couple>` with other codes, or
driven through a :doc:`Python interface <Python_head>`.

All of these methodologies use a C-style interface to LAMMPS that is
provided in the files ``src/library.cpp`` and ``src/library.h``.  The
functions therein have a C-style argument list, but contain C++ code you
could write yourself in a C++ application that was invoking LAMMPS
directly.  The C++ code in the functions illustrates how to invoke
internal LAMMPS operations.  Note that LAMMPS classes are defined
within a C++ namespace (LAMMPS_NS) if you use them from another C++
application.

The ``examples/COUPLE`` and ``python/examples`` directories have example
C++, Fortran, C, and Python codes which show how a driver code can link
to LAMMPS as a library, run LAMMPS on a subset of processors, grab data
from LAMMPS, change it, and put it back into LAMMPS.

Elsewhere in the manual you can find a description of the
:doc:`C++ base classes <pg_base>` and their public API.

Thread-safety
-------------

LAMMPS has not initially been conceived as a thread-safe program, but
over the years changes have been applied to replace operations that
collide with creating multiple LAMMPS instances from multiple-threads
of the same process with thread-safe alternatives.  This primarily
applies to the core LAMMPS code and less so on add-on packages, especially
when those packages require additional code in the *lib* folder,
interface LAMMPS to Fortran libraries, or the code uses static variables
(like the USER-COLVARS package.


Documented functions
====================

.. doxygenfile:: library.h
   :project: progguide

                 

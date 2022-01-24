LAMMPS Library Interfaces
*************************

As described on the :doc:`library interface to LAMMPS <Howto_library>`
page, LAMMPS can be built as a library (static or shared), so that
it can be called by another code, used in a :doc:`coupled manner
<Howto_couple>` with other codes, or driven through a :doc:`Python
script <Python_head>`.  Even the LAMMPS standalone executable is
essentially a thin wrapper on top of the LAMMPS library, creating a
LAMMPS instance, processing input and then existing.

Most of the APIs described below are based on C language wrapper
functions in the files ``src/library.h`` and ``src/library.cpp``, but
it is also possible to use C++ directly.  The basic procedure is
always the same: you create one or more instances of
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`, pass commands as strings or
from files to that LAMMPS instance to execute calculations, and/or
call functions that read, manipulate, and update data from the active
class instances inside LAMMPS to do analysis or perform operations
that are not possible with existing input script commands.

.. _thread-safety:

.. admonition:: Thread-safety
   :class: note

   LAMMPS was initially not conceived as a thread-safe program, but over
   the years changes have been applied to replace operations that
   collide with creating multiple LAMMPS instances from multiple-threads
   of the same process with thread-safe alternatives.  This primarily
   applies to the core LAMMPS code and less so on add-on packages,
   especially when those packages require additional code in the *lib*
   folder, interface LAMMPS to Fortran libraries, or the code uses
   static variables (like the COLVARS package).

   Another major issue to deal with is to correctly handle MPI.
   Creating a LAMMPS instance requires passing an MPI communicator, or
   it assumes the ``MPI_COMM_WORLD`` communicator, which spans all MPI
   processor ranks.  When creating multiple LAMMPS object instances
   from different threads, this communicator has to be different for
   each thread or else collisions can happen.  Or it has to be
   guaranteed, that only one thread at a time is active.  MPI
   communicators, however, are not a problem, if LAMMPS is compiled
   with the MPI STUBS library, which implies that there is no MPI
   communication and only 1 MPI rank.

----------

.. _lammps_c_api:

LAMMPS C Library API
====================

The C library interface is the most commonly used path to manage LAMMPS
instances from a compiled code and it is the basis for the :doc:`Python
<Python_module>` and :doc:`Fortran <Fortran>` modules.  Almost all
functions of the C language API require an argument containing a
"handle" in the form of a ``void *`` type variable, which points to the
location of a LAMMPS class instance.

The ``library.h`` header file by default does not include the ``mpi.h``
header file and thus hides the :cpp:func:`lammps_open` function which
requires the declaration of the ``MPI_comm`` data type.  This is only
a problem when the communicator that would be passed is different from
``MPI_COMM_WORLD``.  Otherwise calling :cpp:func:`lammps_open_no_mpi`
will work just as well.  To make :cpp:func:`lammps_open` available,
you need to compile the code with ``-DLAMMPS_LIB_MPI`` or add the line
``#define LAMMPS_LIB_MPI`` before ``#include "library.h"``.

Please note the ``mpi.h`` file must usually be the same (and thus the
MPI library in use) for the LAMMPS code and library and the calling code.
The exception is when LAMMPS was compiled in serial mode using the
``STUBS`` MPI library.  In that case the calling code may be compiled
with a different MPI library so long as :cpp:func:`lammps_open_no_mpi`
is called to create a LAMMPS instance.  In that case each MPI rank will
run LAMMPS in serial mode.

.. admonition:: Errors versus exceptions
   :class: note

   If the LAMMPS executable encounters an error condition, it will abort
   after printing an error message. For a library interface this is
   usually not desirable.  Thus LAMMPS can be compiled to to :ref:`throw
   a C++ exception <exceptions>` instead.  If enabled, the library
   functions will catch those exceptions and return.  The error status
   :cpp:func:`can be queried <lammps_has_error>` and an :cpp:func:`error
   message retrieved <lammps_get_last_error_message>`.  We thus
   recommend enabling C++ exceptions when using the library interface,

.. warning::

   No checks are made on the arguments of the function calls of the C
   library interface.  *All* function arguments must be non-NULL unless
   *explicitly* allowed, and must point to consistent and valid data.  Buffers
   for storing returned data must be allocated to a suitable size.
   Passing invalid or unsuitable information will likely cause crashes
   or corrupt data.

------------------------------

.. toctree::
   :maxdepth: 1

   Library_create
   Library_execute
   Library_properties
   Library_atoms
   Library_objects
   Library_scatter
   Library_neighbor
   Library_config
   Library_utility
   Library_add

--------------------

.. _lammps_python_api:

LAMMPS Python APIs
==================

The LAMMPS Python module enables calling the LAMMPS C library API from
Python by dynamically loading functions in the LAMMPS shared library through
the `Python ctypes module <https://docs.python.org/3/library/ctypes.html>`_.
Because of the dynamic loading, it is **required** that LAMMPS is compiled
in :ref:`"shared" mode <exe>`.  The Python interface is object oriented, but
otherwise tries to be very similar to the C library API.  Three different
Python classes to run LAMMPS are available and they build on each other.
More information on this is in the :doc:`Python_head`
section of the manual.  Use of the LAMMPS Python module is described in
:doc:`Python_module`.

-------------------

.. _lammps_fortran_api:

LAMMPS Fortran API
==================

The LAMMPS Fortran module is a wrapper around calling functions from the
LAMMPS C library API.  This is done using the ISO_C_BINDING feature in
Fortran 2003.  The interface is object oriented but otherwise tries to
be very similar to the C library API and the basic Python module.

.. toctree::
   :maxdepth: 1

   Fortran

-------------------

.. _lammps_cplusplus_api:

LAMMPS C++ API
==============

It is also possible to invoke the LAMMPS C++ API directly in your code.
It lacks some of the convenience of the C library API, but it allows
more direct access to simulation data and thus more low-level manipulations.
The following links provide some examples and references to the C++ API.

.. toctree::
   :maxdepth: 1

   Cplusplus



LAMMPS Library Interfaces
*************************

As described on the :doc:`library interface to LAMMPS <Howto_library>`
doc page, LAMMPS can be built as a library (static or shared), so that
it can be called by another code, used in a :doc:`coupled manner
<Howto_couple>` with other codes, or driven through a :doc:`Python
script <Python_head>`.  Even the LAMMPS standalone executable is
essentially a thin wrapper on top of the LAMMPS library, creating a
LAMMPS instance, processing input and then existing.

Several of these approaches are based on C language wrapper functions
in the files ``src/library.h`` and ``src/library.cpp``, but it is also
possible to use C++ directly.  The basic procedure is always the same:
you create one or more instances of the
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` and then pass commands as
strings or from files to that LAMMPS instance to execute calculations,
or read, manipulate, and update data from the active class instances
inside the LAMMPS to do analysis or perform operations that are not
possible with existing commands.

.. _lammps_c_api:

LAMMPS C Library API
====================

The C library interface is most commonly used path to manage LAMMPS
instances from a compiled code and it is the bases for the
:doc:`Python <pg_python>` and :doc:`Fortran <pg_fortran>` modules.
Almost all functions of the C language API require an argument containing
a "handle" in the form of a ``void *`` type variable, which points to the
location of a LAMMPS class instance.

The ``library.h`` header file by default includes the ``mpi.h`` header
for an MPI library, so it must be present when compiling code using the
library interface.  This usually must be the header from the same MPI
library as the LAMMPS library was compiled with.  The exception is when
LAMMPS was compiled in serial mode using the ``STUBS`` MPI library.  In
that case the calling code may be compiled with a different MPI library
for as long as :cpp:func:`lammps_open_no_mpi` is called to create a
LAMMPS instance. Then you may set the define ``-DLAMMPS_LIB_NO_MPI``
when compiling your code and the inclusion of ``mpi.h`` will be skipped
and consequently the function :cpp:func:`lammps_open` may not be used.

If any of the function calls in the LAMMPS library API will trigger
an error inside LAMMPS, this will result in an abort of the entire
program.  This is not always desirable.  Instead, LAMMPS can be
compiled to instead :ref:`throw a C++ exception <exceptions>`.

.. _thread-safety:

.. note:: Thread-safety

   LAMMPS was initially not conceived as a thread-safe program, but over
   the years changes have been applied to replace operations that
   collide with creating multiple LAMMPS instances from multiple-threads
   of the same process with thread-safe alternatives.  This primarily
   applies to the core LAMMPS code and less so on add-on packages,
   especially when those packages require additional code in the *lib*
   folder, interface LAMMPS to Fortran libraries, or the code uses
   static variables (like the USER-COLVARS package).

   Another major issue to deal with is to correctly handle MPI.
   Creating a LAMMPS instance requires passing an MPI communicator, or
   it assumes the ``MPI_COMM_WORLD`` communicator, which spans all MPI
   processor ranks.  When creating multiple LAMMPS object instances from
   different threads, this communicator has to be different for each
   thread or else collisions can happen.  or it has to be guaranteed,
   that only one thread at a time is active.  MPI communicators,
   however, are not a problem, if LAMMPS is compiled with the MPI STUBS
   library, which implies that there is no MPI communication and only 1
   MPI rank.

.. warning::
   No checks are made on the arguments of the function calls of the C
   library interface.  All arguments must be non-NULL unless explicitly
   allowed and point to consistent and valid data.  Buffers for storing
   returned data must be allocated to a suitable size.  Passing invalid
   or unsuitable information will likely cause crashes or corrupt data.

------------------------------

.. toctree::
   :maxdepth: 1

   pg_lib_create
   pg_lib_execute
   pg_lib_properties
   pg_lib_objects
   pg_lib_neighbor
   pg_lib_config
   pg_lib_utility

------------------- 

.. _lammps_python_api:

LAMMPS Python Module API
========================

.. toctree::
   :maxdepth: 1

   pg_python

-------------------

.. _lammps_fortran_api:

LAMMPS Fortran Module API
=========================

.. toctree::
   :maxdepth: 1

   pg_fortran

-------------------

TODO: this part still needs to be edited/adapted

.. note::

   You can write code for additional functions as needed to define
   how your code talks to LAMMPS and add them to src/library.cpp and
   src/library.h, as well as to the :doc:`Python interface <Python_head>`.
   The added functions can access or change any internal LAMMPS data you
   wish.


.. code-block:: c

   void lammps_gather_atoms(void *, char *, int, int, void *)
   void lammps_gather_atoms_concat(void *, char *, int, int, void *)
   void lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *)
   void lammps_scatter_atoms(void *, char *, int, int, void *)
   void lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *)

The gather functions collect peratom info of the requested type (atom
coords, atom types, forces, etc) from all processors, and returns the
same vector of values to each calling processor.  The scatter
functions do the inverse.  They distribute a vector of peratom values,
passed by all calling processors, to individual atoms, which may be
owned by different processors.

.. warning::

   These functions are not compatible with the
   -DLAMMPS_BIGBIG setting when compiling LAMMPS.  Dummy functions
   that result in an error message and abort will be substituted
   instead of resulting in random crashes and memory corruption.

The lammps_gather_atoms() function does this for all N atoms in the
system, ordered by atom ID, from 1 to N.  The
lammps_gather_atoms_concat() function does it for all N atoms, but
simply concatenates the subset of atoms owned by each processor.  The
resulting vector is not ordered by atom ID.  Atom IDs can be requested
by the same function if the caller needs to know the ordering.  The
lammps_gather_subset() function allows the caller to request values
for only a subset of atoms (identified by ID).
For all 3 gather function, per-atom image flags can be retrieved in 2 ways.
If the count is specified as 1, they are returned
in a packed format with all three image flags stored in a single integer.
If the count is specified as 3, the values are unpacked into xyz flags
by the library before returning them.

The lammps_scatter_atoms() function takes a list of values for all N
atoms in the system, ordered by atom ID, from 1 to N, and assigns
those values to each atom in the system.  The
lammps_scatter_atoms_subset() function takes a subset of IDs as an
argument and only scatters those values to the owning atoms.

.. removed from Build_link.rst

See the sample codes in examples/COUPLE/simple for examples of C++ and
C and Fortran codes that invoke LAMMPS through its library interface.
Other examples in the COUPLE directory use coupling ideas discussed on
the :doc:`Howto couple <Howto_couple>` doc page.

.. _ctypes: https://docs.python.org/3/library/ctypes.html

.. removed from Howto_couple.rst

Examples of driver codes that call LAMMPS as a library are included in
the examples/COUPLE directory of the LAMMPS distribution; see
examples/COUPLE/README for more details:

* simple: simple driver programs in C++ and C which invoke LAMMPS as a
  library
* plugin: simple driver program in C which invokes LAMMPS as a plugin
  from a shared library.
* lammps_quest: coupling of LAMMPS and `Quest <quest_>`_, to run classical
  MD with quantum forces calculated by a density functional code
* lammps_spparks: coupling of LAMMPS and `SPPARKS <spparks_>`_, to couple
  a kinetic Monte Carlo model for grain growth using MD to calculate
  strain induced across grain boundaries

.. _quest: http://dft.sandia.gov/Quest

.. _spparks: http://www.sandia.gov/~sjplimp/spparks.html

The :doc:`Build basics <Build_basics>` doc page describes how to build
LAMMPS as a library.  Once this is done, you can interface with LAMMPS
either via C++, C, Fortran, or Python (or any other language that
supports a vanilla C-like interface).  For example, from C++ you could
create one (or more) "instances" of LAMMPS, pass it an input script to
process, or execute individual commands, all by invoking the correct
class methods in LAMMPS.  From C or Fortran you can make function
calls to do the same things.  See the :doc:`Python <Python_head>` doc
pages for a description of the Python wrapper provided with LAMMPS
that operates through the LAMMPS library interface.

The files src/library.cpp and library.h contain the C-style interface
to LAMMPS.  See the :doc:`Howto library <Howto_library>` doc page for a
description of the interface and how to extend it for your needs.

Note that the lammps_open() function that creates an instance of
LAMMPS takes an MPI communicator as an argument.  This means that
instance of LAMMPS will run on the set of processors in the
communicator.  Thus the calling code can run LAMMPS on all or a subset
of processors.  For example, a wrapper script might decide to
alternate between LAMMPS and another code, allowing them both to run
on all the processors.  Or it might allocate half the processors to
LAMMPS and half to the other code and run both codes simultaneously
before syncing them up periodically.  Or it might instantiate multiple
instances of LAMMPS to perform different calculations.





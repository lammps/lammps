LAMMPS C library API
********************

Overview
========

As described on the :doc:`library interface to LAMMPS <Howto_library>`
doc page, LAMMPS can be built as a library, so that it can be called by
another code, used in a :doc:`coupled manner <Howto_couple>` with other
codes, or driven through a :doc:`Python interface <Python_head>`.
Several of these approaches are based on a C language wrapper functions
in the files ``src/library.h`` and ``src/library.cpp``, which are
documented below.

Behind the scenes this will create, delete, or modify an instance of the
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class.  Thus almost all functions
require an argument containing a "handle" in the form of a ``void *``
type variable, which points to the location this LAMMPS class instance.


The ``library.h`` header file includes the ``mpi.h`` header for an MPI
library, so it must be present when compiling code using the library
interface.  This usually must be the header from the same MPI library as
the LAMMPS library was compiled with.  The exception is when LAMMPS was
compiled in serial mode using the ``STUBS`` MPI library.  In that case
the calling code may be compiled with a different MPI library for as
long as :cpp:func:`lammps_open_no_mpi` is called to
create a LAMMPS instance.

If any of the function calls in the LAMMPS library API will trigger
an error inside LAMMPS, this will result in an abort of the entire
program.  This is not always desirable.  Instead, LAMMPS can be
compiled to instead :ref:`throw a C++ exception <exceptions>`. 

.. seealso::

   Please see the :ref:`note about thread-safety <thread-safety>`
   in the library Howto doc page.

.. warning::
   No checks are made on the arguments of the function calls of the C
   library interface.  All arguments must be non-NULL unless explicitly
   allowed and point to consistent and valid data.  Buffers for storing
   returned data must be allocated to a suitable size.  Passing invalid
   or unsuitable information will likely cause crashes or corrupt data.


Creating or deleting a LAMMPS object
====================================

The :cpp:func:`lammps_open` and :cpp:func:`lammps_open_no_mpi`
functions are used to create and initialize a
:cpp:func:`LAMMPS` instance.  The calling program has to
provide a handle where a reference to this instance can be stored and
which has to be used in all subsequent function calls until that
instance is destroyed by calling :cpp:func:`lammps_close`.
Here is a simple example demonstrating its use:

.. code-block:: C

   #include "library.h"
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;
     int version;
     const char *lmpargv[] { "liblammps", "-log", "none"};
     int lmpargc = sizeof(lmpargv)/sizeof(const char *);

     /* create LAMMPS instance */
     lammps_open_no_mpi(lmpargc, lmpargv, &handle);
     if (handle == NULL) {
       printf("LAMMPS initialization failed");
       lammps_finalize();
       return 1;
     }

     /* get and print numerical version */
     version = lammps_version(handle);
     printf("LAMMPS Version: %d\n",version);

     /* delete LAMMPS instance and shut down MPI */
     lammps_close(handle);
     lammps_finalize();
     return 0;
   }

The LAMMPS library will be using the MPI library it was compiled with
and will either run on all processors in the ``MPI_COMM_WORLD``
communicator or on the set of processors in the communicator given in
the ``comm`` argument of :cpp:func:`lammps_open`.  This means
the calling code can run LAMMPS on all or a subset of processors.  For
example, a wrapper code might decide to alternate between LAMMPS and
another code, allowing them both to run on all the processors.  Or it
might allocate part of the processors to LAMMPS and the rest to the
other code by creating a custom communicator with ``MPI_Comm_split()``
and running both codes concurrently before syncing them up periodically.
Or it might instantiate multiple instances of LAMMPS to perform
different calculations and either alternate between them, run them
concurrently on split communicators, or run them one after the other.
The :cpp:func:`lammps_open` function may be called multiple
times for this latter purpose.

The :cpp:func:`lammps_close` function is used to shut down
the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class pointed to by the handle
passed as an argument and free all its memory. This has to be called for
every instance created with any of the :cpp:func:`lammps_open` functions.  It will, however, **not** call
``MPI_Finalize()``, since that may only be called once.  See
:cpp:func:`lammps_finalize` for an alternative to calling
``MPI_Finalize()`` explicitly in the calling program.

The :cpp:func:`lammps_free` function is a clean-up
function to free memory that the library allocated previously
via other function calls.  See below for notes in the descriptions
of the individual commands where such memory buffers were allocated.

-----------------------

.. doxygenfunction:: lammps_open
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_open_no_mpi
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_close
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_finalize
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_free
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_version
   :project: progguide

-----------------------

Executing LAMMPS commands
=========================

Once a LAMMPS instance is created, there are multiple ways to "drive" a
simulation.  In most cases it is easiest to process single or multiple
LAMMPS commands like in an input file.  This can be done through reading
a file or passing single commands or lists of commands or blocks of
commands with the following functions.

Via these functions, the calling code can have the LAMMPS instance act
on a series of :doc:`input file commands <Commands_all>` that are either
read from a file or passed as strings.  This for, for example, allows to
setup a problem from a template file and then run it in stages while
performing other operations in between or concurrently.  The caller can
interleave the LAMMPS function calls with operations it performs, calls
to extract information from or set information within LAMMPS, or calls
to another code's library.

Also equivalent to regular :doc:`input script parsing <Commands_parse>`
is the handling of comments and expansion of variables with ``${name}``
or ``$(expression)`` syntax before the commands are parsed and
executed. Below is a short example using some of these functions.

.. code-block:: C

   #include <library.h>
   #include <mpi.h>
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;
     int i;

     MPI_Init(&argc, &argv);
     lammps_open(0, NULL, MPI_COMM_WORLD, &handle);
     lammps_file(handle,"in.sysinit");
     lammps_command(handle,"run 1000 post no");

     for (i=0; i < 100; ++i) {
       lammps_commands_string(handle,"run 100 pre no post no\n"
                                     "print 'PE = $(pe)'\n"
                                     "print 'KE = $(ke)'\n")
     }
     lammps_close(handle);
     MPI_Finalize();
     return 0;
   }

-----------------------

.. doxygenfunction:: lammps_file
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_command
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_commands_list
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_commands_string
   :project: progguide

-----------------------

Retrieving or setting system properties
=======================================

The library interface allows to extract all kinds of information
about the active simulation instance and also modify it.  This
allows to combine MD simulation steps with other processing and
simulation methods computed in the calling code or another code
that is coupled to LAMMPS via the library interface.

-----------------------

.. doxygenfunction:: lammps_extract_setting
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_box
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_extract_atom
   :project: progguide

-------------------

.. doxygenfunction:: lammps_get_natoms
   :project: progguide

-------------------

TODO: this part still needs to be edited/adapted

.. note::

   You can write code for additional functions as needed to define
   how your code talks to LAMMPS and add them to src/library.cpp and
   src/library.h, as well as to the :doc:`Python interface <Python_head>`.
   The added functions can access or change any internal LAMMPS data you
   wish.


The file src/library.cpp also contains these functions for extracting
information from LAMMPS and setting value within LAMMPS.  Again, see
the documentation in the src/library.cpp file for details, including
which quantities can be queried by name:

.. code-block:: c

   int lammps_extract_setting(void *, char *)
   void *lammps_extract_global(void *, char *)
   void lammps_extract_box(void *, double *, double *,
                           double *, double *, double *, int *, int *)
   void *lammps_extract_atom(void *, char *)
   void *lammps_extract_compute(void *, char *, int, int)
   void *lammps_extract_fix(void *, char *, int, int, int, int)
   void *lammps_extract_variable(void *, char *, char *)

The extract_setting() function returns info on the size
of data types (e.g. 32-bit or 64-bit atom IDs) used
by the LAMMPS executable (a compile-time choice).

The other extract functions return a pointer to various global or
per-atom quantities stored in LAMMPS or to values calculated by a
compute, fix, or variable.  The pointer returned by the
extract_global() function can be used as a permanent reference to a
value which may change.  For the extract_atom() method, see the
extract() method in the src/atom.cpp file for a list of valid per-atom
properties.  New names could easily be added if the property you want
is not listed.  For the other extract functions, the underlying
storage may be reallocated as LAMMPS runs, so you need to re-call the
function to assure a current pointer or returned value(s).

.. code-block:: c

   double lammps_get_thermo(void *, char *)
   int lammps_get_natoms(void *)

   int lammps_set_variable(void *, char *, char *)
   void lammps_reset_box(void *, double *, double *, double, double, double)

The lammps_get_thermo() function returns the current value of a thermo
keyword as a double precision value.

The lammps_get_natoms() function returns the total number of atoms in
the system and can be used by the caller to allocate memory for the
lammps_gather_atoms() and lammps_scatter_atoms() functions.

The lammps_set_variable() function can set an existing string-style
variable to a new string value, so that subsequent LAMMPS commands can
access the variable.

The lammps_reset_box() function resets the size and shape of the
simulation box, e.g. as part of restoring a previously extracted and
saved state of a simulation.

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

.. code-block:: c

   void lammps_create_atoms(void *, int, tagint *, int *, double *, double *,
                            imageint *, int)

The lammps_create_atoms() function takes a list of N atoms as input
with atom types and coords (required), an optionally atom IDs and
velocities and image flags.  It uses the coords of each atom to assign
it as a new atom to the processor that owns it.  This function is
useful to add atoms to a simulation or (in tandem with
lammps_reset_box()) to restore a previously extracted and saved state
of a simulation.  Additional properties for the new atoms can then be
assigned via the lammps_scatter_atoms() or lammps_extract_atom()
functions.

.. removed from Build_link.rst

**Calling the LAMMPS library**\ :

Either flavor of library (static or shared) allows one or more LAMMPS
objects to be instantiated from the calling program. When used from a
C++ program, most of the symbols and functions in LAMMPS are wrapped
in a LAMMPS_NS namespace; you can safely use any of its classes and
methods from within the calling code, as needed, and you will not incur
conflicts with functions and variables in your code that share the name.
This, however, does not extend to all additional libraries bundled with
LAMMPS in the lib folder and some of the low-level code of some packages.

To be compatible with C, Fortran, Python programs, the library has a simple
C-style interface, provided in src/library.cpp and src/library.h.

See the :doc:`Python library <Python_library>` doc page for a
description of the Python interface to LAMMPS, which wraps the C-style
interface from a shared library through the `ctypes python module <ctypes_>`_.

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
   


                 

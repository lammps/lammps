Creating or deleting a LAMMPS object
====================================

This section documents the following functions:

- :cpp:func:`lammps_open`
- :cpp:func:`lammps_open_no_mpi`
- :cpp:func:`lammps_open_fortran`
- :cpp:func:`lammps_close`
- :cpp:func:`lammps_mpi_init`
- :cpp:func:`lammps_mpi_finalize`
- :cpp:func:`lammps_kokkos_finalize`
- :cpp:func:`lammps_python_finalize`
- :cpp:func:`lammps_error`

--------------------

The :cpp:func:`lammps_open` and :cpp:func:`lammps_open_no_mpi` functions
are used to create and initialize a :cpp:func:`LAMMPS` instance.  They
return a reference to this instance as a ``void *`` pointer to be used
as the "handle" argument in subsequent function calls until that
instance is destroyed by calling :cpp:func:`lammps_close`.  Here is a
simple example demonstrating its use:

.. code-block:: c

   #include "library.h"
   #include <stdio.h>

   int main(int argc, char **argv)
   {
     void *handle;
     int version;
     const char *lmpargv[] = { "liblammps", "-log", "none"};
     int lmpargc = sizeof(lmpargv)/sizeof(const char *);

     /* create LAMMPS instance */
     handle = lammps_open_no_mpi(lmpargc, (char **)lmpargv, NULL);
     if (handle == NULL) {
       printf("LAMMPS initialization failed");
       lammps_mpi_finalize();
       return 1;
     }

     /* get and print numerical version code */
     version = lammps_version(handle);
     printf("LAMMPS Version: %d\n",version);

     /* delete LAMMPS instance and shut down MPI */
     lammps_close(handle);
     lammps_mpi_finalize();
     return 0;
   }

The LAMMPS library uses the MPI library it was compiled with and will
either run on all processors in the ``MPI_COMM_WORLD`` communicator or
on the set of processors in the communicator passed as the ``comm``
argument of :cpp:func:`lammps_open`.  This means the calling code can
run LAMMPS on all or a subset of processors.  For example, a wrapper
code might decide to alternate between LAMMPS and another code, allowing
them both to run on all the processors.  Or it might allocate part of
the processors to LAMMPS and the rest to the other code by creating a
custom communicator with ``MPI_Comm_split()`` and running both codes
concurrently before syncing them up periodically.  Or it might
instantiate multiple instances of LAMMPS to perform different
calculations and either alternate between them, run them concurrently on
split communicators, or run them one after the other.  The
:cpp:func:`lammps_open` function may be called multiple times for this
latter purpose.

The :cpp:func:`lammps_close` function is used to shut down the
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class pointed to by the handle
passed as an argument and free all its memory. This has to be called
for every instance created with one of the :cpp:func:`lammps_open`
functions.  It will, however, **not** call ``MPI_Finalize()``, since
that may only be called once.  See :cpp:func:`lammps_mpi_finalize` for
an alternative to invoking ``MPI_Finalize()`` explicitly from the
calling program.

-----------------------

.. doxygenfunction:: lammps_open
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_open_no_mpi
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_open_fortran
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_close
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_mpi_init
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_mpi_finalize
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_kokkos_finalize
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_python_finalize
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_error
   :project: progguide

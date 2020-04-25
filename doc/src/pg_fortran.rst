LAMMPS Fortran Module
*********************

The LAMMPS module provides an interface to call LAMMPS from a Fortran code.
It is based on the C-library interface and similar to the LAMMPS Python
wrapper provides an object oriented interface. Below is a minimal example:

.. code-block:: fortran

   PROGRAM testlib
     USE LIBLAMMPS
     TYPE(lammps)     :: lmp
     INTEGER          :: argc,ierr
     CHARACTER(len=200),ALLOCATABLE :: argv(:)

     argc = 3
     ALLOCATE(argv(argc))
     argv(1) = 'liblammps'
     argv(2) = '-log'
     argv(3) = 'log.fortran'

     ! create a LAMMPS instance
     lmp = lammps(argc,argv)
     call lmp%file('in.melt')
     call lmp%command('run 200')
     call lmp%close(.true.)
     deallocate(argv)

   END PROGRAM testlib

--------------------

.. f:type:: lammps

   Derived type that is the general class of the Fortran interface.
   It holds a reference to the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class instance
   that any of the included calls are forwarded to.

   :f c_ptr handle: reference to the LAMMPS class
   :f close: :f:func:`~LIBLAMMPS/close`

.. f:function:: lammps(argc,argv[,comm])

   This is the constructor for the Fortran class and will forward
   the arguments to a call to :cpp:func:`lammps_open_fortran` or :cpp:func:`lammps_open_no_mpi`.

   :p integer argc: number of arguments
   :p character(len=*) argv(): arguments as list of strings
   :o integer comm [optional]: MPI communicator
   :r lammps: an instance of the :f:type:`lammps` derived type

.. f:subroutine:: close([finalize])

   This method will close down the LAMMPS instance through calling :cpp:func:`lammps_close`.
   If the *finalize* argument is present and has a value of ``.true.``, then this subroutine
   also called ``MPI_Finalize()``.

   :o logical finalize [optional]: call ``MPI_Finalize()`` after deleting the LAMMPS instance

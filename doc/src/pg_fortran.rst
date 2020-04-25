LAMMPS Fortran Module
*********************

The LAMMPS module provides an interface to call LAMMPS from a Fortran code.
It is based on the C-library interface and requires a Fortran 2003 compatible
compiler to be compiled.  Similar to the LAMMPS Python wrapper it provides an
object oriented interface.  Below is a minimal example:

.. code-block:: fortran

   PROGRAM testlib
     USE LIBLAMMPS                 ! include the LAMMPS library interface
     TYPE(lammps)     :: lmp       ! derived type to hold LAMMPS instance
     CHARACTER(len=*), DIMENSION(5), PARAMETER :: args = &
         [ CHARACTER(len=12) :: 'liblammps',&
          '-echo', 'both', '-log', 'log.fortran' ]

     ! create a LAMMPS instance (and initialize MPI)
     lmp = lammps(args)
     ! read commands from a file
     CALL lmp%file('in.melt')
     ! execute a single command
     CALL lmp%command('run 100 post no')
     ! delete LAMMPS instance (and shut down MPI)
     CALL lmp%close(.true.)

   END PROGRAM testlib

To compile and link Fortran code with the LAMMPS library, you need to compile
the Fortran library module (located in ``examples/COUPLE/fortran/lammps.f90``)
and link it with your code and also :doc:`link to the LAMMPS library <Build_link>`.
A typical command line would be:

.. code-block:: bash

   gfortran -o testlib.x  lammps.f90 testlib.f90 -L. -llammps
 
--------------------

.. f:type:: lammps

   Derived type that is the general class of the Fortran interface.
   It holds a reference to the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class instance
   that any of the included calls are forwarded to.

   :f c_ptr handle: reference to the LAMMPS class
   :f close: :f:func:`close`
   :f file: :f:func:`file`
   :f command: :f:func:`command`

.. f:function:: lammps(args[,comm])

   This is the constructor for the Fortran class and will forward
   the arguments to a call to either :cpp:func:`lammps_open_fortran`
   or :cpp:func:`lammps_open_no_mpi`. If the LAMMPS library has been
   compiled with MPI support, it will also initialize MPI, if it has
   not already been initialized before.

   The *comm* argument with the MPI communicator is optional. If it
   is not provided, ``MPI_COMM_WORLD`` is assumed. For more details
   please see the documentation of :cpp:func:`lammps_open`.

   :p character(len=*) args(): arguments as list of strings
   :o integer comm [optional]: MPI communicator
   :r lammps: an instance of the :f:type:`lammps` derived type

.. f:subroutine:: close([finalize])

   This method will close down the LAMMPS instance through calling
   :cpp:func:`lammps_close`.  If the *finalize* argument is present and
   has a value of ``.true.``, then this subroutine also calls
   :cpp:func:`lammps_finalize`.

   :o logical finalize [optional]: shut down the MPI environment of the LAMMPS library if true.

.. f:subroutine:: file(filename)

   This method will call :cpp:func:`lammps_file` to have LAMMPS read
   and process commands from a file.

   :p character(len=*) filename: name of file with LAMMPS commands

.. f:subroutine:: command(cmd)

   This method will call :cpp:func:`lammps_command` to have LAMMPS
   execute a single command.

   :p character(len=*) cmd: LAMMPS command

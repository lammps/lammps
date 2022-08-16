The ``LIBLAMMPS`` Fortran Module
********************************

The ``LIBLAMMPS`` module provides an interface to call LAMMPS from a
Fortran code.  It is based on the LAMMPS C-library interface and
requires a Fortran 2003 compatible compiler to be compiled.

While C libraries have a defined binary interface (ABI) and can thus be
used from multiple compiler versions from different vendors for as long
as they are compatible with the hosting operating system, the same is
not true for Fortran codes.  Thus the LAMMPS Fortran module needs to be
compiled alongside the code using it from the source code in
``fortran/lammps.f90``.  When linking, you also need to
:doc:`link to the LAMMPS library <Build_link>`.  A typical command line
for a simple program using the Fortran interface would be:

.. code-block:: bash

   mpifort -o testlib.x  lammps.f90 testlib.f90 -L. -llammps

Please note, that the MPI compiler wrapper is only required when the
calling the library from an MPI parallel code.  Please also note the
order of the source files: the ``lammps.f90`` file needs to be compiled
first, since it provides the ``LIBLAMMPS`` module that is imported by
the Fortran code using the interface.  A working example code can be
found together with equivalent examples in C and C++ in the
``examples/COUPLE/simple`` folder of the LAMMPS distribution.

.. versionadded:: 9Oct2020

.. admonition:: Work in Progress
   :class: note

   This Fortran module is work in progress and only the documented
   functionality is currently available. The final implementation should
   cover the entire range of functionality available in the C and
   Python library interfaces.

.. note::

   A contributed (and more complete!) Fortran interface that more
   closely resembles the C-library interface is available in the
   ``examples/COUPLE/fortran2`` folder.  Please see the ``README`` file
   in that folder for more information about it and how to contact its
   author and maintainer.

----------

Creating or deleting a LAMMPS object
************************************

With the Fortran interface the creation of a :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` instance is included in the constructor for
creating the :f:func:`lammps` derived type.  To import the definition of
that type and its type bound procedures you need to add a ``USE
LIBLAMMPS`` statement.  Internally it will call either
:cpp:func:`lammps_open_fortran` or :cpp:func:`lammps_open_no_mpi` from
the C library API to create the class instance.  All arguments are
optional and :cpp:func:`lammps_mpi_init` will be called automatically,
if it is needed.  Similarly, a possible call to :cpp:func:`lammps_finalize`
is integrated into the :f:func:`close` function and triggered with
the optional logical argument set to ``.true.``. Here is a simple example:

.. code-block:: fortran

   PROGRAM testlib
     USE LIBLAMMPS                 ! include the LAMMPS library interface
     IMPLICIT NONE
     TYPE(lammps)     :: lmp       ! derived type to hold LAMMPS instance
     CHARACTER(len=*), PARAMETER :: args(3) = &
         [ CHARACTER(len=12) :: 'liblammps', '-log', 'none' ]

     ! create a LAMMPS instance (and initialize MPI)
     lmp = lammps(args)
     ! get and print numerical version code
     PRINT*, 'LAMMPS Version: ', lmp%version()
     ! delete LAMMPS instance (and shuts down MPI)
     CALL lmp%close(.true.)

   END PROGRAM testlib

It is also possible to pass command line flags from Fortran to C/C++ and
thus make the resulting executable behave similar to the standalone
executable (it will ignore the `-in/-i` flag, though).  This allows to
use the command line to configure accelerator and suffix settings,
configure screen and logfile output, or to set index style variables
from the command line and more. Here is a correspondingly adapted
version of the previous example:

.. code-block:: fortran

   PROGRAM testlib2
     USE LIBLAMMPS                 ! include the LAMMPS library interface
     IMPLICIT NONE
     TYPE(lammps)     :: lmp       ! derived type to hold LAMMPS instance
     CHARACTER(len=128), ALLOCATABLE :: command_args(:)
     INTEGER :: i, argc

     ! copy command line flags to `command_args()`
     argc = COMMAND_ARGUMENT_COUNT()
     ALLOCATE(command_args(0:argc))
     DO i=0, argc
       CALL GET_COMMAND_ARGUMENT(i, command_args(i))
     END DO

     ! create a LAMMPS instance (and initialize MPI)
     lmp = lammps(command_args)
     ! get and print numerical version code
     PRINT*, 'Program name:   ', command_args(0)
     PRINT*, 'LAMMPS Version: ', lmp%version()
     ! delete LAMMPS instance (and shuts down MPI)
     CALL lmp%close(.TRUE.)
     DEALLOCATE(command_args)

   END PROGRAM testlib2

--------------------

Executing LAMMPS commands
=========================

Once a LAMMPS instance is created, it is possible to "drive" the LAMMPS
simulation by telling LAMMPS to read commands from a file, or pass
individual or multiple commands from strings or lists of strings.  This
is done similar to how it is implemented in the `C-library
<pg_lib_execute>` interface. Before handing off the calls to the
C-library interface, the corresponding Fortran versions of the calls
(:f:func:`file`, :f:func:`command`, :f:func:`commands_list`, and
:f:func:`commands_string`) have to make a copy of the strings passed as
arguments so that they can be modified to be compatible with the
requirements of strings in C without affecting the original strings.
Those copies are automatically deleted after the functions return.
Below is a small demonstration of the uses of the different functions:

.. code-block:: fortran

   PROGRAM testcmd
     USE LIBLAMMPS
     TYPE(lammps)     :: lmp
     CHARACTER(len=512) :: cmds
     CHARACTER(len=40), ALLOCATABLE :: cmdlist(:)
     CHARACTER(len=10) :: trimmed
     INTEGER :: i

     lmp = lammps()
     CALL lmp%file('in.melt')
     CALL lmp%command('variable zpos index 1.0')
     ! define 10 groups of 10 atoms each
     ALLOCATE(cmdlist(10))
     DO i=1, 10
         WRITE(trimmed,'(I10)') 10*i
         WRITE(cmdlist(i),'(A,I1,A,I10,A,A)')       &
             'group g', i-1, ' id ', 10*(i-1)+1, ':', ADJUSTL(trimmed)
     END DO
     CALL lmp%commands_list(cmdlist)
     ! run multiple commands from multi-line string
     cmds = 'clear' // NEW_LINE('A') //                       &
         'region  box block 0 2 0 2 0 2' // NEW_LINE('A') //  &
         'create_box 1 box' // NEW_LINE('A') //               &
         'create_atoms 1 single 1.0 1.0 ${zpos}'
     CALL lmp%commands_string(cmds)
     CALL lmp%close(.TRUE.)

   END PROGRAM testcmd

---------------

The ``LIBLAMMPS`` module API
****************************

Below are the detailed descriptions of definitions and interfaces
of the contents of the ``LIBLAMMPS`` Fortran interface to LAMMPS.

.. f:type:: lammps

   Derived type that is the general class of the Fortran interface.  It
   holds a reference to the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`
   class instance that any of the included calls are forwarded to.

   :f c_ptr handle: reference to the LAMMPS class
   :f close: :f:func:`close`
   :f version: :f:func:`version`
   :f file: :f:func:`file`
   :f command: :f:func:`command`
   :f commands_list: :f:func:`commands_list`
   :f commands_string: :f:func:`commands_string`

.. f:function:: lammps(args[,comm])

   This is the constructor for the Fortran class and will forward
   the arguments to a call to either :cpp:func:`lammps_open_fortran`
   or :cpp:func:`lammps_open_no_mpi`. If the LAMMPS library has been
   compiled with MPI support, it will also initialize MPI, if it has
   not already been initialized before.

   The *args* argument with the list of command line parameters is
   optional and so it the *comm* argument with the MPI communicator.
   If *comm* is not provided, ``MPI_COMM_WORLD`` is assumed. For
   more details please see the documentation of :cpp:func:`lammps_open`.

   :p character(len=*) args(*) [optional]: arguments as list of strings
   :o integer comm [optional]: MPI communicator
   :r lammps: an instance of the :f:type:`lammps` derived type

.. f:subroutine:: close([finalize])

   This method will close down the LAMMPS instance through calling
   :cpp:func:`lammps_close`.  If the *finalize* argument is present and
   has a value of ``.true.``, then this subroutine also calls
   :cpp:func:`lammps_mpi_finalize`.

   :o logical finalize [optional]: shut down the MPI environment of the LAMMPS library if true.

.. f:function:: version()

   This method returns the numeric LAMMPS version like :cpp:func:`lammps_version`

   :r integer: LAMMPS version

--------

.. f:subroutine:: file(filename)

   This method will call :cpp:func:`lammps_file` to have LAMMPS read
   and process commands from a file.

   :p character(len=*) filename: name of file with LAMMPS commands

.. f:subroutine:: command(cmd)

   This method will call :cpp:func:`lammps_command` to have LAMMPS
   execute a single command.

   :p character(len=*) cmd: single LAMMPS command

.. f:subroutine:: commands_list(cmds)

   This method will call :cpp:func:`lammps_commands_list` to have LAMMPS
   execute a list of input lines.

   :p character(len=*) cmd(:): list of LAMMPS input lines

.. f:subroutine:: commands_string(str)

   This method will call :cpp:func:`lammps_commands_string` to have LAMMPS
   execute a block of commands from a string.

   :p character(len=*) str: LAMMPS input in string

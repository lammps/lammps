The ``LIBLAMMPS`` Fortran Module
********************************

The ``LIBLAMMPS`` module provides an interface to call LAMMPS from a
Fortran code.  It is based on the LAMMPS C-library interface and
requires a Fortran 2003 compatible compiler to be compiled.  It is
designed to be self-contained and not require any support functions
written in C, C++, or Fortran.

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
calling the library from an MPI parallel code.  Otherwise, using the
fortran compiler (gfortran, ifort, flang, etc.) will suffice.  It may be
necessary to link to additional libraries depending on how LAMMPS was
configured and whether the LAMMPS library :doc:`was compiled as a static
or shared library <Build_link>`.

If the LAMMPS library itself has been compiled with MPI support, the
resulting executable will still be able to run LAMMPS in parallel with
``mpirun`` or equivalent.  Please also note that the order of the source
files matters: the ``lammps.f90`` file needs to be compiled first, since
it provides the ``LIBLAMMPS`` module that is imported by the Fortran
code using the interface.  A working example code can be found together
with equivalent examples in C and C++ in the ``examples/COUPLE/simple``
folder of the LAMMPS distribution.

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

With the Fortran interface, the creation of a :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` instance is included in the constructor for
creating the :f:func:`lammps` derived type.  To import the definition of
that type and its type bound procedures, you need to add a ``USE
LIBLAMMPS`` statement.  Internally it will call either
:cpp:func:`lammps_open_fortran` or :cpp:func:`lammps_open_no_mpi` from
the C library API to create the class instance.  All arguments are
optional and :cpp:func:`lammps_mpi_init` will be called automatically,
if it is needed.  Similarly, a possible call to
:cpp:func:`lammps_mpi_finalize` is integrated into the :f:func:`close`
function and triggered with the optional logical argument set to
``.true.``. Here is a simple example:

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
thus make the resulting executable behave similarly to the standalone
executable (it will ignore the `-in/-i` flag, though).  This allows
using the command line to configure accelerator and suffix settings,
configure screen and logfile output, or to set index style variables
from the command line and more.  Here is a correspondingly adapted
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
*************************

Once a LAMMPS instance is created, it is possible to "drive" the LAMMPS
simulation by telling LAMMPS to read commands from a file or to pass
individual or multiple commands from strings or lists of strings.  This
is done similarly to how it is implemented in the :doc:`C-library
interface <Library_execute>`. Before handing off the calls to the
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

Accessing system properties
***************************

The C-library interface allows the :doc:`extraction of different kinds
of information <Library_properties>` about the active simulation
instance and also - in some cases - to apply modifications to it.  In
some cases, the C-library interface makes pointers to internal data
structures accessible, thus when accessing them from Fortran, special
care is needed to avoid data corruption and crashes.  Thus please see
the documentation of the individual type bound procedures for details.

Below is an example demonstrating some of the possible uses.

.. code-block:: fortran

  PROGRAM testprop
    USE LIBLAMMPS
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int64_t
    TYPE(lammps)     :: lmp
    INTEGER(kind=8)  :: natoms
    REAL(c_double), POINTER :: dt
    INTEGER(c_int64_t), POINTER :: ntimestep
    REAL(kind=8) :: pe, ke

    lmp = lammps()
    CALL lmp%file('in.sysinit')
    natoms = INT(lmp%get_natoms(),8)
    WRITE(6,'(A,I8,A)') 'Running a simulation with', natoms, ' atoms'
    WRITE(6,'(I8,A,I8,A,I3,A)') lmp%extract_setting('nlocal'), ' local and', &
        lmp%extract_setting('nghost'), ' ghost atom. ',                      &
        lmp%extract_setting('ntypes'), ' atom types'

    CALL lmp%command('run 2 post no')
    dt = lmp%extract_global('dt')
    ntimestep = lmp%extract_global('ntimestep')
    WRITE(6,'(A,I4,A,F4.1,A)') 'At step:', ntimestep, '  Changing timestep from', dt, ' to 0.5'
    dt = 0.5
    CALL lmp%command('run 2 post no')

    WRITE(6,'(A,I4)') 'At step:', ntimestep
    pe = lmp%get_thermo('pe')
    ke = lmp%get_thermo('ke')
    PRINT*, 'PE = ', pe
    PRINT*, 'KE = ', ke

    CALL lmp%close(.TRUE.)

  END PROGRAM testprop

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
   :f subroutine close: :f:func:`close`
   :f subroutine error: :f:func:`error`
   :f subroutine file: :f:func:`file`
   :f subroutine command: :f:func:`command`
   :f subroutine commands_list: :f:func:`commands_list`
   :f subroutine commands_string: :f:func:`commands_string`
   :f function get_natoms: :f:func:`get_natoms`
   :f function get_thermo: :f:func:`get_thermo`
   :f subroutine extract_box: :f:func:`extract_box`
   :f subroutine reset_box: :f:func:`reset_box`
   :f subroutine memory_usage: :f:func:`memory_usage`
   :f function extract_setting: :f:func:`extract_setting`
   :f function extract_global: :f:func:`extract_global`
   :f function version: :f:func:`version`
   :f function is_running: :f:func:`is_running`

--------

.. f:function:: lammps([args][,comm])

   This is the constructor for the Fortran class and will forward
   the arguments to a call to either :cpp:func:`lammps_open_fortran`
   or :cpp:func:`lammps_open_no_mpi`. If the LAMMPS library has been
   compiled with MPI support, it will also initialize MPI, if it has
   not already been initialized before.

   The *args* argument with the list of command line parameters is
   optional and so it the *comm* argument with the MPI communicator.
   If *comm* is not provided, ``MPI_COMM_WORLD`` is assumed. For
   more details please see the documentation of :cpp:func:`lammps_open`.

   :o character(len=\*) args(\*) [optional]: arguments as list of strings
   :o integer comm [optional]: MPI communicator
   :r lammps: an instance of the :f:type:`lammps` derived type

   .. note::

      The ``MPI_F08`` module, which defines Fortran 2008 bindings for MPI,
      is not directly supported by this interface due to the complexities of
      supporting both the ``MPI_F08`` and ``MPI`` modules at the same time.
      However, you should be able to use the ``MPI_VAL`` member of the
      ``MPI_comm`` derived type to access the integer value of the
      communicator, such as in

      .. code-block:: Fortran

         PROGRAM testmpi
           USE LIBLAMMPS
           USE MPI_F08
           TYPE(lammps) :: lmp
           lmp = lammps(MPI_COMM_SELF%MPI_VAL)
         END PROGRAM testmpi

Procedures Bound to the lammps Derived Type
===========================================

.. f:subroutine:: close([finalize])

   This method will close down the LAMMPS instance through calling
   :cpp:func:`lammps_close`.  If the *finalize* argument is present and
   has a value of ``.TRUE.``, then this subroutine also calls
   :cpp:func:`lammps_mpi_finalize`.

   :o logical finalize [optional]: shut down the MPI environment of the LAMMPS
    library if true.

--------

.. f:subroutine:: error(error_type, error_text)

   This method is a wrapper around the :cpp:func:`lammps_error` function and
   will dispatch an error through the LAMMPS Error class.

   .. versionadded:: TBD

   :p integer error_type: constant to select which Error class function to call
   :p character(len=\*) error_text: error message

--------

.. f:subroutine:: file(filename)

   This method will call :cpp:func:`lammps_file` to have LAMMPS read
   and process commands from a file.

   :p character(len=\*) filename: name of file with LAMMPS commands

--------

.. f:subroutine:: command(cmd)

   This method will call :cpp:func:`lammps_command` to have LAMMPS
   execute a single command.

   :p character(len=\*) cmd: single LAMMPS command

--------

.. f:subroutine:: commands_list(cmds)

   This method will call :cpp:func:`lammps_commands_list` to have LAMMPS
   execute a list of input lines.

   :p character(len=\*) cmd(:): list of LAMMPS input lines

--------

.. f:subroutine:: commands_string(str)

   This method will call :cpp:func:`lammps_commands_string` to have LAMMPS
   execute a block of commands from a string.

   :p character(len=\*) str: LAMMPS input in string

--------

.. f:function:: get_natoms()

   This function will call :cpp:func:`lammps_get_natoms` and return the number
   of atoms in the system.

   :r real(c_double): number of atoms

--------

.. f:function:: get_thermo(name)

   This function will call :cpp:func:`lammps_get_thermo` and return the value
   of the corresponding thermodynamic keyword.

   .. versionadded:: TBD

   :p character(len=\*) name: string with the name of the thermo keyword
   :r real(c_double): value of the requested thermo property or `0.0_c_double`

--------

.. f:subroutine:: extract_box([boxlo][, boxhi][, xy][, yz][, xz][, pflags][, boxflag])

   This subroutine will call :cpp:func:`lammps_extract_box`. All
   parameters are optional, though obviously at least one should be
   present. The parameters *pflags* and *boxflag* are stored in LAMMPS
   as integers, but should be declared as ``LOGICAL`` variables when
   calling from Fortran.

   .. versionadded:: TBD

   :o real(c_double) boxlo [dimension(3),optional]: vector in which to store
    lower-bounds of simulation box
   :o real(c_double) boxhi [dimension(3),optional]: vector in which to store
    upper-bounds of simulation box
   :o real(c_double) xy [optional]: variable in which to store *xy* tilt factor
   :o real(c_double) yz [optional]: variable in which to store *yz* tilt factor
   :o real(c_double) xz [optional]: variable in which to store *xz* tilt factor
   :o logical pflags [dimension(3),optional]: vector in which to store
    periodicity flags (``.TRUE.`` means periodic in that dimension)
   :o logical boxflag [optional]: variable in which to store boolean denoting
    whether the box will change during a simulation
    (``.TRUE.`` means box will change)

.. note::

   Note that a frequent use case of this function is to extract only one or
   more of the options rather than all seven. For example, assuming "lmp"
   represents a properly-initialized LAMMPS instance, the following code will
   extract the periodic box settings into the variable "periodic":

   .. code-block:: Fortran

      ! code to start up
      logical :: periodic(3)
      ! code to initialize LAMMPS / run things / etc.
      call lmp%extract_box(pflags = periodic)

--------

.. f:subroutine:: reset_box(boxlo, boxhi, xy, yz, xz)

   This subroutine will call :cpp:func:`lammps_reset_box`. All parameters
   are required.

   .. versionadded:: TBD

   :p real(c_double) boxlo [dimension(3)]: vector of three doubles containing
    the lower box boundary
   :p real(c_double) boxhi [dimension(3)]: vector of three doubles containing
    the upper box boundary
   :p real(c_double) xy: *x--y* tilt factor
   :p real(c_double) yz: *y--z* tilt factor
   :p real(c_double) xz: *x--z* tilt factor

--------

.. f:subroutine:: memory_usage(meminfo)

   This subroutine will call :cpp:func:`lammps_memory_usage` and store the
   result in the three-element array *meminfo*.

   .. versionadded:: TBD

   :p real(c_double) meminfo [dimension(3)]: vector of three doubles in which
    to store memory usage data

--------

.. f:function:: get_mpi_comm()

   This function returns a Fortran representation of the LAMMPS "world"
   communicator.

   .. versionadded:: TBD

   :r integer: Fortran integer equivalent to the MPI communicator LAMMPS is
    using

   .. note::

       The C library interface currently returns type ``int`` instead of
       type ``MPI_Fint``, which is the C type corresponding to Fortran
       ``INTEGER`` types of the default kind.  On most compilers, these
       are the same anyway, but this interface exchanges values this way
       to avoid warning messages.

   .. note::

      The ``MPI_F08`` module, which defines Fortran 2008 bindings for MPI,
      is not directly supported by this function.  However, you should be
      able to convert between the two using the `MPI_VAL` member of the
      communicator.  For example,

      .. code-block:: fortran

         USE MPI_F08
         USE LIBLAMMPS
         TYPE (lammps) :: lmp
         TYPE (MPI_Comm) :: comm
         ! ... [commands to set up LAMMPS/etc.]
         comm%MPI_VAL = lmp%get_mpi_comm()

      should assign an ``MPI_F08`` communicator properly.

--------

.. f:function:: extract_setting(keyword)

   Query LAMMPS about global settings. See the documentation for the
   :cpp:func:`lammps_extract_setting` function from the C library.

   .. versionadded:: TBD

   :p character(len=\*) keyword: string containing the name of the thermo keyword
   :r integer(c_int): value of the queried setting or :math:`-1` if unknown

--------

.. f:function:: extract_global(name)

   This function calls :cpp:func:`lammps_extract_global` and returns
   either a string or a pointer to internal global LAMMPS data,
   depending on the data requested through *name*.

   .. versionadded:: TBD

   Note that this function actually does not return a value, but rather
   associates the pointer on the left side of the assignment to point to
   internal LAMMPS data (with the exception of string data, which are
   copied and returned as ordinary Fortran strings). Pointers must be of
   the correct data type to point to said data (typically
   ``INTEGER(c_int)``, ``INTEGER(c_int64_t)``, or ``REAL(c_double)``)
   and have compatible kind and rank.  The pointer being associated with
   LAMMPS data is type-, kind-, and rank-checked at run-time via an
   overloaded assignment operator.  The pointers returned by this
   function are generally persistent; therefore it is not necessary to
   call the function again, unless a :doc:`clear` command has been
   issued, which wipes out and recreates the contents of the
   :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class.

   For example,

   .. code-block:: fortran

      PROGRAM demo
        USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t, c_int, c_double
        USE LIBLAMMPS
        TYPE(lammps) :: lmp
        INTEGER(c_int), POINTER :: nlocal => NULL()
        INTEGER(c_int64_t), POINTER :: ntimestep => NULL()
        REAL(c_double), POINTER :: dt => NULL()
        CHARACTER(LEN=10) :: units
        lmp = lammps()
        ! other commands
        nlocal = lmp%extract_global('nlocal')
        ntimestep = lmp%extract_global('ntimestep')
        dt = lmp%extract_global('dt')
        units = lmp%extract_global('units')
        ! more commands
        lmp.close(.TRUE.)
      END PROGRAM demo

   would extract the number of atoms on this processor, the current time step,
   the size of the current time step, and the units being used into the
   variables *nlocal*, *ntimestep*, *dt*, and *units*, respectively.

   .. note::

      If :f:func:`extract_global` returns a string, the string must have length
      greater than or equal to the length of the string (not including the
      terminal ``NULL`` character) that LAMMPS returns. If the variable's
      length is too short, the string will be truncated. As usual in Fortran,
      strings are padded with spaces at the end.

   :p character(len=\*) name: string with the name of the property to extract
   :r polymorphic: pointer to LAMMPS data. The left-hand side of the assignment
    should be either a string (if expecting string data) or a C-compatible
    pointer (e.g., ``INTEGER (c_int), POINTER :: nlocal``) to the extracted
    property. If expecting vector data, the pointer should have dimension ":".

   .. warning::

       Modifying the data in the location pointed to by the returned pointer
       may lead to inconsistent internal data and thus may cause failures,
       crashes, or bogus simulations.  In general, it is much better
       to use a LAMMPS input command that sets or changes these parameters.
       Using an input command will take care of all side effects and necessary
       updates of settings derived from such settings.

--------

.. f:function:: extract_atom(name)

   This function calls :c:func:`lammps_extract_atom` and returns a pointer to
   LAMMPS data tied to the :cpp:class:`Atom` class, depending on the data
   requested through *name*.

   Note that this function actually does not return a pointer, but rather
   associates the pointer on the left side of the assignment to point
   to internal LAMMPS data. Pointers must be of the correct type, kind, and
   rank (e.g., ``INTEGER(c_int), DIMENSION(:)`` for "type", "mask", or "tag";
   ``INTEGER(c_int64_t), DIMENSION(:)`` for "tag" if LAMMPS was compiled
   with the ``-DLAMMPS_BIGBIG`` flag; ``REAL(c_double), DIMENSION(:,:)`` for
   "x", "v", or "f"; and so forth). The pointer being associated with LAMMPS
   data is type-, kind-, and rank-checked at run-time. Pointers returned by
   this function are generally persistent; therefore, it is not necessary to
   call the function again unless the underlying LAMMPS data are destroyed,
   such as through the :doc:`clear` command.

   :p character(len=\*) name: string with the name of the property to extract
   :r polymorphic: pointer to LAMMPS data. The left-hand side of the assignment
    should be a C-interoperable pointer of appropriate kind and rank
    (e.g., ``INTEGER (c_int), POINTER :: mask(:)``) to the extracted
    property. If expecting vector data, the pointer should have dimension ":";
    if expecting matrix data, the pointer should have dimension ":,:".

    .. admonition:: Array index order

       Two-dimensional arrays returned from :f:func:`extract_atom` will be
       **transposed** from equivalent arrays in C, and they will be indexed
       from 1 instead of 0. For example, in C,
       
       .. code-block:: C

          void *lmp;
          double **x;
          /* more code to setup, etc. */
          x = lammps_extract_atom(lmp, "x");
          printf("%f\n", x[5][1]);

       will print the *y*-coordinate of the sixth atom on this processor.
       Conversely,

       .. code-block:: Fortran

          TYPE(lammps) :: lmp
          REAL(c_double), DIMENSION(:,:), POINTER :: x => NULL()
          ! more code to setup, etc.
          x = lmp%extract_atom("x")
          print '(f0.6)', x(2,6)

       will print the *y*-coordinate of the sixth atom on this processor
       (note the transposition of the two indices). This is not a choice, but
       rather a consequence of the different conventions adopted by the Fortran
       and C standards decades ago.

       If you would like the indices to start at 0 instead of 1 (which follows
       typical notation in C and C++, but not Fortran), you can create another
       pointer and associate it thus:
       
       .. code-block:: Fortran

          REAL(c_double), DIMENSION(:,:), POINTER :: x, x0
          x = lmp%extract_atom("x")
          x0(0:,0:) => x
  
       The above would cause the dimensions of *x* to be (1:3, 1:nmax)
       and those of *x0* to be (0:2, 0:nmax-1).

--------

.. f:function:: extract_compute(id, style, type)

   This function calls :c:func:`lammps_extract_compute` and returns a pointer
   to LAMMPS data tied to the :cpp:class:`Compute` class, specifically data
   provided by the compute identified by *id*. Computes may provide global,
   per-atom, or local data, and those data may be a scalar, a vector, or an
   array. Since computes may provide multiple kinds of data, the user is
   required to specify which set of data is to be returned through the
   *style* and *type* variables.

   Note that this function actually does not return a value, but rather
   associates the pointer on the left side of the assignment to point to
   internal LAMMPS data. Pointers must be of the correct data type to point to
   said data (i.e., ``REAL(c_double)``) and have compatible rank.  The pointer
   being associated with LAMMPS data is type-, kind-, and rank-checked at
   run-time via an overloaded assignment operator.

   For example,

   .. code-block:: Fortran

      TYPE(lammps) :: lmp
      REAL(c_double), DIMENSION(:), POINTER :: COM
      ! code to setup, create atoms, etc.
      CALL lmp%compute('compute COM all com')
      COM = lmp%extract_compute('COM', lmp%style%global, lmp%style%type)

   will bind the variable *COM* to the center of mass of the atoms created in
   your simulation. The vector in this case has length 3; the length (or, in
   the case of array data, the number of rows and columns) is determined for
   you based on data from the :cpp:class:`Compute` class.

   .. admonition:: Array index order

      Two-dimensional arrays returned from :f:func:`extract_compute` will be
      **transposed** from equivalent arrays in C, and they will be indexed
      from 1 instead of 0. See the similar note under
      :f:func:`extract_atom` for further details.

   The following combinations are possible (assuming ``lmp`` is the name of
   your LAMMPS instance):

   .. list-table::
      :header-rows: 1
      :widths: auto

      * - Style
        - Type
        - Pointer type to assign to
        - Returned data
      * - ``lmp%style%global``
        - ``lmp%type%scalar``
        - ``REAL(c_double), POINTER``
        - Global scalar
      * - ``lmp%style%global``
        - ``lmp%type%vector``
        - ``REAL(c_double), DIMENSION(:), POINTER``
        - Global vector
      * - ``lmp%style%global``
        - ``lmp%type%array``
        - ``REAL(c_double), DIMENSION(:,:), POINTER``
        - Global array
      * - ``lmp%style%atom``
        - ``lmp%type%vector``
        - ``REAL(c_double), DIMENSION(:), POINTER``
        - Per-atom vector
      * - ``lmp%style%atom``
        - ``lmp%type%array``
        - ``REAL(c_double), DIMENSION(:,:), POINTER``
        - Per-atom array
      * - ``lmp%style%local``
        - ``lmp%type%vector``
        - ``REAL(c_double), DIMENSION(:), POINTER``
        - Local vector
      * - ``lmp%style%local``
        - ``lmp%type%array``
        - ``REAL(c_double), DIMENSION(:,:), POINTER``
        - Local array

   :p character(len=\*) id: compute ID from which to extract data
   :p integer(c_int) style: value indicating the style of data to extract
    (global, per-atom, or local)
   :p integer(c_int) type: value indicating the type of data to extract
    (scalar, vector, or array)

   .. note::

      If the compute's data are not already computed for the current step, the
      compute will be invoked. LAMMPS cannot easily check at that time if it is
      valid to invoke a compute, so it may fail with an error. The caller has
      to check to avoid such an error.

   .. warning::

      The pointers returned by this function are generally not persistent,
      since the computed data may be re-distributed, re-allocated, and
      re-ordered at every invocation. It is advisable to re-invoke this
      function before the data are accessed or make a copy if the data are to
      be used after other LAMMPS commands have been issued. Do **not** modify
      the data returned by this function.

--------

.. f:function:: version()

   This method returns the numeric LAMMPS version like
   :cpp:func:`lammps_version` does.

   :r integer: LAMMPS version

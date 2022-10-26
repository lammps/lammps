The :f:mod:`LIBLAMMPS` Fortran Module
*************************************

The :f:mod:`LIBLAMMPS` module provides an interface to call LAMMPS from
Fortran.  It is based on the LAMMPS C library interface and requires a
fully Fortran 2003-compatible compiler to be compiled.  It is designed
to be self-contained and not require any support functions written in C,
C++, or Fortran other than those in the C library interface and the module
itself.

While C libraries have a defined binary interface (ABI) and can thus be
used from multiple compiler versions from different vendors as long
as they are compatible with the hosting operating system, the same is
not true for Fortran programs.  Thus, the LAMMPS Fortran module needs to be
compiled alongside the code using it from the source code in
``fortran/lammps.f90``.  When linking, you also need to
:doc:`link to the LAMMPS library <Build_link>`.  A typical command line
for a simple program using the Fortran interface would be:

.. code-block:: bash

   mpifort -o testlib.x lammps.f90 testlib.f90 -L. -llammps

Please note that the MPI compiler wrapper is only required when the
calling the library from an MPI-parallelized program.  Otherwise, using
the plain Fortran compiler (gfortran, ifort, flang, etc.) will suffice.
It may be necessary to link to additional libraries, depending on how
LAMMPS was configured and whether the LAMMPS library :doc:`was compiled
as a static or dynamic library <Build_link>`.

If the LAMMPS library itself has been compiled with MPI support, the
resulting executable will still be able to run LAMMPS in parallel with
``mpirun``, ``mpiexec`` or equivalent.  Please also note that the order
of the source files matters: the ``lammps.f90`` file needs to be
compiled first, since it provides the :f:mod:`LIBLAMMPS` module that is
imported by the Fortran code that uses the interface.  A working example
can be found together with equivalent examples in C and C++ in the
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
   closely resembles the C library interface is available in the
   ``examples/COUPLE/fortran2`` folder.  Please see the ``README`` file
   in that folder for more information about it and how to contact its
   author and maintainer.

----------

Creating or deleting a LAMMPS object
************************************

With the Fortran interface, the creation of a :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` instance is included in the constructor for
creating the :f:func:`lammps` derived type.  To import the definition of
that type and its type-bound procedures, you need to add a ``USE LIBLAMMPS``
statement.  Internally, it will call either
:cpp:func:`lammps_open_fortran` or :cpp:func:`lammps_open_no_mpi` from
the C library API to create the class instance.  All arguments are
optional and :cpp:func:`lammps_mpi_init` will be called automatically
if it is needed.  Similarly, a possible call to
:cpp:func:`lammps_mpi_finalize` is integrated into the :f:func:`close`
function and triggered with the optional logical argument set to
``.TRUE.``. Here is a simple example:

.. code-block:: fortran

   PROGRAM testlib
     USE LIBLAMMPS                 ! include the LAMMPS library interface
     IMPLICIT NONE
     TYPE(lammps) :: lmp           ! derived type to hold LAMMPS instance
     CHARACTER(LEN=*), PARAMETER :: args(3) = &
         [ CHARACTER(LEN=12) :: 'liblammps', '-log', 'none' ]

     ! create a LAMMPS instance (and initialize MPI)
     lmp = lammps(args)
     ! get and print numerical version code
     PRINT*, 'LAMMPS Version: ', lmp%version()
     ! delete LAMMPS instance (and shutdown MPI)
     CALL lmp%close(.TRUE.)
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
     TYPE(lammps) :: lmp           ! derived type to hold LAMMPS instance
     CHARACTER(LEN=128), ALLOCATABLE :: command_args(:)
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
is done similarly to how it is implemented in the :doc:`C library
interface <Library_execute>`. Before handing off the calls to the
C library interface, the corresponding Fortran versions of the calls
(:f:func:`file`, :f:func:`command`, :f:func:`commands_list`, and
:f:func:`commands_string`) have to make a copy of the strings passed as
arguments so that they can be modified to be compatible with the
requirements of strings in C without affecting the original strings.
Those copies are automatically deleted after the functions return.
Below is a small demonstration of the uses of the different functions:

.. code-block:: fortran

   PROGRAM testcmd
     USE LIBLAMMPS
     TYPE(lammps) :: lmp
     CHARACTER(LEN=512) :: cmds
     CHARACTER(LEN=40), ALLOCATABLE :: cmdlist(:)
     CHARACTER(LEN=10) :: trimmed
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

The C library interface allows the :doc:`extraction of different kinds
of information <Library_properties>` about the active simulation
instance and also---in some cases---to apply modifications to it, and the
Fortran interface provides access to the same data using Fortran-style,
C-interoperable data types.  In some cases, the Fortran library interface makes
pointers to internal LAMMPS data structures accessible; when accessing them
through the library interfaces, special care is needed to avoid data corruption
and crashes.  Please see the documentation of the individual type-bound
procedures for details.

Below is an example demonstrating some of the possible uses.

.. code-block:: fortran

  PROGRAM testprop
    USE LIBLAMMPS
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int64_t
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT
    TYPE(lammps) :: lmp
    INTEGER(KIND=c_int64_t), POINTER :: natoms
    REAL(KIND=c_double), POINTER :: dt
    INTEGER(KIND=c_int64_t), POINTER :: ntimestep
    REAL(KIND=c_double) :: pe, ke

    lmp = lammps()
    CALL lmp%file('in.sysinit')
    natoms = lmp%extract_global('natoms')
    WRITE(OUTPUT_UNIT,'(A,I0,A)') 'Running a simulation with ', natoms, ' atoms'
    WRITE(OUTPUT_UNIT,'(I0,A,I0,A,I0,A)') lmp%extract_setting('nlocal'), &
        ' local and ', lmp%extract_setting('nghost'), ' ghost atoms. ', &
        lmp%extract_setting('ntypes'), ' atom types'

    CALL lmp%command('run 2 post no')
    dt = lmp%extract_global('dt')
    ntimestep = lmp%extract_global('ntimestep')
    WRITE(OUTPUT_UNIT,'(A,I0,A,F4.1,A)') 'At step: ', ntimestep, &
        '  Changing timestep from', dt, ' to 0.5'
    dt = 0.5_c_double
    CALL lmp%command('run 2 post no')

    WRITE(OUTPUT_UNIT,'(A,I0)') 'At step: ', ntimestep
    pe = lmp%get_thermo('pe')
    ke = lmp%get_thermo('ke')
    PRINT*, 'PE = ', pe
    PRINT*, 'KE = ', ke

    CALL lmp%close(.TRUE.)
  END PROGRAM testprop

---------------

The :f:mod:`LIBLAMMPS` module API
*********************************

.. f:module:: LIBLAMMPS

Below are the detailed descriptions of definitions and interfaces
of the contents of the :f:mod:`LIBLAMMPS` Fortran interface to LAMMPS.

.. f:currentmodule:: None

.. f:type:: lammps

   Derived type that is the general class of the Fortran interface.
   It holds a reference to the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`
   class instance to which any of the included calls are forwarded.

   :f handle: reference to the LAMMPS class
   :ftype handle: c_ptr
   :f style: derived type to access lammps style constants
   :ftype style: type(lammps_style)
   :f type: derived type to access lammps type constants
   :ftype type: type(lammps_type) 
   :f close: :f:subr:`close`
   :ftype close: subroutine
   :f subroutine error: :f:subr:`error`
   :ftype error: subroutine
   :f file: :f:subr:`file`
   :ftype file: subroutine
   :f command: :f:subr:`command`
   :ftype command: subroutine
   :f commands_list: :f:subr:`commands_list`
   :ftype commands_list: subroutine
   :f commands_string: :f:subr:`commands_string`
   :ftype commands_string: subroutine
   :f get_natoms: :f:func:`get_natoms`
   :ftype get_natoms: function
   :f get_thermo: :f:func:`get_thermo`
   :ftype get_thermo: function
   :f extract_box: :f:subr:`extract_box`
   :ftype extract_box: subroutine
   :f reset_box: :f:subr:`reset_box`
   :ftype reset_box: subroutine
   :f memory_usage: :f:subr:`memory_usage`
   :ftype memory_usage: subroutine
   :f get_mpi_comm: :f:func:`get_mpi_comm`
   :ftype get_mpi_comm: function
   :f extract_setting: :f:func:`extract_setting`
   :ftype extract_setting: function
   :f extract_global: :f:func:`extract_global`
   :ftype extract_global: function
   :f extract_atom: :f:func:`extract_atom`
   :ftype extract_atom: function
   :f extract_compute: :f:func:`extract_compute`
   :ftype extract_compute: function
   :f extract_fix: :f:func:`extract_fix`
   :ftype extract_fix: function
   :f extract_variable: :f:func:`extract_variable`
   :ftype extract_variable: function
   :f gather_atoms: :f:subr:`gather_atoms`
   :ftype gather_atoms: subroutine
   :f gather_atoms_concat: :f:subr:`gather_atoms_concat`
   :ftype gather_atoms_concat: subroutine
   :f gather_atoms_subset: :f:subr:`gather_atoms_subset`
   :ftype gather_atoms_subset: subroutine
   :f scatter_atoms: :f:subr:`scatter_atoms`
   :ftype scatter_atoms: subroutine
   :f scatter_atoms_subset: :f:subr:`scatter_atoms_subset`
   :ftype scatter_atoms_subset: subroutine
   :f gather_bonds: :f:subr:`gather_bonds`
   :ftype gather_bonds: subroutine
   :f create_atoms: :f:subr:`create_atoms`
   :ftype create_atoms: subroutine
   :f find_pair_neighlist: :f:func:`find_pair_neighlist`
   :ftype find_pair_neighlist: function
   :f find_fix_neighlist: :f:func:`find_fix_neighlist`
   :ftype find_fix_neighlist: function
   :f find_compute_neighlist: :f:func:`find_compute_neighlist`
   :ftype find_compute_neighlist: function
   :f neighlist_num_elements: :f:func:`neighlist_num_elements`
   :ftype neighlist_num_elements: function
   :f version: :f:func:`version`
   :ftype version: function
   :f get_os_info: :f:subr:`get_os_info`
   :ftype get_os_info: subroutine
   :f config_has_mpi_support: :f:func:`config_has_mpi_support`
   :ftype config_has_mpi_support: function
   :f config_has_gzip_support: :f:func:`config_has_gzip_support`
   :ftype config_has_gzip_support: function
   :f config_has_png_support: :f:func:`config_has_png_support`
   :ftype config_has_png_support: function
   :f config_has_jpeg_support: :f:func:`config_has_jpeg_support`
   :ftype config_has_jpeg_support: function
   :f config_has_ffmpeg_support: :f:func:`config_has_ffmpeg_support`
   :ftype config_has_ffmpeg_support: function
   :f config_has_exceptions: :f:func:`config_has_exceptions`
   :ftype config_has_exceptions: function
   :f config_has_package: :f:func:`config_has_package`
   :ftype config_has_package: function
   :f config_package_count: :f:func:`config_package_count`
   :ftype config_package_count: function
   :f config_package_name: :f:func:`config_package_name`
   :ftype config_package_name: function
   :f installed_packages: :f:subr:`installed_packages`
   :ftype installed_packages: subroutine
   :f config_accelerator: :f:func:`config_accelerator`
   :ftype config_accelerator: function
   :f has_gpu_device: :f:func:`has_gpu_device`
   :ftype has_gpu_device: function
   :f get_gpu_device_info: :f:subr:`get_gpu_device_info`
   :ftype get_gpu_device_info: subroutine
   :f has_style: :f:func:`has_style`
   :ftype has_style: function
   :f style_count: :f:func:`style_count`
   :ftype style_count: function
   :f style_name: :f:func:`style_name`
   :ftype style_name: function
   :f has_id: :f:func:`has_id`
   :ftype has_id: function
   :f id_count: :f:func:`id_count`
   :ftype id_count: function
   :f id_name: :f:subr:`id_name`
   :ftype id_name: subroutine
   :f plugin_count: :f:func:`plugin_count`
   :ftype plugin_count: function
   :f plugin_name: :f:subr:`plugin_name`
   :ftype plugin_count: subroutine
   :f encode_image_flags: :f:func:`encode_image_flags`
   :ftype encode_image_flags: function
   :f decode_image_flags: :f:subr:`decode_image_flags`
   :ftype decode_image_flags: subroutine
   :f flush_buffers: :f:subr:`flush_buffers`
   :ftype flush_buffers: subroutine
   :f is_running: :f:func:`is_running`
   :ftype is_running: function
   :f force_timeout: :f:subr:`force_timeout`
   :ftype force_timeout: subroutine
   :f has_error: :f:func:`has_error`
   :ftype has_error: function
   :f get_last_error_message: :f:subr:`get_last_error_message`
   :ftype get_last_error_message: subroutine

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

   :o args: arguments as list of strings
   :otype args: character(len=\*),dimension(:),optional
   :o comm: MPI communicator
   :otype comm: integer,optional
   :to: :cpp:func:`lammps_open_fortran`
   :to: :cpp:func:`lammps_open_no_mpi`
   :r lammps: an instance of the :f:type:`lammps` derived type

   .. note::

      The :f:mod:`MPI_F08` module, which defines Fortran 2008 bindings for MPI,
      is not directly supported by this interface due to the complexities of
      supporting both the :f:mod:`MPI_F08` and :f:mod:`MPI` modules at the same
      time. However, you should be able to use the ``MPI_VAL`` member of the
      ``MPI_comm`` derived type to access the integer value of the
      communicator, such as in

      .. code-block:: fortran

         PROGRAM testmpi
           USE LIBLAMMPS
           USE MPI_F08
           TYPE(lammps) :: lmp
           lmp = lammps(MPI_COMM_SELF%MPI_VAL)
         END PROGRAM testmpi

.. f:type:: lammps_style

   This derived type is there to provide a convenient interface for the style
   constants used with :f:func:`extract_compute`, :f:func:`extract_fix`, and
   :f:func:`extract_variable`. Assuming your LAMMPS instance is called ``lmp``,
   these constants will be ``lmp%style%global``, ``lmp%style%atom``,
   and ``lmp%style%local``. These values are identical to the values described
   in :cpp:enum:`_LMP_STYLE_CONST` for the C library interface.

   :f integer(c_int) global: used to request global data
   :f integer(c_int) atom:   used to request per-atom data
   :f integer(c_int) local:  used to request local data

.. f:type:: lammps_type

   This derived type is there to provide a convenient interface for the type
   constants used with :f:func:`extract_compute`, :f:func:`extract_fix`, and
   :f:func:`extract_variable`. Assuming your LAMMPS instance is called ``lmp``,
   these constants will be ``lmp%type%scalar``, ``lmp%type%vector``, and
   ``lmp%type%array``. These values are identical to the values described
   in :cpp:enum:`_LMP_TYPE_CONST` for the C library interface.

   :f integer(c_int) scalar: used to request scalars
   :f integer(c_int) vector: used to request vectors
   :f integer(c_int) array:  used to request arrays (matrices)

Procedures Bound to the :f:type:`lammps` Derived Type
=====================================================

.. f:subroutine:: close([finalize])

   This method will close down the LAMMPS instance through calling
   :cpp:func:`lammps_close`.  If the *finalize* argument is present and
   has a value of ``.TRUE.``, then this subroutine also calls
   :cpp:func:`lammps_kokkos_finalize` and
   :cpp:func:`lammps_mpi_finalize`.

   :o finalize: shut down the MPI environment of the LAMMPS
    library if ``.TRUE.``.
   :otype finalize: logical,optional
   :to: :cpp:func:`lammps_close`
   :to: :cpp:func:`lammps_mpi_finalize`
   :to: :cpp:func:`lammps_kokkos_finalize`

--------

.. f:subroutine:: error(error_type, error_text)

   This method is a wrapper around the :cpp:func:`lammps_error` function and
   will dispatch an error through the LAMMPS Error class.

   .. versionadded:: TBD

   :p error_type: constant to select which Error class function to call
   :ptype error_type: integer(c_int)
   :p character(len=\*) error_text: error message
   :to: :cpp:func:`lammps_error`

--------

.. f:subroutine:: file(filename)

   This method will call :cpp:func:`lammps_file` to have LAMMPS read
   and process commands from a file.

   :p character(len=\*) filename: name of file with LAMMPS commands
   :to: :cpp:func:`lammps_file`

--------

.. f:subroutine:: command(cmd)

   This method will call :cpp:func:`lammps_command` to have LAMMPS
   execute a single command.

   :p character(len=\*) cmd: single LAMMPS command
   :to: :cpp:func:`lammps_command`

--------

.. f:subroutine:: commands_list(cmds)

   This method will call :cpp:func:`lammps_commands_list` to have LAMMPS
   execute a list of input lines.

   :p character(len=\*) cmd [dimension(:)]: list of LAMMPS input lines
   :to: :cpp:func:`lammps_commands_list`

--------

.. f:subroutine:: commands_string(str)

   This method will call :cpp:func:`lammps_commands_string` to have LAMMPS
   execute a block of commands from a string.

   :p character(len=\*) str: LAMMPS input in string
   :to: :cpp:func:`lammps_commands_string`

--------

.. f:function:: get_natoms()

   This function will call :cpp:func:`lammps_get_natoms` and return the number
   of atoms in the system.

   :to: :cpp:func:`lammps_get_natoms`
   :r real(c_double) natoms: number of atoms

   .. note::

      If you would prefer to get the number of atoms in its native format
      (i.e., as a 32- or 64-bit integer, depending on how LAMMPS was compiled),
      this can be extracted with :f:func:`extract_global`.

--------

.. f:function:: get_thermo(name)

   This function will call :cpp:func:`lammps_get_thermo` and return the value
   of the corresponding thermodynamic keyword.

   .. versionadded:: TBD

   :p character(len=\*) name: string with the name of the thermo keyword
   :to: :cpp:func:`lammps_get_thermo`
   :r value [real(c_double)]: value of the requested thermo property or `0.0_c_double`

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
   :o pflags: vector in which to store
    periodicity flags (``.TRUE.`` means periodic in that dimension)
   :otype pflags: logical,dimension(3),optional
   :o boxflag: variable in which to store boolean denoting
    whether the box will change during a simulation
    (``.TRUE.`` means box will change)
   :otype boxflag: logical,optional
   :to: :cpp:func:`lammps_extract_box`

.. note::

   Note that a frequent use case of this function is to extract only one or
   more of the options rather than all seven. For example, assuming "lmp"
   represents a properly-initialized LAMMPS instance, the following code will
   extract the periodic box settings into the variable "periodic":

   .. code-block:: fortran

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
   :to: :cpp:func:`lammps_reset_box`

--------

.. f:subroutine:: memory_usage(meminfo)

   This subroutine will call :cpp:func:`lammps_memory_usage` and store the
   result in the three-element array *meminfo*.

   .. versionadded:: TBD

   :p real(c_double) meminfo [dimension(3)]: vector of three doubles in which
    to store memory usage data
   :to: :cpp:func:`lammps_memory_usage`

--------

.. f:function:: get_mpi_comm()

   This function returns a Fortran representation of the LAMMPS "world"
   communicator.

   .. versionadded:: TBD

   :to: :cpp:func:`lammps_get_mpi_comm`
   :r comm: Fortran integer equivalent to the MPI communicator LAMMPS is
    using
   :rtype comm: integer

   .. note::

       The C library interface currently returns type ``int`` instead of
       type ``MPI_Fint``, which is the C type corresponding to Fortran
       ``INTEGER`` types of the default kind.  On most compilers, these
       are the same anyway, but this interface exchanges values this way
       to avoid warning messages.

   .. note::

      The :f:mod:`MPI_F08` module, which defines Fortran 2008 bindings for MPI,
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

      should assign an :f:mod:`MPI_F08` communicator properly.

--------

.. f:function:: extract_setting(keyword)

   Query LAMMPS about global settings. See the documentation for the
   :cpp:func:`lammps_extract_setting` function from the C library.

   .. versionadded:: TBD

   :p character(len=\*) keyword: string containing the name of the thermo keyword
   :to: :cpp:func:`lammps_extract_setting`
   :r integer(c_int) setting: value of the queried setting or :math:`-1` if
    unknown

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
      strings are padded with spaces at the end. If you use an allocatable
      string, the string **must be allocated** prior to calling this function,
      but you can automatically reallocate it to the correct length after the
      function returns, viz.,

      .. code-block :: fortran

         PROGRAM test
           USE LIBLAMMPS
           TYPE(lammps) :: lmp
           CHARACTER(LEN=:), ALLOCATABLE :: str
           lmp = lammps()
           CALL lmp%command('units metal')
           ALLOCATE ( CHARACTER(LEN=80) :: str )
           str = lmp%extract_global('units')
           str = TRIM(str) ! re-allocates to length len_trim(str) here
           PRINT*, LEN(str), LEN_TRIM(str)
         END PROGRAM test

      will print the number 5 (the length of the word "metal") twice.

   :p character(len=\*) name: string with the name of the property to extract
   :to: :cpp:func:`lammps_extract_global`
   :r pointer [polymorphic]: pointer to LAMMPS data. The left-hand side of the
    assignment should be either a string (if expecting string data) or a
    C-compatible pointer (e.g., ``INTEGER(c_int), POINTER :: nlocal``) to the
    extracted property. If expecting vector data, the pointer should have
    dimension ":".

   .. warning::

       Modifying the data in the location pointed to by the returned pointer
       may lead to inconsistent internal data and thus may cause failures,
       crashes, or bogus simulations.  In general, it is much better
       to use a LAMMPS input command that sets or changes these parameters.
       Using an input command will take care of all side effects and necessary
       updates of settings derived from such settings.

--------

.. f:function:: extract_atom(name)

   This function calls :cpp:func:`lammps_extract_atom` and returns a pointer to
   LAMMPS data tied to the :cpp:class:`Atom` class, depending on the data
   requested through *name*.

   .. versionadded:: TBD

   Note that this function actually does not return a pointer, but rather
   associates the pointer on the left side of the assignment to point
   to internal LAMMPS data. Pointers must be of the correct type, kind, and
   rank (e.g., ``INTEGER(c_int), DIMENSION(:)`` for "type", "mask", or "tag";
   ``INTEGER(c_int64_t), DIMENSION(:)`` for "tag" if LAMMPS was compiled
   with the ``-DLAMMPS_BIGBIG`` flag; ``REAL(c_double), DIMENSION(:,:)`` for
   "x", "v", or "f"; and so forth). The pointer being associated with LAMMPS
   data is type-, kind-, and rank-checked at run-time.

   :p character(len=\*) name: string with the name of the property to extract
   :to: :cpp:func:`lammps_extract_atom`
   :r pointer: pointer to LAMMPS data. The left-hand side of the
    assignment should be a C-interoperable pointer of appropriate kind and rank
    (e.g., ``INTEGER(c_int), POINTER :: mask(:)``) to the extracted
    property. If expecting vector data, the pointer should have dimension ":";
    if expecting matrix data, the pointer should have dimension ":,:".
   :rtype pointer: polymorphic

   .. warning::

      Pointers returned by this function are generally not persistent, as
      per-atom data may be redistributed, reallocated, and reordered at every
      re-neighboring operation. It is advisable to re-bind pointers using
      :f:func:`extract_atom` between runs.

   .. admonition:: Array index order

      Two-dimensional arrays returned from :f:func:`extract_atom` will be
      **transposed** from equivalent arrays in C, and they will be indexed
      from 1 instead of 0. For example, in C,

      .. code-block:: c

         void *lmp;
         double **x;
         /* more code to setup, etc. */
         x = lammps_extract_atom(lmp, "x");
         printf("%f\n", x[5][1]);

      will print the *y*-coordinate of the sixth atom on this processor.
      Conversely,

      .. code-block:: fortran

         TYPE(lammps) :: lmp
         REAL(c_double), DIMENSION(:,:), POINTER :: x => NULL()
         ! more code to setup, etc.
         x = lmp%extract_atom("x")
         PRINT '(f0.6)', x(2,6)

      will print the *y*-coordinate of the sixth atom on this processor
      (note the transposition of the two indices). This is not a choice, but
      rather a consequence of the different conventions adopted by the Fortran
      and C standards decades ago: in C, the block of data

      .. parsed-literal::

         1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

      interpreted as a :math:`4\times4` matrix would be

      .. math::

         \begin{bmatrix}
           1 & 2 & 3 & 4 \\
           5 & 6 & 7 & 8 \\
           9 & 10 & 11 & 12 \\
           13 & 14 & 15 & 16
         \end{bmatrix},

      that is, in row-major order. In Fortran, the same block of data is
      interpreted in column-major order, namely,

      .. math::

         \begin{bmatrix}
           1 & 5 & 9  & 13 \\
           2 & 6 & 10 & 14 \\
           3 & 7 & 11 & 15 \\
           4 & 8 & 12 & 16
         \end{bmatrix}.

      This difference in interpretation of the same block of data by the two
      languages means, in effect, that matrices from C or C++ will be
      transposed when interpreted in Fortran.

   .. note::

      If you would like the indices to start at 0 instead of 1 (which follows
      typical notation in C and C++, but not Fortran), you can create another
      pointer and associate it thus:

      .. code-block:: fortran

         REAL(c_double), DIMENSION(:,:), POINTER :: x, x0
         x = lmp%extract_atom("x")
         x0(0:,0:) => x

      The above would cause the dimensions of *x* to be (1:3, 1:nmax)
      and those of *x0* to be (0:2, 0:nmax\ :math:`-`\ 1).

--------

.. f:function:: extract_compute(id, style, type)

   This function calls :cpp:func:`lammps_extract_compute` and returns a pointer
   to LAMMPS data tied to the :cpp:class:`Compute` class, specifically data
   provided by the compute identified by *id*. Computes may provide global,
   per-atom, or local data, and those data may be a scalar, a vector, or an
   array. Since computes may provide multiple kinds of data, the user is
   required to specify which set of data is to be returned through the
   *style* and *type* variables.

   .. versionadded:: TBD

   Note that this function actually does not return a value, but rather
   associates the pointer on the left side of the assignment to point to
   internal LAMMPS data. Pointers must be of the correct data type to point to
   said data (i.e., ``REAL(c_double)``) and have compatible rank.  The pointer
   being associated with LAMMPS data is type-, kind-, and rank-checked at
   run-time via an overloaded assignment operator.

   For example,

   .. code-block:: fortran

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
      from 1 instead of 0. See the note at :f:func:`extract_atom` for
      further details.

   The following combinations are possible (assuming ``lmp`` is the name of
   your LAMMPS instance):

   .. list-table::
      :header-rows: 1
      :widths: auto

      * - Style
        - Type
        - Type to assign to
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
   :to: :cpp:func:`lammps_extract_compute`
   :r pointer: pointer to LAMMPS data. The left-hand side of the assignment
    should be a C-compatible pointer (e.g., ``REAL(c_double), POINTER :: x``)
    to the extracted property. If expecting vector data, the pointer should
    have dimension ":"; if expecting array (matrix) data, the pointer should
    have dimension ":,:".
   :rtype pointer: polymorphic

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

.. f:function:: extract_fix(id, style, type[, nrow][, ncol])

   This function calls :cpp:func:`lammps_extract_fix` and returns a pointer to
   LAMMPS data tied to the :cpp:class:`Fix` class, specifically data provided
   by the fix identified by *id*. Fixes may provide global, per-atom, or
   local data, and those data may be a scalar, a vector, or an array. Since
   many fixes provide multiple kinds of data, the user is required to specify
   which set of data is to be returned through the *style* and *type*
   variables.

   .. versionadded:: TBD

   Global data are calculated at the time they are requested and are only
   available element-by-element. As such, the user is expected to provide
   the *nrow* variable to specify which element of a global vector or the
   *nrow* and *ncol* variables to specify which element of a global array the
   user wishes LAMMPS to return. The *ncol* variable is optional for global
   scalar or vector data, and both *nrow* and *ncol* are optional when a
   global scalar is requested, as well as when per-atom or local data are
   requested. The following combinations are possible (assuming ``lmp`` is the
   name of your LAMMPS instance):

   .. list-table::
      :header-rows: 1
      :widths: auto

      * - Style
        - Type
        - nrow
        - ncol
        - Type to assign to
        - Returned data
      * - ``lmp%style%global``
        - ``lmp%type%scalar``
        - Ignored
        - Ignored
        - ``REAL(c_double)``
        - Global scalar
      * - ``lmp%style%global``
        - ``lmp%type%vector``
        - Required
        - Ignored
        - ``REAL(c_double)``
        - Element of global vector
      * - ``lmp%style%global``
        - ``lmp%type%array``
        - Required
        - Required
        - ``REAL(c_double)``
        - Element of global array
      * - ``lmp%style%atom``
        - ``lmp%type%scalar``
        -
        -
        -
        - (not allowed)
      * - ``lmp%style%atom``
        - ``lmp%type%vector``
        - Ignored
        - Ignored
        - ``REAL(c_double), DIMENSION(:), POINTER``
        - Per-atom vector
      * - ``lmp%style%atom``
        - ``lmp%type%array``
        - Ignored
        - Ignored
        - ``REAL(c_double), DIMENSION(:,:), POINTER``
        - Per-atom array
      * - ``lmp%style%local``
        - ``lmp%type%scalar``
        -
        -
        -
        - (not allowed)
      * - ``lmp%style%local``
        - ``lmp%type%vector``
        - Ignored
        - Ignored
        - ``REAL(c_double), DIMENSION(:), POINTER``
        - Per-atom vector
      * - ``lmp%style%local``
        - ``lmp%type%array``
        - Ignored
        - Ignored
        - ``REAL(c_double), DIMENSION(:,:), POINTER``
        - Per-atom array

   In the case of global data, this function returns a value of type
   ``real(c_double)``. For per-atom or local data, this function does not
   return a value but instead associates the pointer on the left side of the
   assignment to point to internal LAMMPS data. Pointers must be of the correct
   data type to point to said data (i.e., ``REAL(c_double)``) and have
   compatible rank.  The pointer being associated with LAMMPS data is type-,
   kind-, and rank-checked at run-time via an overloaded assignment operator.

   For example,

   .. code-block:: fortran

      TYPE(lammps) :: lmp
      REAL(c_double) :: dr, dx, dy, dz
      ! more code to set up, etc.
      lmp%command('fix george all recenter 2 2 2')
      ! more code
      dr = lmp%extract_fix("george", lmp%style%global, lmp%style%scalar)
      dx = lmp%extract_fix("george", lmp%style%global, lmp%style%vector, 1)
      dy = lmp%extract_fix("george", lmp%style%global, lmp%style%vector, 2)
      dz = lmp%extract_fix("george", lmp%style%global, lmp%style%vector, 3)

   will extract the global scalar calculated by
   :doc:`fix recenter <fix_recenter>` into the variable *dr* and the
   three elements of the global vector calculated by fix recenter into the
   variables *dx*, *dy*, and *dz*, respectively.

   If asked for per-atom or local data, :f:func:`extract_compute` returns a
   pointer to actual LAMMPS data. The pointer so returned will have the
   appropriate size to match the internal data, and will be
   type/kind/rank-checked at the time of the assignment. For example,

   .. code-block:: fortran

      TYPE(lammps) :: lmp
      REAL(c_double), DIMENSION(:), POINTER :: r
      ! more code to set up, etc.
      lmp%command('fix state all store/state 0 x y z')
      ! more code
      r = lmp%extract_fix('state', lmp%style%atom, lmp%type%array)

   will bind the pointer *r* to internal LAMMPS data representing the per-atom
   array computed by :doc:`fix store/state <fix_store_state>` when three
   inputs are specified. Similarly,

   .. code-block:: fortran

      TYPE(lammps) :: lmp
      REAL(c_double), DIMENSION(:), POINTER :: x
      ! more code to set up, etc.
      lmp%command('fix state all store/state 0 x')
      ! more code
      x = lmp%extract_fix('state', lmp%style%atom, lmp%type%vector)

   will associate the pointer *x* with internal LAMMPS data corresponding to
   the per-atom vector computed by :doc:`fix store/state <fix_store_state>`
   when only one input is specified. Similar examples with ``lmp%style%atom``
   replaced by ``lmp%style%local`` will extract local data from fixes that
   define local vectors and/or arrays.

   .. warning::

      The pointers returned by this function for per-atom or local data are
      generally not persistent, since the computed data may be redistributed,
      reallocated, and reordered at every invocation of the fix.  It is thus
      advisable to re-invoke this function before the data are accessed or to
      make a copy if the data are to be used after other LAMMPS commands have
      been issued.

   .. note::

      LAMMPS cannot easily check if it is valid to access the data, so it
      may fail with an error.  The caller has to avoid such an error.

   :p character(len=\*) id: string with the name of the fix from which
    to extract data
   :p integer(c_int) style: value indicating the style of data to extract
    (global, per-atom, or local)
   :p integer(c_int) type: value indicating the type of data to extract
    (scalar, vector, or array)
   :p integer(c_int) nrow: row index (used only for global vectors and arrays)
   :p integer(c_int) ncol: column index (only used for global arrays)
   :to: :cpp:func:`lammps_extract_fix`
   :r data: LAMMPS data (for global data) or a pointer to LAMMPS data
    (for per-atom or local data). The left-hand side of the assignment should
    be of type ``REAL(c_double)`` and have appropriate rank (i.e.,
    ``DIMENSION(:)`` if expecting per-atom or local vector data and
    ``DIMENSION(:,:)`` if expecting per-atom or local array data). If expecting
    local or per-atom data, it should have the ``POINTER`` attribute, but
    if expecting global data, it should be an ordinary (non-``POINTER``)
    variable.
   :rtype data: polymorphic

   .. admonition:: Array index order

      Two-dimensional global, per-atom, or local array data from
      :f:func:`extract_fix` will be **transposed** from equivalent arrays in
      C (or in the ordinary LAMMPS interface accessed through thermodynamic
      output), and they will be indexed from 1, not 0. This is true even for
      global data, which are returned as scalars---this is done primarily so
      the interface is consistent, as there is no choice but to transpose the
      indices for per-atom or local array data. See the similar note under
      :f:func:`extract_atom` for further details.

--------

.. f:function:: extract_variable(name[,group])

   This function calls :cpp:func:`lammps_extract_variable` and returns a scalar,
   vector, or string containing the value of the variable identified by
   *name*. When the variable is an *equal*-style variable (or one compatible
   with that style such as *internal*), the variable is evaluated and the
   corresponding value returned. When the variable is an *atom*-style variable,
   the variable is evaluated and a vector of values is returned. With all
   other variables, a string is returned. The *group* argument is only used
   for *atom* style variables and is ignored otherwise. If *group* is absent
   for *atom*-style variables, the group is assumed to be "all".

   .. versionadded:: TBD

   This function returns the values of the variables, not pointers to them.
   Vectors pointing to *atom*-style variables should be of type
   ``REAL(c_double)``, be of rank 1 (i.e., ``DIMENSION(:)``), and have the
   ``ALLOCATABLE`` attribute.

   .. note::

      Unlike the C library interface, the Fortran interface does not require
      you to deallocate memory when you are through; this is done for you,
      behind the scenes.

   For example,

   .. code-block:: fortran

      TYPE(lammps) :: lmp
      REAL(c_double) :: area
      ! more code to set up, etc.
      lmp%command('variable A equal lx*ly')
      ! more code
      area = lmp%extract_variable("A")

   will extract the *x*\ --*y* cross-sectional area of the simulation into the
   variable *area*.

   :p character(len=\*) name: variable name to evaluate
   :o character(len=\*) group [optional]: group for which to extract per-atom
    data (if absent, use "all")
   :to: :cpp:func:`lammps_extract_variable`
   :r data: scalar of type ``REAL(c_double)`` (for *equal*-style
    variables and others that are *equal*-compatible), vector of type
    ``REAL(c_double), DIMENSION(:), ALLOCATABLE`` for *atom*- or *vector*-style
    variables, or ``CHARACTER(LEN=*)`` for *string*-style and compatible
    variables. Strings whose length is too short to hold the result will be
    truncated. Allocatable strings must be allocated before this function is
    called; see note at :f:func:`extract_global` regarding allocatable strings.
    Allocatable arrays (for *atom*- and *vector*-style data) will be
    reallocated on assignment.
   :rtype data: polymorphic

.. note::

   LAMMPS cannot easily check if it is valid to access the data
   referenced by the variables (e.g., computes, fixes, or thermodynamic
   info), so it may fail with an error.  The caller has to make certain
   that the data are extracted only when it is safe to evaluate the variable
   and thus an error and crash are avoided.

--------

.. f:subroutine:: gather_atoms(name, count, data)

   This function calls :cpp:func:`lammps_gather_atoms` to gather the named
   atom-based entity for all atoms on all processors and return it in the
   vector *data*. The vector *data* will be ordered by atom
   ID, which requires consecutive atom IDs (1 to *natoms*).

   .. versionadded:: TBD

   If you need a similar array but have non-consecutive atom IDs, see
   :f:func:`gather_atoms_concat`; for a similar array but for a subset
   of atoms, see :f:func:`gather_atoms_subset`.

   The *data* array will be ordered in groups of *count* values, sorted by atom
   ID (e.g., if *name* is *x* and *count* = 3, then *data* = x[1][1], x[2][1],
   x[3][1], x[1][2], x[2][2], x[3][2], x[1][3], :math:`\dots`);
   *data* must be ``ALLOCATABLE`` and will be allocated to length
   (*count* :math:`\times` *natoms*), as queried by
   :f:func:`extract_setting`.

   :p character(len=\*) name: desired quantity (e.g., *x* or *mask*)
   :p integer(c_int) count: number of per-atom values you expect per atom
    (e.g., 1 for *type*, *mask*, or *charge*; 3 for *x*, *v*, or *f*). Use
    *count* = 3 with *image* if you want a single image flag unpacked into
    *x*/*y*/*z* components.
   :p real(c_double) data [dimension(:),allocatable]: array into which to store
    the data. Array *must* have the ``ALLOCATABLE`` attribute and be of rank 1
    (i.e., ``DIMENSION(:)``). If this array is already allocated, it will be
    reallocated to fit the length of the incoming data.
   :to: :cpp:func:`lammps_gather_atoms`

   .. note::

      If you want data from this function to be accessible as a two-dimensional
      array, you can declare a rank-2 pointer and reassign it, like so:

      .. code-block:: fortran

         USE, INTRINSIC :: ISO_C_BINDING
         USE LIBLAMMPS
         TYPE(lammps) :: lmp
         REAL(c_double), DIMENSION(:), ALLOCATABLE, TARGET :: xdata
         REAL(c_double), DIMENSION(:,:), POINTER :: x
         ! other code to set up, etc.
         CALL lmp%gather_atoms('x',3,xdata)
         x(1:3,1:size(xdata)/3) => xdata

      You can then access the *y*\ -component of atom 3 with ``x(2,3)``.
      See the note about array index order at :f:func:`extract_atom`.

--------

.. f:subroutine:: gather_atoms_concat(name, count, data)

   This function calls :cpp:func:`lammps_gather_atoms_concat` to gather the
   named atom-based entity for all atoms on all processors and return it in the
   vector *data*.

   .. versionadded:: TBD

   The vector *data* will not be ordered by atom ID, and there is no
   restriction on the IDs being consecutive. If you need the IDs, you can do
   another :f:func:`gather_atoms_concat` with *name* set to ``id``.

   If you need a similar array but have consecutive atom IDs, see
   :f:func:`gather_atoms`; for a similar array but for a subset of atoms, see
   :f:func:`gather_atoms_subset`.

   :p character(len=\*) name: desired quantity (e.g., *x* or *mask*)
   :p integer(c_int) count: number of per-atom values you expect per atom
    (e.g., 1 for *type*, *mask*, or *charge*; 3 for *x*, *v*, or *f*). Use
    *count* = 3 with *image* if you want a single image flag unpacked into
    *x*/*y*/*z* components.
   :p real(c_double) data [dimension(:),allocatable]: array into which to store
    the data. Array *must* have the ``ALLOCATABLE`` attribute and be of rank 1
    (i.e., ``DIMENSION(:)``). If this array is already allocated, it will be
    reallocated to fit the length of the incoming data.
   :to: :cpp:func:`lammps_gather_atoms_concat`

--------

.. f:subroutine:: gather_atoms_subset(name, count, ids, data)

   This function calls :cpp:func:`lammps_gather_atoms_subset` to gather the
   named atom-based entity for the atoms in the array *ids* from all processors
   and return it in the vector *data*.

   .. versionadded: TBD

   This subroutine gathers data for the requested atom IDs and stores them in a
   one-dimensional array allocated by the user. The data will be ordered by
   atom ID, but there is no requirement that the IDs be consecutive. If you
   wish to return a similar array for *all* the atoms, use
   :f:func:`gather_atoms` or :f:func:`gather_atoms_concat`.

   The *data* array will be in groups of *count* values, sorted by atom ID
   in the same order as the array *ids* (e.g., if *name* is *x*, *count* = 3,
   and *ids* is [100, 57, 210], then *data* might look like
   [x(1,100), x(2,100), x(3,100), x(1,57), x(2,57), x(3,57), x(1,210),
   :math:`\dots`]; *ids* must be provided by the user, and *data* must be
   of rank 1 (i.e., ``DIMENSION(:)``) and have the ``ALLOCATABLE`` attribute.

   :p character(len=\*) name: desired quantity (e.g., *x* or *mask*)
   :p integer(c_int) count: number of per-atom values you expect per atom
    (e.g., 1 for *type*, *mask*, or *charge*; 3 for *x*, *v*, or *f*). Use
    *count* = 3 with *image* if you want a single image flag unpacked into
    *x*/*y*/*z* components.
   :p integer(c_int) ids [dimension(:)]: atom IDs corresponding to the atoms
    to be gathered
   :p real(c_double) data [dimension(:),allocatable]: array into which to store
    the data. Array *must* have the ``ALLOCATABLE`` attribute and be of rank 1
    (i.e., ``DIMENSION(:)``). If this array is already allocated, it will be
    reallocated to fit the length of the incoming data.
   :to: :cpp:func:`lammps_gather_atoms_subset`

--------

.. f:subroutine:: scatter_atoms(name, data)

   This function calls :cpp:func:`lammps_scatter_atoms` to scatter the named
   atom-based entities in *data* to all processors.

   .. versionadded:: TBD

   This subroutine takes data stored in a one-dimensional array supplied by the
   user and scatters them to all atoms on all processors. The data must be
   ordered by atom ID, with the requirement that the IDs be consecutive.
   Use :f:func:`scatter_atoms_subset` to scatter data for some (or all)
   atoms, in any order.

   The *data* array needs to be ordered in groups of *count* values, sorted by
   atom ID (e.g., if *name* is *x* and *count* = 3, then
   *data* = [x(1,1) x(2,1) x(3,1) x(1,2) x(2,2) x(3,2) x(1,3) :math:`\dots`];
   *data* must be of length (*count* :math:`\times` *natoms*).

   :p character(len=\*) name: quantity to be scattered (e.g., *x* or *charge*)
   :p data: per-atom values packed in a one-dimensional array
    containing the data to be scattered. This array must have length *natoms*
    (e.g., for *type* or *charge*) or length *natoms*\ :math:`\times 3`
    (e.g., for *x* or *f*). The array *data* must be rank 1 (i.e.,
    ``DIMENSION(:)``) and be of type ``INTEGER(c_int)`` (e.g., for *mask* or
    *type*) or of type ``REAL(c_double)`` (e.g., for *x* or *charge* or *f*).
   :ptype data: polymorphic,dimension(:)
   :to: :cpp:func:`lammps_scatter_atoms`

--------

.. f:subroutine:: scatter_atoms_subset(name, ids, data)

   This function calls :cpp:func:`lammps_scatter_atoms_subset` to scatter the
   named atom-based entities in *data* to all processors.

   .. versionadded:: TBD

   This subroutine takes data stored in a one-dimensional array supplied by the
   user and scatters them to a subset of atoms on all processors. The array
   *data* contains data associated with atom IDs, but there is no requirement
   that the IDs be consecutive, as they are provided in a separate array,
   *ids*. Use :f:func:`scatter_atoms` to scatter data for all atoms, in order.

   The *data* array needs to be organized in groups of 1 or 3 values,
   depending on which quantity is being scattered, with the groups in the same
   order as the array *ids*. For example, if you want *data* to be the array
   [x(1,1) x(2,1) x(3,1) x(1,100) x(2,100) x(3,100) x(1,57) x(2,57) x(3,57)],
   then *ids* would be [1 100 57] and *name* would be *x*.

   :p character(len=\*) name: quantity to be scattered (e.g., *x* or *charge*)
   :p integer(c_int) ids [dimension(:)]: atom IDs corresponding to the atoms
    being scattered
   :p data: per-atom values packed into a
    one-dimensional array containing the data to be scattered. This array must
    have either the same length as *ids* (for *mask*, *type*, etc.) or three
    times its length (for *x*, *f*, etc.); the array must be rank 1
    and be of type ``INTEGER(c_int)`` (e.g., for *mask* or *type*) or of type
    ``REAL(c_double)`` (e.g., for *charge*, *x*, or *f*).
   :ptype data: polymorphic,dimension(:)
   :to: :cpp:func:`lammps_scatter_atoms_subset`

--------

.. f:subroutine:: gather_bonds(data)

   Gather type and constituent atom information for all bonds.

   .. versionadded:: TBD

   This function copies the list of all bonds into an allocated array.
   The array will be filled with (bond type, bond atom 1, bond atom 2) for each
   bond. The array is allocated to the right length (i.e., three times the
   number of bonds). The array *data* must be of the same type as the LAMMPS
   ``tagint`` type, which is equivalent to either ``INTEGER(c_int)`` or
   ``INTEGER(c_int64_t)``, depending on whether ``-DLAMMPS_BIGBIG`` was used
   when LAMMPS was built. If the supplied array does not match, an error will
   result at run-time.

   An example of how to use this routine is below:

   .. code-block:: fortran

      PROGRAM bonds
        USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT
        USE LIBLAMMPS
        IMPLICIT NONE
        INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: bonds_array
        INTEGER(c_int), DIMENSION(:,:), POINTER :: bonds
        TYPE(lammps) :: lmp
        INTEGER :: i
        ! other commands to initialize LAMMPS, create bonds, etc.
        CALL lmp%gather_bonds(bonds)
        bonds(1:3,1:size(bonds)/3) => bonds_array
        DO i = 1, size(bonds)/3
          WRITE(OUTPUT_UNIT,'(A,1X,I4,A,I4,1X,I4)') 'bond', bonds(1,i), &
            '; type = ', bonds(2,i), bonds(3,i)
        END DO
      END PROGRAM bonds

   :p data: array into which to copy the result. \*The ``KIND`` parameter is
    either ``c_int`` or, if LAMMPS was compiled with ``-DLAMMPS_BIGBIG``,
    kind ``c_int64_t``.
   :ptype data: integer(kind=\*),allocatable
   :to: :cpp:func:`lammps_gather_bonds`

--------

.. f:subroutine:: create_atoms([id,] type, x, [v,] [image,] [bexpand])

   This method calls :cpp:func:`lammps_create_atoms` to create additional atoms
   from a given list of coordinates and a list of atom types. Additionally,
   the atom IDs, velocities, and image flags may be provided.

   .. versionadded:: TBD

   :p integer(c_int) type [dimension(N)]: vector of :math:`N` atom types
    (required/see note below)
   :p real(c_double) x [dimension(3N)]: vector of :math:`3N` x/y/z positions
    of the new atoms, arranged as :math:`[x_1,y_1,z_1,x_2,y_2,\dotsc]`
    (required/see note below)
   :o integer(kind=\*) id [dimension(N)]: vector of :math:`N` atom IDs; if
    absent, LAMMPS will generate them for you. \*The ``KIND`` parameter should
    be ``c_int`` unless LAMMPS was compiled with ``-DLAMMPS_BIGBIG``, in which
    case it should be ``c_int64_t``.
   :o real(c_double) v [dimension(3N)]: vector of :math:`3N` x/y/z velocities
    of the new atoms, arranged as :math:`[v_{1,x},v_{1,y},v_{1,z},v_{2,x},
    \dotsc]`; if absent, they will be set to zero
   :o integer(kind=\*) image [dimension(N)]: vector of :math:`N` image flags;
    if absent, they are set to zero. \*The ``KIND`` parameter should be
    ``c_int`` unless LAMMPS was compiled with ``-DLAMMPS_BIGBIG``, in which
    case it should be ``c_int64_t``. See note below.
   :o bexpand: if ``.TRUE.``, atoms outside of shrink-wrap boundaries
    will be created, not dropped, and the box dimensions will be extended.
    Default is ``.FALSE.``
   :otype bexpand: logical
   :to: :cpp:func:`lammps_create_atoms`

   .. note::

      The *type* and *x* arguments are required, but they are declared
      ``OPTIONAL`` in the module because making them mandatory would require
      *id* to be present as well. To have LAMMPS generate the ids for you,
      use a call something like

      .. code-block:: fortran

         lmp%create_atoms(type=new_types, x=new_xs)

   .. note::

      When LAMMPS has been compiled with ``-DLAMMPS_BIGBIG``, it is not
      possible to include the *image* parameter but omit the *id* parameter.
      Either *id* must be present, or both *id* and *image* must be absent.
      This is required because having all arguments be optional in both
      generic functions creates an ambiguous interface. This limitation does
      not exist if LAMMPS was not compiled with ``-DLAMMPS_BIGBIG``.

--------

.. f:function:: find_pair_neighlist(style, exact, nsub, reqid)

   Find index of a neighbor list requested by a pair style.

   .. versionadded:: TBD

   This function determines which of the available neighbor lists for pair
   styles matches the given conditions.  It first matches the style name.
   If *exact* is ``.TRUE.``, the name must match exactly; if ``.FALSE.``, a
   regular expression or sub-string match is done.  If the pair style is
   *hybrid* or *hybrid/overlay*, the style is matched against the sub-styles
   instead. If the same pair style is used multiple times as a sub-style, the
   *nsub* argument must be :math:`> 0`; this argument represents the *n*\ th
   instance of the sub-style (same as for the :doc:`pair_coeff <pair_coeff>`
   command, for example). In that case, *nsub*\ :math:`{} = 0` will not
   produce a match, and the function will return :math:`-1`.

   The final condition to be checked is the request ID (\ *reqid*\ ). This
   will usually be zero, but some pair styles request multiple neighbor
   lists and set the request ID to a value greater than zero.

   :p character(len=\*) style: String used to search for pair style instance.
   :p exact: Flag to control whther style should match exactly or only a
    regular expression/sub-string match is applied.
   :ptype exact: logical
   :p integer(c_int) nsub:     Match *nsub*\ th hybrid sub-style instance of
    the same style
   :p integer(c_int) reqid:    Request ID to identify the neighbor list in
    case there are multiple requests from the same pair style instance.
   :to: :cpp:func:`lammps_find_pair_neighlist`
   :r integer(c_int) index:    Neighbor list index if found, otherwise
    :math:`-1`.

--------

.. f:function:: find_fix_neighlist()

   Find index of a neighbor list requested by a fix.

   .. versionadded:: TBD

   The neighbor list request from a fix is identified by the fix ID and the
   request ID. The request ID is typically zero, but will be :math:`>0` for
   fixes with multiple neighbor list requests.

   :p character(len=\*) id: Identifier of fix instance
   :p integer(c_int) reqid: request ID to identify the neighbor list in cases
    in which there are multiple requests from the same fix.
   :to: :cpp:func:`lammps_find_fix_neighlist`
   :r index: neighbor list index if found, otherwise :math:`-1`
   :rtype index: integer(c_int)

--------

.. f:function:: find_compute_neighlist()

   Find index of a neighbor list requested by a compute.

   .. versionadded:: TBD

   The neighbor list request from a compute is identified by the compute ID and
   the request ID.  The request ID is typically zero, but will be :math:`> 0`
   in case a compute has multiple neighbor list requests.

   :p character(len=\*) id: Identifier of compute instance
   :p integer(c_int) reqid: request ID to identify the neighbor list in cases
    in which there are multiple requests from the same compute
   :to: :cpp:func:`lammps_find_compute_neighlist`
   :r index: neighbor list index if found, otherwise :math:`-1`
   :rtype index: integer(c_int)

--------

.. f:function:: neighlist_num_elements(idx)

   Return the number of entries in the neighbor list with the given index.

   .. versionadded:: TBD

   :p integer(c_int) idx: neighbor list index
   :to: :cpp:func:`lammps_neighlist_num_elements`
   :r inum: number of entries in neighbor list, or :math:`-1` if *idx* is not
    a valid index.
   :rtype inum: integer(c_int)

--------

.. f:function:: version()

   This method returns the numeric LAMMPS version like
   :cpp:func:`lammps_version` does.

   :to: :cpp:func:`lammps_version`
   :r version: LAMMPS version
   :rtype version: integer

--------

.. f:subroutine:: get_os_info(buffer)

   This function can be used to retrieve detailed information about the hosting
   operating system and compiler/runtime environment.

   .. versionadded:: TBD

   A suitable buffer has to be provided. The assembled text will be truncated
   so as not to overflow this buffer. The string is typically a few hundred
   bytes long.

   :p character(len=\*) buffer: string that will house the information.
   :to: :cpp:func:`lammps_get_os_info`


--------

.. f:function:: config_has_mpi_support()

   This function is used to query whether LAMMPS was compiled with a real MPI
   library or in serial.

   .. versionadded:: TBD

   :to: :cpp:func:`lammps_config_has_mpi_support`
   :r has_mpi: ``.FALSE.`` when compiled with STUBS, ``.TRUE.`` if
    complied with MPI.
   :rtype has_mpi: logical

--------

.. f:function:: config_has_gzip_support()

   Check if the LAMMPS library supports reading or writing compressed
   files via a pipe to gzip or similar compression programs.

   .. versionadded:: TBD

   Several LAMMPS commands (e.g., :doc:`read_data`, :doc:`write_data`,
   :doc:`dump styles atom, custom, and xyz <dump>`) support reading and writing
   compressed files via creating a pipe to the ``gzip`` program.  This function
   checks whether this feature was :ref:`enabled at compile time <gzip>`.
   It does **not** check whether ``gzip`` or any other supported compression
   programs themselves are installed and usable.

   :to: :cpp:func:`lammps_config_has_gzip_support`
   :r has_gzip:
   :rtype has_gzip: logical

--------

.. f:function:: config_has_png_support()

   Check if the LAMMPS library supports writing PNG format images.

   .. versionadded:: TBD

   The LAMMPS :doc:`dump style image <dump_image>` supports writing multiple
   image file formats.  Most of them, however, need support from an external
   library, and using that has to be :ref:`enabled at compile time <graphics>`.
   This function checks whether support for the `PNG image file format
   <https://en.wikipedia.org/wiki/Portable_Network_Graphics>`_ is available
   in the current LAMMPS library.

   :to: :cpp:func:`lammps_config_has_png_support`
   :r has_png:
   :rtype has_png: logical

--------

.. f:function:: config_has_jpeg_support()

   Check if the LAMMPS library supports writing JPEG format images.

   .. versionadded:: TBD

   The LAMMPS :doc:`dump style image <dump_image>` supports writing multiple
   image file formats.  Most of them, however, need support from an external
   library, and using that has to be :ref:`enabled at compile time <graphics>`.
   This function checks whether support for the `JPEG image file format
   <https://jpeg.org/jpeg/>`_ is available in the current LAMMPS library.

   :to: :cpp:func:`lammps_config_has_jpeg_support`
   :r has_jpeg:
   :rtype has_jpeg: logical

--------

.. f:function:: config_has_ffmpeg_support()

   Check if the LAMMPS library supports creating movie files via a pipe to
   ffmpeg.

   .. versionadded:: TBD

   The LAMMPS :doc:`dump style movie <dump_image>` supports generating movies
   from images on-the-fly via creating a pipe to the
   `ffmpeg <https://ffmpeg.org/ffmpeg/>`_ program.
   This function checks whether this feature was
   :ref:`enabled at compile time <graphics>`.
   It does **not** check whether the ``ffmpeg`` itself is installed and usable.

   :to: :cpp:func:`lammps_config_has_ffmpeg_support`
   :r has_ffmpeg:
   :rtype has_ffmpeg: logical

--------

.. f:function:: config_has_exceptions()

   Check whether LAMMPS errors will throw C++ exceptions.

   .. versionadded:: TBD

   In case of an error, LAMMPS will either abort or throw a C++ exception.
   The latter has to be :ref:`enabled at compile time <exceptions>`.
   This function checks if exceptions were enabled.

   When using the library interface with C++ exceptions enabled, the library
   interface functions will "catch" them, and the error status can then be
   checked by calling :f:func:`has_error`. The most recent error message can be
   retrieved via :f:func:`get_last_error_message`.
   This can allow one to restart a calculation or delete and recreate
   the LAMMPS instance when a C++ exception occurs.  One application
   of using exceptions this way is the :ref:`lammps_shell`.  If C++
   exceptions are disabled and an error happens during a call to
   LAMMPS or the Fortran API, the application will terminate.

   :to: :cpp:func:`lammps_config_has_exceptions`
   :r has_exceptions:
   :rtype has_exceptions: logical

--------

.. f:function:: config_has_package(name)

   Check whether a specific package has been included in LAMMPS

   .. versionadded:: TBD

   This function checks whether the LAMMPS library in use includes the specific
   :doc:`LAMMPS package <Packages>` provided as argument.

   :to: :cpp:func:`lammps_config_has_package`
   :r has_package:
   :rtype has_package: logical

--------

.. f:function:: config_package_count()

   Count the number of installed packages in the LAMMPS library.

   .. versionadded:: TBD

   This function counts how many :doc:`LAMMPS packages <Packages>` are
   included in the LAMMPS library in use. It directly calls the C library
   function :cpp:func:`lammps_config_package_count`.

   :to: :cpp:func:`lammps_config_package_count`
   :r integer(c_int) npackages: number of packages installed

--------

.. f:subroutine:: config_package_name(idx, buffer)

   Get the name of a package in the list of installed packages in the LAMMPS
   library.

   .. versionadded:: TBD

   This subroutine copies the name of the package with the index *idx* into the
   provided string *buffer*. If the name of the package exceeds the length of
   the buffer, it will be truncated accordingly.  If the index is out of range,
   *buffer* is set to an empty string.

   :p integer(c_int) idx: index of the package in the list of included packages
    :math:`(0 \le idx < \text{package count})`
   :p character(len=\*) buffer: string to hold the name of the package
   :to: :cpp:func:`lammps_config_package_name`

--------

.. f:subroutine:: installed_packages(package[, length])

   Obtain a list of the names of enabled packages in the LAMMPS shared library
   and store it in *package*.

   .. versionadded:: TBD

   This function is analogous to the :py:func`installed_packages` function in
   the Python API. The optional argument *length* sets the length of each
   string in the vector *package* (default: 31).

   :p character(len=:) package [dimension(:),allocatable]: list of packages;
    *must* have the ``ALLOCATABLE`` attribute and be of rank 1
    (i.e., ``DIMENSION(:)``) with allocatable length.
   :o length: length of each string in the list. Default: 31.
   :otype length: integer,optional
   :to: :cpp:func:`lammps_config_package_count`
   :to: :cpp:func:`lammps_config_package_name`

--------

.. f:function:: config_accelerator(package, category, setting)

   This function calls :cpp:func:`lammps_config_accelerator` to check the
   availability of compile time settings of included
   :doc:`accelerator packages <Speed_packages>` in LAMMPS.

   .. versionadded:: TBD

   Supported packages names are "GPU", "KOKKOS", "INTEL", and "OPENMP".
   Supported categories are "api" with possible settings "cuda", "hip", "phi",
   "pthreads", "opencl", "openmp", and "serial"; and "precision" with
   possible settings "double", "mixed", and "single".

   :p character(len=\*) package:   string with the name of the accelerator
    package
   :p character(len=\*) category:  string with the name of the setting
   :p character(len=\*) setting:   string with the name of the specific
    setting
   :to: :cpp:func:`lammps_config_accelerator`
   :r available: ``.TRUE.`` if the combination of package, category,
    and setting is available, otherwise ``.FALSE.``.
   :rtype available: logical

--------

.. f:function:: has_gpu_device()

   Checks for the presence of a viable GPU package device.

   .. versionadded:: TBD

   This function calls :cpp:func:`lammps_has_gpu_device`, which checks at
   runtime whether an accelerator device is present that can be used with the
   :doc:`GPU package <Speed_gpu>`.

   More detailed information about the available device or devices can
   be obtained by calling the
   :f:subr:`get_gpu_device_info` subroutine.

   :to: :cpp:func:`lammps_has_gpu_device`
   :r available: ``.TRUE.`` if a viable device is available, ``.FALSE.`` if not.
   :rtype available: logical

--------

.. f:subroutine:: get_gpu_device_info(buffer)

   Get GPU package device information.

   .. versionadded:: TBD

   Calls :cpp:func:`lammps_get_gpu_device_info` to retrieve detailed
   information about any accelerator devices that are viable for use with the
   :doc:`GPU package <Speed_gpu>`. It will fill *buffer* with a string that is
   equivalent to the output of the ``nvc_get_device`` or ``ocl_get_device`` or
   ``hip_get_device`` tools that are compiled alongside LAMMPS if the GPU
   package is enabled.

   A suitable-length Fortran string has to be provided. The assembled text will
   be truncated so as not to overflow this buffer.  This string can be several
   kilobytes long if multiple devices are present.

   :p character(len=\*) buffer: string into which to copy the information.
   :to: :cpp:func:`lammps_get_gpu_device_info`

--------

.. f:function:: has_style(category, name)

   Check whether a specific style has been included in LAMMPS.

   .. versionadded:: TBD

   This function calls :cpp:func:`lammps_has_style` to check whether the
   LAMMPS library in use includes the specific style *name* associated with a
   specific *category* provided as arguments.  Please see
   :cpp:func:`lammps_has_style` for a list of valid categories.

   :p character(len=\*) category: category of the style
   :p character(len=\*) name:     name of the style
   :to: :cpp:func:`lammps_has_style`
   :r has_style: ``.TRUE.`` if included, ``.FALSE.`` if not.
   :rtype has_style: logical

--------

.. f:function:: style_count(category)

   Count the number of styles of *category* in the LAMMPS library.

   .. versionadded:: TBD

   This function counts how many styles in the provided *category* are
   included in the LAMMPS library currently in use. Please see
   :cpp:func:`lammps_has_style` for a list of valid categories.

   :p character(len=\*) category: category of styles to count
   :to: :cpp:func:`lammps_style_count`
   :r integer(c_int) count: number of styles in *category*

--------

.. f:subroutine:: style_name(category, idx, buffer)

   Look up the name of a style by index in the list of styles of a given
   category in the LAMMPS library.

   .. versionadded:: TBD

   This function calls :cpp:func:`lammps_style_name` and copies the name of
   the *category* style with index *idx* into the provided string *buffer*.
   The length of *buffer* must be long enough to contain the name of the
   style; if it is too short, the name will be truncated accordingly.
   If *idx* is out of range, *buffer* will be the empty string and a warning
   will be issued.

   :p character(len=\*) category: category of styles
   :p integer(c_int) idx:         index of the style in the list of *category*
    styles :math:`(1 \leq idx \leq \text{style count})`
   :p character(len\*) buffer:    string buffer to copy the name of the style
    into
   :to: :cpp:func:`lammps_style_name`

--------

.. f:function:: has_id(category, name)

   This function checks if the current LAMMPS instance a *category* ID of
   the given *name* exists.  Valid categories are: *compute*\ , *dump*\ ,
   *fix*\ , *group*\ , *molecule*\ , *region*\ , and *variable*\ .

   .. versionadded:: TBD

   :p character(len=\*) category: category of the ID
   :p character(len=\*) name:     name of the ID
   :to: :cpp:func:`lammps_has_id`
   :r has_id: ``.TRUE.`` if *category* style *name* exists, ``.FALSE.`` if not.
   :rtype has_id: logical

--------

.. f:function:: id_count(category)

   This function counts how many IDs in the provided *category* are defined in
   the current LAMMPS instance. Please see :f:func:`has_id` for a list of
   valid categories.

   .. versionadded:: TBD

   :p character(len=\*) category: category of the ID
   :to: :cpp:func:`lammps_id_count`
   :r count: number of IDs in *category*
   :rtype count: integer(c_int)

--------

.. f:subroutine:: id_name(category, idx, buffer)

   Look up the name of an ID by index in the list of IDs of a given category.

   .. versionadded:: TBD

   This function copies the name of the *category* ID with the index *idx* into
   the provided string *buffer*\ .  The length of the buffer must be long
   enough to hold the string; if the name of the style exceeds the length of
   the buffer, it will be truncated accordingly. If *buffer* is
   ``ALLOCATABLE``, it must be allocated *before* the function is called.
   If *idx* is out of range, *buffer* is set to an empty string and a warning
   is issued.

   :p character(len=\*) category: category of IDs
   :p integer(c_int) idx:         index of the ID in the list of *category*
    styles (:math:`0 \leq idx < count`)
   :p character(len=\*) buffer:   string into which to copy the name of the
    style
   :to: :cpp:func:`lammps_id_name`

--------

.. f:function:: plugin_count()

   This function counts the number of loaded plugins.

   .. versionadded:: TBD

   :to: :cpp:func:`lammps_plugin_count`
   :r n: number of loaded plugins
   :rtype n: integer(c_int)

--------

.. f:subroutine:: plugin_name(idx, stylebuf, namebuf)

   Look up the style and name of a plugin by its index in the list of plugins.

   .. versionadded:: TBD

   This function copies the name of the *style* plugin with the index *idx*
   into the provided C-style string buffer.  The length of the buffer must be
   provided as *buf_size* argument.  If the name of the style exceeds the
   length of the buffer, it will be truncated accordingly.  If the index is out
   of range, both strings are set to the empty string and a warning is
   printed.

   :p integer(c_int) idx:         index of the plugin in the list all or
    *style* plugins
   :p character(len=\*) stylebuf: string into which to copy the style of the
    plugin
   :p character(len=\*) namebuf:  string into which to copy the style of the
    plugin
   :to: :cpp:func:`lammps_plugin_name`

--------

.. f:function:: encode_image_flags(ix, iy, iz)

   Encodes three integer image flags into a single imageint.

   .. versionadded:: TBD

   This function performs the bit-shift, addition, and bit-wise OR operations
   necessary to combine the values of three integers representing the image
   flags in the :math:`x`-, :math:`y`-, and :math:`z`-directions. Unless LAMMPS
   is compiled with ``-DLAMMPS_BIGBIG``, those integers are limited to 10-bit
   signed integers :math:`[-512,512]`. If ``-DLAMMPS_BIGBIG`` was used when
   compiling, then the return value is of kind ``c_int64_t`` instead of
   kind ``c_int``, and the valid range for the individual image flags becomes
   :math:`[-1048576,1048575]` (i.e., the range of a 21-bit signed integer).
   There is no check on whether the arguments conform to these requirements.

   :p integer(c_int) ix: image flag in :math:`x`-direction
   :p integer(c_int) iy: image flag in :math:`y`-direction
   :p integer(c_int) iz: image flag in :math:`z`-direction
   :r integer(kind=\*) imageint: encoded image flag. \*The ``KIND`` parameter
    is ``c_int`` unless LAMMPS was built with ``-DLAMMPS_BIGBIG``, in which
    case it is ``c_int64_t``.

   .. note::

     The fact that the programmer does not know the ``KIND`` parameter of the
     return value until compile time means that it is impossible to define an
     interface that works for both sizes of ``imageint``. One side effect of
     this is that you must assign the return value of this function to a
     variable; it cannot be used as the argument to another function or as part
     of an array constructor. For example,

     .. code-block:: fortran

        my_images = [lmp%encode_image_flags(0,0,0), lmp%encode_image_flags(1,0,0)]

     will *not* work; instead, do something like

     .. code-block:: fortran

        my_images(1) = lmp%encode_image_flags(0,0,0)
        my_images(2) = lmp%encode_image_flags(1,0,0)

--------

.. f:subroutine:: decode_image_flags(image, flags)

   This function does the reverse operation of :f:func:`encode_image_flags`:
   it takes the image flag and performs the bit-shift and bit-masking
   operations to decode it and stores the resulting three integers into the
   array *flags*.

   .. versionadded:: TBD

   :p integer(kind=\*) image: encoded image flag. \*The ``KIND`` paremeter is
    either ``c_int`` or, if LAMMPS was compiled with ``-DLAMMPS_BIGBIG``,
    ``c_int64_t``. Kind compatibility is checked at run-time.
   :p integer(c_int) flags [dimension(3)]: three-element vector where the
    decoded image flags will be stored.

--------

.. f:subroutine:: flush_buffers()

   This function calls :cpp:func:`lammps_flush_buffers`, which flushes buffered
   output to be written to screen and logfile. This can simplify capturing
   output from LAMMPS library calls.

   .. versionadded:: TBD

   :to: :cpp:func:`lammps_flush_buffers`

--------

.. f:function:: is_running()

   Check if LAMMPS is currently inside a run or minimization.

   .. versionadded:: TBD

   This function can be used from signal handlers or multi-threaded
   applications to determine if the LAMMPS instance is currently active.

   :to: :cpp:func:`lammps_is_running`
   :r is_running: ``.FALSE.`` if idle or ``.TRUE.`` if active
   :rtype is_running: logical

--------

.. f:subroutine:: force_timeout()

   Force a timeout to stop an ongoing run cleanly.

   .. versionadded:: TBD

   This function can be used from signal handlers or multi-threaded
   applications to cleanly terminate an ongoing run.

   :to: :cpp:func:`lammps_force_timeout`

--------

.. f:function:: has_error()

   Check if there is a (new) error message available.

   .. versionadded:: TBD

   This function can be used to query if an error inside of LAMMPS
   has thrown a :ref:`C++ exception <exceptions>`.

   .. note::

      This function will always report "no error" when the LAMMPS library
      has been compiled without ``-DLAMMPS_EXCEPTIONS``, which turns fatal
      errors aborting LAMMPS into C++ exceptions. You can use the library
      function :cpp:func:`lammps_config_has_exceptions` to check if this is
      the case.

   :to: :cpp:func:`lammps_has_error`
   :r has_error: ``.TRUE.`` if there is an error.
   :rtype has_error: logical

--------

.. f:subroutine:: get_last_error_message(buffer[,status])

   Copy the last error message into the provided buffer.

   .. versionadded:: TBD

   This function can be used to retrieve the error message that was set
   in the event of an error inside of LAMMPS that resulted in a
   :ref:`C++ exception <exceptions>`.  A suitable buffer for a string has
   to be provided.  If the internally-stored error message is longer than the
   string and the string does not have ``ALLOCATABLE`` length, it will be
   truncated accordingly.  The optional argument *status* indicates the
   kind of error: a "1" indicates an error that occurred on all MPI ranks and
   is often recoverable, while a "2" indicates an abort that would happen only
   in a single MPI rank and thus may not be recoverable, as other MPI ranks may
   be waiting on the failing MPI rank(s) to send messages.

   .. note::

      This function will do nothing when the LAMMPS library has been
      compiled without ``-DLAMMPS_EXCEPTIONS``, which turns errors aborting
      LAMMPS into C++ exceptions.  You can use the function
      :f:func:`config_has_exceptions` to check whether this is the case.

   :p character(len=\*) buffer: string buffer to copy the error message into
   :o integer(c_int) status [optional]: 1 when all ranks had the error,
    2 on a single-rank error.
   :to: :cpp:func:`lammps_get_last_error_message`

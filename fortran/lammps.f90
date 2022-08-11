! -------------------------------------------------------------------------
!   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
!   https://www.lammps.org/ Sandia National Laboratories
!   Steve Plimpton, sjplimp@sandia.gov
!
!   Copyright (2003) Sandia Corporation.  Under the terms of Contract
!   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
!   certain rights in this software.  This software is distributed under
!   the GNU General Public License.
!
!   See the README file in the top-level LAMMPS directory.
! -------------------------------------------------------------------------
!
! Fortran interface to the LAMMPS library implemented as a Fortran 2003
! style module that wraps the C-style library interface in library.cpp
! and library.h using the ISO_C_BINDING module of the Fortran compiler.
!
! Based on the LAMMPS Fortran 2003 module contributed by:
!   Karl D. Hammond <hammondkd@missouri.edu>
!   University of Missouri, 2012-2020
!
! The Fortran module tries to follow the API of the C-library interface
! closely, but like the Python wrapper it employs an object oriented
! approach.  To accommodate the object oriented approach, all exported
! subroutine and functions have to be implemented in Fortran to then
! call the interfaced C style functions with adapted calling conventions
! as needed.  The C-library interfaced functions retain their names
! starting with "lammps_" while the Fortran versions start with "lmp_".
!
MODULE LIBLAMMPS

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_null_ptr, c_loc, &
      c_int, c_int64_t, c_char, c_null_char, c_double, c_size_t, c_f_pointer

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lammps

!  ! These are public-interface constants that have the same purpose as the
!  ! constants in library.h, except that their types match the type of the
!  ! constant in question. Their purpose is to specify the type of the
!  ! return value without something akin to a C/C++ type cast
!  INTEGER (c_int), PUBLIC, PARAMETER :: LMP_INT = 0_c_int
!  INTEGER (c_int), PUBLIC, DIMENSION(3), PARAMETER :: LMP_INT_1D = 0_c_int
!  INTEGER (c_int), PUBLIC, DIMENSION(3,3), PARAMETER :: LMP_INT_2D = 1_c_int
!  REAL (c_double), PUBLIC, PARAMETER :: LMP_DOUBLE = 2.0_c_double
!  REAL (c_double), PUBLIC, DIMENSION(3), PARAMETER :: &
!    LMP_DOUBLE_1D = 2.0_c_double
!  REAL (c_double), PUBLIC, DIMENSION(3,3), PARAMETER :: &
!    LMP_DOUBLE_2D = 3.0_c_double
!  INTEGER (c_int64_t), PUBLIC, PARAMETER :: LMP_INT64 = 4_c_int64_t
!  INTEGER (c_int64_t), PUBLIC, DIMENSION(3), PARAMETER :: &
!    LMP_INT64_1D = 4_c_int64_t
!  INTEGER (c_int64_t), PUBLIC, DIMENSION(3,3), PARAMETER :: &
!    LMP_INT64_2D = 5_c_int64_t
!  CHARACTER(LEN=*), PUBLIC, PARAMETER :: LMP_STRING = 'six'

  ! Data type constants for extracting data from global, atom, compute, and fix
  !
  ! Must be kept in sync with the equivalent declarations in
  ! src/library.h and python/lammps/constants.py
  !
  ! NOT part of the API (the part the user sees)
  INTEGER (c_int), PARAMETER :: &
    LAMMPS_INT = 0_c_int, &    ! 32-bit integer (array)
    LAMMPS_INT_2D = 1, &    ! two-dimensional 32-bit integer array
    LAMMPS_DOUBLE = 2, &    ! 64-bit double (array)
    LAMMPS_DOUBLE_2D = 3, & ! two-dimensional 64-bit double array
    LAMMPS_INT64 = 4, &     ! 64-bit integer (array)
    LAMMPS_INT64_2D = 5, &  ! two-dimensional 64-bit integer array
    LAMMPS_STRING = 6       ! C-String

  TYPE lammps
      TYPE(c_ptr) :: handle
    CONTAINS
      PROCEDURE :: close              => lmp_close
      PROCEDURE :: file               => lmp_file
      PROCEDURE :: command            => lmp_command
      PROCEDURE :: commands_list      => lmp_commands_list
      PROCEDURE :: commands_string    => lmp_commands_string
      PROCEDURE :: get_natoms         => lmp_get_natoms
      PROCEDURE :: get_thermo         => lmp_get_thermo
      PROCEDURE :: extract_box        => lmp_extract_box
      PROCEDURE :: reset_box          => lmp_reset_box
      PROCEDURE :: memory_usage       => lmp_memory_usage
      PROCEDURE :: get_mpi_comm       => lmp_get_mpi_comm
      PROCEDURE :: extract_setting    => lmp_extract_setting
!      PROCEDURE :: extract_global     => lmp_extract_global
!      PROCEDURE, PRIVATE :: lmp_extract_global_int
!      PROCEDURE, PRIVATE :: lmp_extract_global_int64_t
!      PROCEDURE, PRIVATE :: lmp_extract_global_double
!      PROCEDURE, PRIVATE :: lmp_extract_global_str
!      GENERIC :: extract_global       => lmp_extract_global_int, &
!        lmp_extract_global_int64_t, lmp_extract_global_double, &
!        lmp_extract_global_str
!      PROCEDURE, PRIVATE :: lmp_extract_global_scalar
!      !PROCEDURE, PRIVATE :: lmp_extract_global_string
!      GENERIC :: extract_global       => lmp_extract_global_scalar
      PROCEDURE :: version            => lmp_version
  END TYPE lammps

  INTERFACE lammps
      MODULE PROCEDURE lmp_open
  END INTERFACE lammps

  ! interface definitions for calling functions in library.cpp
  INTERFACE
      FUNCTION lammps_open(argc, argv, comm) BIND(C, name='lammps_open_fortran')
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(IN)     :: argc, comm
        TYPE(c_ptr), DIMENSION(*), INTENT(IN) :: argv
        TYPE(c_ptr)                           :: lammps_open
      END FUNCTION lammps_open

      FUNCTION lammps_open_no_mpi(argc, argv, handle) BIND(C)
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(IN)     :: argc
        TYPE(c_ptr), DIMENSION(*), INTENT(IN) :: argv
        TYPE(c_ptr), VALUE, INTENT(in)        :: handle
        TYPE(c_ptr)                           :: lammps_open_no_mpi
      END FUNCTION lammps_open_no_mpi

      SUBROUTINE lammps_close(handle) BIND(C)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
      END SUBROUTINE lammps_close

      SUBROUTINE lammps_mpi_init() BIND(C)
      END SUBROUTINE lammps_mpi_init

      SUBROUTINE lammps_mpi_finalize() BIND(C)
      END SUBROUTINE lammps_mpi_finalize

      SUBROUTINE lammps_kokkos_finalize() BIND(C)
      END SUBROUTINE lammps_kokkos_finalize

      SUBROUTINE lammps_file(handle, filename) BIND(C)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: filename
      END SUBROUTINE lammps_file

      SUBROUTINE lammps_command(handle, cmd) BIND(C)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: cmd
      END SUBROUTINE lammps_command

      SUBROUTINE lammps_commands_list(handle, ncmd, cmds) BIND(C)
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE :: handle
        INTEGER(c_int), VALUE, INTENT(IN)     :: ncmd
        TYPE(c_ptr), DIMENSION(*), INTENT(IN) :: cmds
      END SUBROUTINE lammps_commands_list

      SUBROUTINE lammps_commands_string(handle, str) BIND(C)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: str
      END SUBROUTINE lammps_commands_string

      FUNCTION lammps_get_natoms(handle) BIND(C)
        IMPORT :: c_ptr, c_double
        TYPE(c_ptr), VALUE :: handle
        REAL(c_double) :: lammps_get_natoms
      END FUNCTION lammps_get_natoms

      FUNCTION lammps_get_thermo(handle,name) BIND(C)
        IMPORT :: c_ptr, c_double
        IMPLICIT NONE
        REAL(c_double) :: lammps_get_thermo
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: name
      END FUNCTION lammps_get_thermo

      SUBROUTINE lammps_extract_box(handle,boxlo,boxhi,xy,yz,xz,pflags, &
            boxflag) BIND(C)
         IMPORT :: c_ptr, c_double, c_int
         IMPLICIT NONE
         TYPE(c_ptr), VALUE :: handle, boxlo, boxhi, xy, yz, xz, pflags, &
            boxflag
      END SUBROUTINE lammps_extract_box

      SUBROUTINE lammps_reset_box(handle,boxlo,boxhi,xy,yz,xz) BIND(C)
        IMPORT :: c_ptr, c_double
        IMPLICIT NONE
        TYPE (c_ptr), VALUE :: handle
        REAL (c_double), DIMENSION(3) :: boxlo, boxhi
        REAL (c_double), VALUE :: xy, yz, xz
      END SUBROUTINE lammps_reset_box

      SUBROUTINE lammps_memory_usage(handle,meminfo) BIND(C)
        IMPORT :: c_ptr, c_double
        IMPLICIT NONE
        TYPE(c_ptr), VALUE :: handle
        REAL(c_double), DIMENSION(*) :: meminfo
      END SUBROUTINE lammps_memory_usage

      FUNCTION lammps_get_mpi_comm(handle) BIND(C)
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE (c_ptr), VALUE :: handle
        INTEGER (c_int) :: lammps_get_mpi_comm
      END FUNCTION lammps_get_mpi_comm

      FUNCTION lammps_extract_setting(handle,keyword) BIND(C)
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE :: handle, keyword
        INTEGER (c_int) :: lammps_extract_setting
      END FUNCTION lammps_extract_setting

      FUNCTION lammps_extract_global_datatype(handle,name) BIND(C)
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE :: handle, name
        INTEGER (c_int) :: lammps_extract_global_datatype
      END FUNCTION lammps_extract_global_datatype

      FUNCTION c_strlen (str) bind(C,name='strlen')
        IMPORT :: c_ptr, c_size_t
        IMPLICIT NONE
        TYPE(c_ptr) :: str
        INTEGER(c_size_t) :: c_strlen
      END FUNCTION c_strlen

      FUNCTION lammps_extract_global(handle, name) BIND(C)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle, name
        TYPE(c_ptr) :: lammps_extract_global
      END FUNCTION lammps_extract_global

      !INTEGER (c_int) FUNCTION lammps_extract_atom_datatype

      !(generic) lammps_extract_atom

      !(generic) lammps_extract_compute

      !(generic) lammps_extract_fix

      !(generic) lammps_extract_variable

      !INTEGER (c_int) lammps_set_variable

      !SUBROUTINE lammps_gather_atoms

      !SUBROUTINE lammps_gather_atoms_concat

      !SUBROUTINE lammps_gather_atoms_subset

      !SUBROUTINE lammps_scatter_atoms

      !SUBROUTINE lammps_scatter_atoms_subset

      !SUBROUTINE lammps_gather_bonds

      !SUBROUTINE lammps_gather

      !SUBROUTINE lammps_gather_concat

      !SUBROUTINE lammps_gather_subset

      !SUBROUTINE lammps_scatter_subset

      !(generic / id, type, and image are special) / requires LAMMPS_BIGBIG
      !INTEGER (C_int) FUNCTION lammps_create_atoms

      !INTEGER (C_int) FUNCTION lammps_find_pair_neighlist

      !INTEGER (C_int) FUNCTION lammps_find_fix_neighlist

      !INTEGER (C_int) FUNCTION lammps_find_compute_neighlist

      !INTEGER (C_int) FUNCTION lammps_neighlist_num_elements

      !SUBROUTINE lammps_neighlist_element_neighbors

      FUNCTION lammps_version(handle) BIND(C)
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE :: handle
        INTEGER(c_int) :: lammps_version
      END FUNCTION lammps_version

      !SUBROUTINE lammps_get_os_info

      !LOGICAL FUNCTION lammps_config_has_mpi_support
      !LOGICAL FUNCTION lammps_config_has_gzip_support
      !LOGICAL FUNCTION lammps_config_has_png_support
      !LOGICAL FUNCTION lammps_config_has_jpeg_support
      !LOGICAL FUNCTION lammps_config_has_ffmpeg_support
      !LOGICAL FUNCTION lammps_config_has_exceptions
      !LOGICAL FUNCTION lammps_config_has_package
      !INTEGER (C_int) FUNCTION lammps_config_package_count
      !SUBROUTINE lammps_config_package_name

      !LOGICAL FUNCTION lammps_config_accelerator
      !LOGICAL FUNCTION lammps_has_gpu_device
      !SUBROUTINE lammps_get_gpu_device

      !LOGICAL FUNCTION lammps_has_id
      !INTEGER (C_int) FUNCTION lammps_id_count
      !SUBROUTINE lammps_id_name

      !INTEGER (C_int) FUNCTION lammps_plugin_count
      !SUBROUTINE lammps_plugin_name

      !Both of these use LAMMPS_BIGBIG
      !INTEGER (LAMMPS_imageint) FUNCTION lammps_encode_image_flags
      !SUBROUTINE lammps_decode_image_flags

      !SUBROUTINE lammps_set_fix_external_callback ! may have trouble....
      !FUNCTION lammps_fix_external_get_force() ! returns real(c_double) (:)

      !SUBROUTINE lammps_fix_external_set_energy_global
      !SUBROUTINE lammps_fix_external_set_energy_peratom
      !SUBROUTINE lammps_fix_external_set_virial_global
      !SUBROUTINE lammps_fix_external_set_virial_peratom
      !SUBROUTINE lammps_fix_external_set_vector_length
      !SUBROUTINE lammps_fix_external_set_vector

      !SUBROUTINE lammps_flush_buffers

      FUNCTION lammps_malloc(size) BIND(C, name='malloc')
        IMPORT :: c_ptr, c_size_t
        INTEGER(c_size_t), VALUE :: size
        TYPE(c_ptr) :: lammps_malloc
      END FUNCTION lammps_malloc

      SUBROUTINE lammps_free(ptr) BIND(C)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: ptr
      END SUBROUTINE lammps_free


      !LOGICAL FUNCTION lammps_is_running
      !SUBROUTINE lammps_force_timeout
      !LOGICAL FUNCTION lammps_has_error
      !INTEGER (c_int) FUNCTION lammps_get_last_error_message
  END INTERFACE

CONTAINS
  ! Fortran wrappers and helper functions.

  ! Constructor for the LAMMPS class.
  ! Combined wrapper around lammps_open_fortran() and lammps_open_no_mpi()
  TYPE(lammps) FUNCTION lmp_open(args, comm)
    IMPLICIT NONE
    INTEGER, INTENT(in), OPTIONAL :: comm
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: args(:)
    TYPE(c_ptr), ALLOCATABLE     :: argv(:)
    INTEGER(c_int)               :: i, c_comm, argc

    IF (PRESENT(args)) THEN
        ! convert fortran argument list to c style
        argc = SIZE(args)
        ALLOCATE(argv(argc))
        DO i=1, argc
           argv(i) = f2c_string(args(i))
        END DO
    ELSE
        argc = 1
        ALLOCATE(argv(1))
        argv(1) = f2c_string("liblammps")
    ENDIF

    IF (PRESENT(comm)) THEN
        c_comm = comm
        lmp_open%handle = lammps_open(argc, argv, c_comm)
    ELSE
        lmp_open%handle = lammps_open_no_mpi(argc, argv, c_null_ptr)
    END IF

    ! Clean up allocated memory
    DO i=1, argc
        CALL lammps_free(argv(i))
    END DO
    DEALLOCATE(argv)
  END FUNCTION lmp_open

  ! Combined Fortran wrapper around lammps_close() and lammps_mpi_finalize()
  SUBROUTINE lmp_close(self, finalize)
    IMPLICIT NONE
    CLASS(lammps) :: self
    LOGICAL, INTENT(IN), OPTIONAL :: finalize

    CALL lammps_close(self%handle)

    IF (PRESENT(finalize)) THEN
        IF (finalize) THEN
            CALL lammps_kokkos_finalize()
            CALL lammps_mpi_finalize()
        END IF
    END IF
  END SUBROUTINE lmp_close

  SUBROUTINE lmp_file(self, filename)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*) :: filename
    TYPE(c_ptr) :: str

    str = f2c_string(filename)
    CALL lammps_file(self%handle, str)
    CALL lammps_free(str)
  END SUBROUTINE lmp_file

  ! equivalent function to lammps_command()
  SUBROUTINE lmp_command(self, cmd)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*) :: cmd
    TYPE(c_ptr) :: str

    str = f2c_string(cmd)
    CALL lammps_command(self%handle, str)
    CALL lammps_free(str)
  END SUBROUTINE lmp_command

  ! equivalent function to lammps_commands_list()
  SUBROUTINE lmp_commands_list(self, cmds)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cmds(:)
    TYPE(c_ptr), ALLOCATABLE     :: cmdv(:)
    INTEGER :: i, ncmd

    ! convert command list to c style
    ncmd = SIZE(cmds)
    ALLOCATE(cmdv(ncmd))
    DO i=1, ncmd
        cmdv(i) = f2c_string(cmds(i))
    END DO

    CALL lammps_commands_list(self%handle, ncmd, cmdv)

    ! Clean up allocated memory
    DO i=1, ncmd
        CALL lammps_free(cmdv(i))
    END DO
    DEALLOCATE(cmdv)
  END SUBROUTINE lmp_commands_list

  ! equivalent function to lammps_commands_string()
  SUBROUTINE lmp_commands_string(self, str)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*) :: str
    TYPE(c_ptr) :: tmp

    tmp = f2c_string(str)
    CALL lammps_commands_string(self%handle, tmp)
    CALL lammps_free(tmp)
  END SUBROUTINE lmp_commands_string

  ! equivalent function to lammps_get_natoms
  DOUBLE PRECISION FUNCTION lmp_get_natoms(self)
    IMPLICIT NONE
    CLASS(lammps) :: self

    lmp_get_natoms = lammps_get_natoms(self%handle)
  END FUNCTION lmp_get_natoms

  ! equivalent function to lammps_get_thermo
  REAL (C_double) FUNCTION lmp_get_thermo(self,name)
    CLASS(lammps), INTENT(IN) :: self
    CHARACTER(LEN=*) :: name
    TYPE(C_ptr) :: Cname

    Cname = f2c_string(name)
    lmp_get_thermo = lammps_get_thermo(self%handle, Cname)
    CALL lammps_free(Cname)
  END FUNCTION lmp_get_thermo

  ! equivalent subroutine to lammps_extract_box
  SUBROUTINE lmp_extract_box (self, boxlo, boxhi, xy, yz, xz, pflags, boxflag)
    CLASS(lammps), INTENT(IN) :: self
    REAL(c_double), INTENT(OUT), TARGET, OPTIONAL :: boxlo(3), boxhi(3)
    REAL(c_double), INTENT(OUT), TARGET, OPTIONAL :: xy, yz, xz
    LOGICAL, INTENT(OUT), OPTIONAL :: pflags(3), boxflag
    INTEGER(c_int), TARGET :: C_pflags(3), C_boxflag
    TYPE (c_ptr) :: ptr(7)

    ptr = c_null_ptr
    IF ( PRESENT(boxlo) ) ptr(1) = C_LOC(boxlo(1))
    IF ( PRESENT(boxhi) ) ptr(2) = C_LOC(boxhi(1))
    IF ( PRESENT(xy) ) ptr(3) = C_LOC(xy)
    IF ( PRESENT(yz) ) ptr(4) = C_LOC(yz)
    IF ( PRESENT(xz) ) ptr(5) = C_LOC(xz)
    IF ( PRESENT(pflags) ) ptr(6) = C_LOC(C_pflags(1))
    IF ( PRESENT(boxflag) ) ptr(7) = C_LOC(C_boxflag)
    CALL lammps_extract_box(self%handle, ptr(1), ptr(2), ptr(3), ptr(4), &
      ptr(5), ptr(6), ptr(7))
    IF ( PRESENT(pflags) ) pflags = ( C_pflags /= 0_C_int )
    IF ( PRESENT(boxflag) ) boxflag = ( C_boxflag /= 0_C_int )
  END SUBROUTINE lmp_extract_box

  ! equivalent function to lammps_reset_box
  SUBROUTINE lmp_reset_box (self, boxlo, boxhi, xy, yz, xz)
    CLASS(lammps), INTENT(IN) :: self
    REAL(C_double), INTENT(IN) :: boxlo(3), boxhi(3), xy, yz, xz

    CALL lammps_reset_box (self%handle, boxlo, boxhi, xy, yz, xz)
  END SUBROUTINE lmp_reset_box

  ! equivalent function to lammps_memory_usage
  SUBROUTINE lmp_memory_usage(self,meminfo)
    CLASS(lammps), INTENT(IN) :: self
    INTEGER, PARAMETER :: MEMINFO_ELEM = 3
    REAL (c_double), DIMENSION(MEMINFO_ELEM), INTENT(OUT) :: meminfo

    CALL lammps_memory_usage(self%handle,meminfo)
  END SUBROUTINE lmp_memory_usage

  ! equivalent function to lammps_get_mpi_comm
  INTEGER FUNCTION lmp_get_mpi_comm (self)
    CLASS(lammps), INTENT(IN) :: self

    lmp_get_mpi_comm = lammps_get_mpi_comm(self%handle)
  END FUNCTION lmp_get_mpi_comm

  ! equivalent function to lammps_extract_setting
  INTEGER (c_int) FUNCTION lmp_extract_setting (self, keyword)
    CLASS(lammps), INTENT(IN) :: self
    CHARACTER(LEN=*), INTENT(IN) :: keyword
    TYPE(c_ptr) :: Ckeyword

    Ckeyword = f2c_string(keyword)
    lmp_extract_setting = lammps_extract_setting(self%handle, Ckeyword)
    CALL lammps_free(Ckeyword)
  END FUNCTION lmp_extract_setting

! FIXME Now that I think about it...do we need this at all?
!  ! equivalent function to lammps_extract_global_datatype
!  ! this function doesn't need to be accessed by the user, but is instead used
!  ! for type checking
!  INTEGER (c_int) FUNCTION lmp_extract_global_datatype (name)
!    CHARACTER(LEN=*), INTENT(IN) :: name
!    TYPE(c_ptr) :: Cname
!
!    Cname = f2c_string(name)
!    lmp_extract_global_datatype 
!      = lammps_extract_global_datatype(c_null_ptr, Cname)
!    CALL lammps_free(Cname)
!  END FUNCTION lmp_extract_global_datatype

  ! equivalent function to lammps_extract_global
  ! the return value should be automatically returned and assigned correctly
  ! based on the information available from LAMMPS
!  SUBROUTINE lmp_extract_global_scalar (self, global_data, name)
!    CLASS(lammps), INTENT(IN) :: self
!    CLASS(*), INTENT(OUT), POINTER :: global_data
!    CHARACTER(LEN=*), INTENT(IN) :: name
!    INTEGER(c_int) :: datatype
!    TYPE(c_ptr) :: Cname, Cptr
!    INTEGER(c_size_t) :: length, i
!    CHARACTER(KIND=c_char, LEN=1), DIMENSION(:), POINTER :: Fptr
!
!    Cname = f2c_string(name)
!    datatype = lammps_extract_global_datatype(self%handle, Cname)
!      ! above could be c_null_ptr in place of self%handle...doesn't matter
!    Cptr = lammps_extract_global(self%handle, Cname)
!    SELECT CASE (datatype)
!      CASE (LAMMPS_INT)
!        SELECT TYPE (global_data)
!          TYPE IS (INTEGER(c_int))
!            CALL C_F_POINTER(Cptr, global_data)
!          CLASS DEFAULT
!            ! FIXME
!            WRITE(0,'(A)') 'ERROR: Incompatible pointer type in extract_global'
!            STOP
!        END SELECT
!      CASE (LAMMPS_INT64)
!        SELECT TYPE (global_data)
!          TYPE IS (INTEGER(c_int64_t))
!            CALL C_F_POINTER(Cptr, global_data)
!          CLASS DEFAULT
!            ! FIXME
!            WRITE(0,'(A)') 'ERROR: Incompatible pointer type in extract_global'
!            STOP
!        END SELECT
!      CASE (LAMMPS_DOUBLE)
!        SELECT TYPE (global_data)
!          TYPE IS (REAL(c_double))
!            CALL C_F_POINTER(Cptr, global_data)
!          CLASS DEFAULT
!            ! FIXME
!            WRITE(0,'(A)') 'ERROR: Incompatible pointer type in extract_global'
!            STOP
!        END SELECT
!      CASE (LAMMPS_STRING)
!        SELECT TYPE (global_data)
!          TYPE IS (CHARACTER(LEN=*))
!            length = c_strlen(Cptr)
!            CALL C_F_POINTER(Cptr, Fptr, [length])
!            IF ( length < len(global_data) ) length = len(global_data)
!            FORALL ( i = 1:length )
!              global_data(i:i) = Fptr(i)
!            END FORALL
!        END SELECT
!      CASE DEFAULT
!        ! FIXME
!        WRITE (0,'(A,1X,I0,1X,A)') 'ERROR: Unknown type', datatype, &
!          'returned from extract_global_datatype'
!        STOP
!    END SELECT
!    CALL lammps_free(Cname)
!  END SUBROUTINE lmp_extract_global_scalar
!
!  SUBROUTINE lmp_extract_global_string (self, global_data, name)
!    CLASS(lammps), INTENT(IN) :: self
!    CHARACTER(LEN=*), INTENT(OUT) :: global_data
!    CHARACTER(LEN=*), INTENT(IN) :: name
!    INTEGER(c_int) :: datatype
!    TYPE(c_ptr) :: Cname, Cptr
!    CHARACTER(KIND=c_char, LEN=1), DIMENSION(:), POINTER :: Fptr
!    INTEGER(c_size_t) :: length
!    INTEGER :: i
!
!    global_data = ''
!    Cname = f2c_string(name)
!    datatype = lammps_extract_global_datatype(self%handle, Cname)
!    IF ( datatype /= LAMMPS_STRING ) THEN
!      ! FIXME
!      WRITE (0,'(A)') 'ERROR: Cannot assign string to non-string variable.'
!      STOP
!    END IF
!    Cptr = lammps_extract_global(self%handle, Cname)
!    length = c_strlen(Cptr)
!    CALL C_F_POINTER(Cptr, Fptr, [length])
!    IF ( length < len(global_data) ) length = len(global_data)
!    FORALL ( i = 1:length )
!      global_data(i:i) = Fptr(i)
!    END FORALL
!    CALL lammps_free(Cname)
!  END SUBROUTINE lmp_extract_global_string

!  ! equivalent functions to lammps_extract_global (overloaded)
!  ! This implementation assumes there are no non-scalar data that can be
!  ! extracted through lammps_extract_global
!  FUNCTION lmp_extract_global_int (self, name, dtype)
!    CLASS(lammps), INTENT(IN) :: self
!    CHARACTER(LEN=*), INTENT(IN) :: name
!    INTEGER(c_int), INTENT(IN) :: dtype
!    INTEGER(c_int) :: lmp_extract_global_int
!    TYPE(c_ptr) :: Cname, Cptr
!    INTEGER(c_int) :: datatype
!    INTEGER(c_int), POINTER :: ptr
!
!    Cname = f2c_string(name)
!    datatype = lammps_extract_global_datatype(c_null_ptr, Cname)
!    IF ( datatype /= LAMMPS_INT ) THEN
!      ! throw an exception or something; data type doesn't match!
!      WRITE(0,*) 'WARNING: global data type is inconsistent (not an int)'
!    END IF
!    Cptr = lammps_extract_global(self%handle, Cname)
!    CALL C_F_POINTER(Cptr, ptr)
!    lmp_extract_global_int = ptr
!    CALL lammps_free(Cname)
!  END FUNCTION lmp_extract_global_int
!  FUNCTION lmp_extract_global_int64_t (self, name, dtype)
!    CLASS(lammps), INTENT(IN) :: self
!    CHARACTER(LEN=*), INTENT(IN) :: name
!    INTEGER(c_int64_t), INTENT(IN) :: dtype
!    INTEGER(c_int64_t) :: lmp_extract_global_int64_t
!    TYPE(c_ptr) :: Cname, Cptr
!    INTEGER(c_int) :: datatype
!    INTEGER(c_int64_t), POINTER :: ptr
!
!    Cname = f2c_string(name)
!    datatype = lammps_extract_global_datatype(c_null_ptr, Cname)
!    IF ( datatype /= LAMMPS_INT64 ) THEN
!      ! throw an exception or something; data type doesn't match!
!      WRITE(0,*) 'WARNING: global data type is inconsistent (not an int64_t)'
!    END IF
!    Cptr = lammps_extract_global(self%handle, Cname)
!    CALL C_F_POINTER(Cptr, ptr)
!    lmp_extract_global_int64_t = ptr
!    CALL lammps_free(Cname)
!  END FUNCTION lmp_extract_global_int64_t
!  FUNCTION lmp_extract_global_double (self, name, dtype)
!    CLASS(lammps), INTENT(IN) :: self
!    CHARACTER(LEN=*), INTENT(IN) :: name
!    REAL(c_double), INTENT(IN) :: dtype
!    REAL(c_double) :: lmp_extract_global_double
!    TYPE(c_ptr) :: Cname, Cptr
!    INTEGER(c_int) :: datatype
!    REAL(c_double), POINTER :: ptr
!
!    Cname = f2c_string(name)
!    datatype = lammps_extract_global_datatype(c_null_ptr, Cname)
!    IF ( datatype /= LAMMPS_DOUBLE ) THEN
!      ! throw an exception or something; data type doesn't match!
!      WRITE(0,*) 'WARNING: global data type is inconsistent (not a double)'
!    END IF
!    Cptr = lammps_extract_global(self%handle, Cname)
!    CALL C_F_POINTER(Cptr, ptr)
!    lmp_extract_global_double = ptr
!    CALL lammps_free(Cname)
!  END FUNCTION lmp_extract_global_double
!  FUNCTION lmp_extract_global_str (self, name, dtype)
!    CLASS(lammps), INTENT(IN) :: self
!    CHARACTER(LEN=*), INTENT(IN) :: name, dtype
!    CHARACTER(LEN=:), ALLOCATABLE :: lmp_extract_global_str
!    TYPE(c_ptr) :: Cname, Cptr
!    INTEGER(c_int) :: datatype
!    CHARACTER(KIND=c_char,LEN=1), dimension(:), POINTER :: ptr
!    INTEGER(c_size_t) :: length
!    INTEGER :: i
!
!    Cname = f2c_string(name)
!    datatype = lammps_extract_global_datatype(c_null_ptr, Cname)
!    IF ( datatype /= LAMMPS_STRING ) THEN
!      ! throw an exception or something; data type doesn't match!
!      WRITE(0,*) 'WARNING: global data type is inconsistent (not a string)'
!    END IF
!    Cptr = lammps_extract_global(self%handle, Cname)
!    length = c_strlen(Cptr)
!    CALL C_F_POINTER(Cptr, ptr, [length])
!    ALLOCATE ( CHARACTER(LEN=length) :: lmp_extract_global_str )
!    FORALL ( I=1:length )
!      lmp_extract_global_str(i:i) = ptr(i)
!    END FORALL
!    CALL lammps_free(Cname)
!    ! the allocatable scalar (return value) gets auto-deallocated on return
!  END FUNCTION lmp_extract_global_str

  ! equivalent function to lammps_version()
  INTEGER FUNCTION lmp_version(self)
    IMPLICIT NONE
    CLASS(lammps) :: self

    lmp_version = lammps_version(self%handle)
  END FUNCTION lmp_version

  ! ----------------------------------------------------------------------
  ! local helper functions
  ! copy fortran string to zero terminated c string
  ! ----------------------------------------------------------------------
  FUNCTION f2c_string(f_string) RESULT(ptr)
    CHARACTER (LEN=*), INTENT(IN)           :: f_string
    CHARACTER (LEN=1, KIND=c_char), POINTER :: c_string(:)
    TYPE(c_ptr) :: ptr
    INTEGER(c_size_t) :: i, n

    n = LEN_TRIM(f_string)
    ptr = lammps_malloc(n+1)
    CALL C_F_POINTER(ptr, c_string, [1])
    DO i=1, n
        c_string(i) = f_string(i:i)
    END DO
    c_string(n+1) = c_null_char
  END FUNCTION f2c_string
END MODULE LIBLAMMPS

! vim: ts=2 sts=2 sw=2 et

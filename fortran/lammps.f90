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
      c_int, c_char, c_null_char, c_double, c_size_t, c_f_pointer
  USE MPI_F08, only : MPI_Comm

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lammps

  TYPE lammps
      TYPE(c_ptr) :: handle
    CONTAINS
      PROCEDURE :: close              => lmp_close
      PROCEDURE :: file               => lmp_file
      PROCEDURE :: command            => lmp_command
      PROCEDURE :: commands_list      => lmp_commands_list
      PROCEDURE :: commands_string    => lmp_commands_string
      PROCEDURE :: version            => lmp_version
      PROCEDURE :: get_natoms         => lmp_get_natoms
  END TYPE lammps

  INTERFACE lammps
      MODULE PROCEDURE lmp_open, lmp_open_str, lmp_open_f08, lmp_open_f08_str
  END INTERFACE lammps

  ! interface definitions for calling functions in library.cpp
  INTERFACE
      FUNCTION lammps_open(argc,argv,comm) &
          BIND(C, name='lammps_open_fortran')
        IMPORT :: c_ptr, c_int
        ! Technically, comm is of type INTEGER (default kind),
        ! not INTEGER (C_int) but...this avoids the warning
        INTEGER(c_int), VALUE, INTENT(in)     :: argc, comm
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: argv
        TYPE(c_ptr)                           :: lammps_open
      END FUNCTION lammps_open

      FUNCTION lammps_open_no_mpi(argc,argv,handle) &
          BIND(C, name='lammps_open_no_mpi')
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(in)     :: argc
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: argv
        ! 06/24/22: bug fix - INTENT(OUT) -> VALUE (used to send C_NULL_PTR)
        TYPE(c_ptr), VALUE                    :: handle
        TYPE(c_ptr)                           :: lammps_open_no_mpi
      END FUNCTION lammps_open_no_mpi

      SUBROUTINE lammps_close(handle) BIND(C, name='lammps_close')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
      END SUBROUTINE lammps_close

      SUBROUTINE lammps_mpi_init() BIND(C, name='lammps_mpi_init')
      END SUBROUTINE lammps_mpi_init

      SUBROUTINE lammps_mpi_finalize() BIND(C, name='lammps_mpi_finalize')
      END SUBROUTINE lammps_mpi_finalize

      SUBROUTINE lammps_kokkos_finalize() BIND(C, name='lammps_kokkos_finalize')
      END SUBROUTINE lammps_kokkos_finalize

      SUBROUTINE lammps_file(handle,filename) BIND(C, name='lammps_file')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: filename
      END SUBROUTINE lammps_file

      SUBROUTINE lammps_command(handle,cmd) BIND(C, name='lammps_command')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: cmd
      END SUBROUTINE lammps_command

      SUBROUTINE lammps_commands_list(handle,ncmd,cmds) &
          BIND(C, name='lammps_commands_list')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE :: handle
        INTEGER(c_int), VALUE, INTENT(in)     :: ncmd
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: cmds
      END SUBROUTINE lammps_commands_list

      SUBROUTINE lammps_commands_string(handle,str) &
          BIND(C, name='lammps_commands_string')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
        TYPE(c_ptr), VALUE :: str
      END SUBROUTINE lammps_commands_string

      FUNCTION lammps_malloc(size) BIND(C, name='malloc')
        IMPORT :: c_ptr, c_size_t
        INTEGER(c_size_t), value :: size
        TYPE(c_ptr) :: lammps_malloc
      END FUNCTION lammps_malloc

      SUBROUTINE lammps_free(ptr) BIND(C, name='lammps_free')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: ptr
      END SUBROUTINE lammps_free

      FUNCTION lammps_version(handle) BIND(C, name='lammps_version')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE :: handle
        INTEGER(c_int) :: lammps_version
      END FUNCTION lammps_version

      FUNCTION lammps_get_natoms(handle) BIND(C, name='lammps_get_natoms')
        IMPORT :: c_ptr, c_double
        TYPE(c_ptr), VALUE :: handle
        REAL(c_double) :: lammps_get_natoms
      END FUNCTION lammps_get_natoms
  END INTERFACE

CONTAINS

  ! Fortran wrappers and helper functions.

  ! Constructor for the LAMMPS class.
  ! Combined wrapper around lammps_open_fortran() and lammps_open_no_mpi()
  TYPE(lammps) FUNCTION lmp_open(args,comm)
    IMPLICIT NONE
    INTEGER,INTENT(in), OPTIONAL :: comm
    CHARACTER(len=*), INTENT(in), OPTIONAL :: args(:)
    TYPE(c_ptr), ALLOCATABLE     :: argv(:)
    TYPE(c_ptr)                  :: dummy=c_null_ptr
    INTEGER (C_int) :: i, argc

    IF (PRESENT(args)) THEN
        ! convert argument list to c style
        argc = SIZE(args)
        ALLOCATE(argv(argc))
        DO i=1,argc
           argv(i) = f2c_string(args(i))
        END DO
    ELSE
        argc = 1
        ALLOCATE(argv(1))
        argv(1) = f2c_string("liblammps")
    ENDIF

    IF (PRESENT(comm)) THEN
        lmp_open%handle = lammps_open(argc,argv,comm)
    ELSE
        lmp_open%handle = lammps_open_no_mpi(argc,argv,dummy)
    END IF

    ! Clean up allocated memory
    DO i=1,argc
        CALL lammps_free(argv(i))
    END DO
    DEALLOCATE(argv)
  END FUNCTION lmp_open

  ! Alternative constructor for the LAMMPS class for users of MPI_F08
  TYPE (lammps) FUNCTION lmp_open_f08(args,comm)
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: args(:)
    TYPE (MPI_Comm), INTENT(IN) :: comm
    IF ( PRESENT(args) ) THEN
      lmp_open_f08 = lmp_open(args,comm%MPI_VAL)
    ELSE
      lmp_open_f08 = lmp_open(comm=comm%MPI_VAL)
    END IF
  END FUNCTION lmp_open_f08

  ! Alternative constructor for the LAMMPS class.
  ! Overloaded version of lmp_open that accepts a string (e.g., the entire
  ! command line as returned by GET_COMMAND) instead of an array
  TYPE (lammps) FUNCTION lmp_open_str(str,comm)
    CHARACTER(LEN=*), INTENT(IN) :: str
    INTEGER,INTENT(IN), OPTIONAL :: comm
    CHARACTER(LEN=:), ALLOCATABLE :: args(:)
    INTEGER :: i, nargs, beginstr, argnum, maxlen
    CHARACTER(LEN=:), ALLOCATABLE :: tstr

    ALLOCATE ( CHARACTER(LEN=LEN_TRIM(str)) :: tstr )
    tstr = TRIM(str)
    ! find number of arguments (non-consecutive, non-escaped spaces)
    nargs = 1
    i = 2
    beginstr = 1
    IF ( INDEX(tstr, ' ') == 0 ) THEN
      ! Only one argument
      maxlen = LEN(tstr)
    ELSE
      maxlen = INDEX(tstr, ' ')
    END IF
    DO WHILE ( i < LEN(tstr) )
      IF (tstr(i:i) == ' ' .AND. tstr(i-1:i-1) /= '\') THEN
        nargs = nargs + 1
        DO WHILE ( i+1 < LEN(tstr) .AND. tstr(i+1:i+1) == ' ' )
          i = i + 1
        END DO
        IF ( maxlen < (i + 1 - beginstr) ) maxlen = i + 1 - beginstr
        beginstr = i + 1
      END IF
      i = i + 1
    END DO

    ALLOCATE ( CHARACTER(LEN=maxlen) :: args(nargs) )
    beginstr = 1
    i = 1
    argnum = 1
    DO WHILE ( i < LEN(tstr) )
      IF (tstr(i:i) == ' ' .AND. tstr(i-1:i-1) /= '\') THEN
        ! found a real space => split string and store in args(argnum)
        args(argnum) = tstr(beginstr:i-1)
        argnum = argnum + 1
        ! find next non-space
        DO WHILE ( i+1 <= LEN(tstr) .AND. tstr(i+1:i+1) == ' ' )
          i = i + 1
        END DO
        beginstr = i + 1
      END IF
      i = i + 1
    END DO
    ! grab last argument
    args(argnum) = tstr(beginstr:)

    IF ( PRESENT(comm) ) THEN
       lmp_open_str = lmp_open (args, comm)
    ELSE
       lmp_open_str = lmp_open (args)
    END IF

    DEALLOCATE (args, tstr)
  END FUNCTION lmp_open_str

  ! Alternative constructor for users of MPI_F08 that accepts a string
  TYPE (lammps) FUNCTION lmp_open_f08_str(str,comm)
    CHARACTER(LEN=*), INTENT(IN) :: str
    TYPE (MPI_Comm), INTENT(IN) :: comm
    lmp_open_f08_str = lmp_open_str(str,comm%MPI_VAL)
  END FUNCTION lmp_open_f08_str

  ! Combined Fortran wrapper around lammps_close() and lammps_mpi_finalize()
  SUBROUTINE lmp_close(self,finalize)
    IMPLICIT NONE
    CLASS(lammps) :: self
    LOGICAL,INTENT(in),OPTIONAL :: finalize

    CALL lammps_close(self%handle)

    IF (PRESENT(finalize)) THEN
        IF (finalize) THEN
            CALL lammps_kokkos_finalize()
            CALL lammps_mpi_finalize()
        END IF
    END IF
  END SUBROUTINE lmp_close

  INTEGER FUNCTION lmp_version(self)
    IMPLICIT NONE
    CLASS(lammps) :: self

    lmp_version = lammps_version(self%handle)
  END FUNCTION lmp_version

  DOUBLE PRECISION FUNCTION lmp_get_natoms(self)
    IMPLICIT NONE
    CLASS(lammps) :: self

    lmp_get_natoms = lammps_get_natoms(self%handle)
  END FUNCTION lmp_get_natoms

  SUBROUTINE lmp_file(self,filename)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*) :: filename
    TYPE(c_ptr) :: str

    str = f2c_string(filename)
    CALL lammps_file(self%handle,str)
    CALL lammps_free(str)
  END SUBROUTINE lmp_file

  ! equivalent function to lammps_command()
  SUBROUTINE lmp_command(self,cmd)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*) :: cmd
    TYPE(c_ptr) :: str

    str = f2c_string(cmd)
    CALL lammps_command(self%handle,str)
    CALL lammps_free(str)
  END SUBROUTINE lmp_command

  ! equivalent function to lammps_commands_list()
  SUBROUTINE lmp_commands_list(self,cmds)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cmds(:)
    TYPE(c_ptr), ALLOCATABLE     :: cmdv(:)
    INTEGER :: i,ncmd

    ! convert command list to c style
    ncmd = SIZE(cmds)
    ALLOCATE(cmdv(ncmd))
    DO i=1,ncmd
        cmdv(i) = f2c_string(cmds(i))
    END DO

    CALL lammps_commands_list(self%handle,ncmd,cmdv)

    ! Clean up allocated memory
    DO i=1,ncmd
        CALL lammps_free(cmdv(i))
    END DO
    DEALLOCATE(cmdv)
  END SUBROUTINE lmp_commands_list

  ! equivalent function to lammps_commands_string()
  SUBROUTINE lmp_commands_string(self,str)
    IMPLICIT NONE
    CLASS(lammps) :: self
    CHARACTER(len=*) :: str
    TYPE(c_ptr) :: tmp

    tmp = f2c_string(str)
    CALL lammps_commands_string(self%handle,tmp)
    CALL lammps_free(tmp)
  END SUBROUTINE lmp_commands_string

  ! ----------------------------------------------------------------------
  ! local helper functions
  ! ----------------------------------------------------------------------

  ! copy fortran string to zero terminated c string
  FUNCTION f2c_string(f_string) RESULT(ptr)
    CHARACTER (len=*), INTENT(in)           :: f_string
    CHARACTER (len=1, kind=c_char), POINTER :: c_string(:)
    TYPE(c_ptr) :: ptr
    INTEGER(c_size_t) :: i, n

    n = LEN_TRIM(f_string)
    ptr = lammps_malloc(n+1)
    CALL C_F_POINTER(ptr,c_string,[1])
    DO i=1,n
        c_string(i) = f_string(i:i)
    END DO
    c_string(n+1) = c_null_char
  END FUNCTION f2c_string
END MODULE LIBLAMMPS

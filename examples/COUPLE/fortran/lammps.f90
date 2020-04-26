! -------------------------------------------------------------------------
!   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
!   http://lammps.sandia.gov, Sandia National Laboratories
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
!   Karl D. Hammond <karlh@ugcs.caltech.edu>
!   University of Tennessee, Knoxville (USA), 2012
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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_loc, &
      c_int, c_char, c_null_char

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lammps

  TYPE lammps
      TYPE(c_ptr) :: handle
    CONTAINS
      PROCEDURE :: close => lmp_close
      PROCEDURE :: file  => lmp_file
      PROCEDURE :: command => lmp_command
      PROCEDURE :: version => lmp_version
  END TYPE lammps

  INTERFACE lammps
      MODULE PROCEDURE lmp_open
  END INTERFACE lammps

  ! interface definitions for calling functions in library.cpp
  INTERFACE
      SUBROUTINE lammps_open(argc,argv,comm,handle) &
          BIND(C, name='lammps_open_fortran')
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(in)     :: argc, comm
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: argv
        TYPE(c_ptr), INTENT(out)              :: handle
      END SUBROUTINE lammps_open
      SUBROUTINE lammps_open_no_mpi(argc,argv,handle) &
          BIND(C, name='lammps_open_no_mpi')
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(in)     :: argc
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: argv
        TYPE(c_ptr), INTENT(out)              :: handle
      END SUBROUTINE lammps_open_no_mpi
      SUBROUTINE lammps_close(handle) BIND(C, name='lammps_close')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
      END SUBROUTINE lammps_close
      SUBROUTINE lammps_finalize(handle) BIND(C, name='lammps_finalize')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
      END SUBROUTINE lammps_finalize
      FUNCTION lammps_version(handle) BIND(C, name='lammps_version')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE :: handle
        INTEGER(c_int) :: lammps_version
      END FUNCTION lammps_version
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
      SUBROUTINE lammps_free(ptr) BIND(C, name='lammps_free')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: ptr
      END SUBROUTINE lammps_free
  END INTERFACE

CONTAINS

  ! Fortran wrappers and helper functions.

  ! Constructor for the LAMMPS class.
  ! Combined wrapper around lammps_open_fortran() and lammps_open_no_mpi()
  TYPE(lammps) FUNCTION lmp_open(args,comm)
    IMPLICIT NONE
    INTEGER,INTENT(in), OPTIONAL :: comm
    CHARACTER(len=*), INTENT(in) :: args(:)
    TYPE(c_ptr), ALLOCATABLE     :: argv(:)
    INTEGER :: i,argc

    ! convert argument list to c style
    argc = SIZE(args)
    ALLOCATE(argv(argc))
    DO i=1,argc
        argv(i) = f2c_string(args(i))
    END DO
    
    IF (PRESENT(comm)) THEN
        CALL lammps_open(argc,argv,comm,lmp_open%handle)
    ELSE
        CALL lammps_open_no_mpi(argc,argv,lmp_open%handle)
    END IF

    ! Clean up allocated memory
    DO i=1,argc
        CALL lammps_free(argv(i))
    END DO
    DEALLOCATE(argv)
  END FUNCTION lmp_open

  ! Combined Fortran wrapper around lammps_close() and lammps_finalize()
  SUBROUTINE lmp_close(self,finalize)
    IMPLICIT NONE
    CLASS(lammps) :: self
    LOGICAL,INTENT(in),OPTIONAL :: finalize

    CALL lammps_close(self%handle)

    IF (PRESENT(finalize)) THEN
        IF (finalize) THEN
            CALL lammps_finalize(self%handle)
        END IF
    END IF
  END SUBROUTINE lmp_close

  INTEGER FUNCTION lmp_version(self)
    IMPLICIT NONE
    CLASS(lammps) :: self

    lmp_version = lammps_version(self%handle)
  END FUNCTION lmp_version

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

  ! ----------------------------------------------------------------------
  ! local helper functions
  ! copy fortran string to zero terminated c string
  FUNCTION f2c_string(f_string) RESULT(ptr)
    CHARACTER (len=*), INTENT(in)           :: f_string
    CHARACTER (len=1, kind=c_char), POINTER :: c_string(:)
    TYPE(c_ptr) :: ptr
    INTEGER :: i, n

    n = LEN_TRIM(f_string)
    ALLOCATE(c_string(n+1))
    DO i=1,n
        c_string(i) = f_string(i:i)
    END DO
    c_string(n+1) = c_null_char
    ptr = c_loc(c_string(1))
  END FUNCTION f2c_string
END MODULE LIBLAMMPS

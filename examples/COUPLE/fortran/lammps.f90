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

! LAMMPS Fortran library interface implemented as a Fortran 03 style
! module that wraps the C-style library interface in library.cpp
! and library.h using the ISO_C_BINDING module of the fortran compilers.

! Based on the LAMMPS module by Karl D. Hammond <karlh@ugcs.caltech.edu>
! University of Tennessee, Knoxville (USA), 2012

!> LAMMPS module providing the Fortran library interface of LAMMPS.
!!
!! The Fortran module tries to follow the API of the C-library interface
!! as closely as possible. However, since there are some differences in
!! the conventions and specifically how strings are passed, most public
!! functions exported by this module are Fortran code wrappers around
!! the C-library functions that first convert the data as needed to fully
!! comply with the requirements of C/C++.  This is most important
!! for strings, which are not '\0' terminated as required by C/C++.
MODULE LIBLAMMPS

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_loc, c_int, &
      c_char, c_null_char

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lammps

  TYPE lammps
      TYPE(c_ptr) :: handle
    CONTAINS
      PROCEDURE :: close => lmp_close
      PROCEDURE :: file  => lmp_file
      PROCEDURE :: command => lmp_command
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

  ! constructor for the LAMMPS class
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

    ! clean up
    DO i=1,argc
        CALL lammps_free(argv(i))
    END DO
    DEALLOCATE(argv)
  END FUNCTION lmp_open

  ! equivalent function to lammps_close() and lammps_finalize()
  ! optional argument finalize is a logical which, if .true.
  ! causes calling MPI_Finalize after the LAMMPS object is destroyed.
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

  ! equivalent function to lammps_file()
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

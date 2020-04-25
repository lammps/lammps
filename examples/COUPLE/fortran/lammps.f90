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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_int, c_loc

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lammps

  TYPE lammps
      PRIVATE
      TYPE(c_ptr) :: handle
    CONTAINS
      PROCEDURE :: close
  END TYPE lammps

  INTERFACE lammps
      MODULE PROCEDURE lmp_open
  END INTERFACE lammps

  ! interface definitions for calling functions in library.cpp
  INTERFACE
      ! Interface for calling lammps_open_fortran()
      SUBROUTINE lammps_open(argc,argv,comm,handle) &
          BIND(C, name='lammps_open_fortran')
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(in)     :: argc, comm
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: argv
        TYPE(c_ptr), INTENT(out)              :: handle
      END SUBROUTINE lammps_open
      ! Interface for calling lammps_open_no_mpi()
      SUBROUTINE lammps_open_no_mpi(argc,argv,handle) &
          BIND(C, name='lammps_open_no_mpi')
        IMPORT :: c_ptr, c_int
        INTEGER(c_int), VALUE, INTENT(in)     :: argc
        TYPE(c_ptr), DIMENSION(*), INTENT(in) :: argv
        TYPE(c_ptr), INTENT(out)              :: handle
      END SUBROUTINE lammps_open_no_mpi
      ! Interface for calling :cpp:func:`lammps_close`
      SUBROUTINE lammps_close(handle) BIND(C, name='lammps_close')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: handle
      END SUBROUTINE lammps_close
  END INTERFACE

CONTAINS

  ! constructor for the LAMMPS class
  TYPE(lammps) FUNCTION lmp_open(argc,args,comm)
    IMPLICIT NONE
    INTEGER,INTENT(in)           :: argc
    INTEGER,INTENT(in), optional :: comm
    CHARACTER(len=*), INTENT(in) :: args(:)
    TYPE(c_ptr), ALLOCATABLE     :: argv(:)
    INTEGER :: i

    ALLOCATE(argv(argc))
    ! TODO: convert args to argv
    i=0
    IF (PRESENT(comm)) THEN
        CALL lammps_open(i,argv,comm,lmp_open%handle)
    ELSE
        CALL lammps_open_no_mpi(i,argv,lmp_open%handle)
    END IF
  END FUNCTION lmp_open

  ! equivalent function to lammps_close() and lammps_finalize()
  ! optional argument finalize is a logical which, if .true.
  ! causes calling MPI_Finalize after the LAMMPS object is destroyed.
  SUBROUTINE close(self,finalize)
    USE MPI, ONLY: mpi_finalize
    IMPLICIT NONE
    CLASS(lammps) :: self
    LOGICAL,INTENT(in),OPTIONAL :: finalize
    INTEGER :: ierr

    CALL lammps_close(self%handle)

    IF (PRESENT(finalize)) THEN
        IF (finalize) THEN
            CALL mpi_finalize(ierr)
        END IF
    END IF
  END SUBROUTINE close

  ! ----------------------------------------------------------------------
  ! local helper functions
  ! convert fortran string to c string
  PURE FUNCTION fstring2cstring(f_string) RESULT(c_string)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char
    CHARACTER (len=*), INTENT(in) :: f_string
    CHARACTER (len=1, kind=c_char) :: c_string(LEN_TRIM(f_string)+1)
    INTEGER :: i, n

    n = LEN_TRIM(f_string)
    FORALL (i = 1:n)
        c_string(i) = f_string(i:i)
    END FORALL
    c_string(n+1) = c_null_char
  END FUNCTION fstring2cstring
END MODULE LIBLAMMPS

!! -----------------------------------------------------------------------
!   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
!   www.cs.sandia.gov/~sjplimp/lammps.html
!   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
!
!   Copyright (2003) Sandia Corporation.  Under the terms of Contract
!   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
!   certain rights in this software.  This software is distributed under 
!   the GNU General Public License.
!
!   See the README file in the top-level LAMMPS directory.
!--------------------------------------------------------------------------

!! ------------------------------------------------------------------------
!   Contributing author:  Karl D. Hammond <karlh@ugcs.caltech.edu>
!                         University of Tennessee, Knoxville (USA), 2012
!--------------------------------------------------------------------------

!! LAMMPS, a Fortran 2003 module containing an interface between Fortran
!! programs and the C-style functions in library.cpp that ship with LAMMPS.
!! This file should be accompanied by LAMMPS-wrapper.cpp and LAMMPS-wrapper.h,
!! which define wrapper functions that ease portability and enforce array
!! dimensions.
!!
!! Everything in this module should be 100% portable by way of Fortran 2003's
!! ISO_C_BINDING intrinsic module.  See the README for instructions for
!! compilation and use.
!!
!! Here are the PUBLIC functions and subroutines included in this module.
!!    subroutine lammps_open (command_line, communicator, ptr)
!!    subroutine lammps_open_no_mpi (command_line, ptr)
!!    subroutine lammps_close (ptr)
!!    subroutine lammps_file (ptr, str)
!!    subroutine lammps_command (ptr, str)
!!    subroutine lammps_free (ptr)
!!    subroutine lammps_extract_global (global, ptr, name)
!!    subroutine lammps_extract_atom (atom, ptr, name)
!!    subroutine lammps_extract_fix (fix, ptr, id, style, type, i, j)
!!    subroutine lammps_extract_compute (compute, ptr, id, style, type)
!!    subroutine lammps_extract_variable (variable, ptr, name, group)
!!    function lammps_get_natoms (ptr)
!!    subroutine lammps_gather_atoms (ptr, name, count, data)
!!    subroutine lammps_scatter_atoms (ptr, name, data)

#define FLERR __FILE__,__LINE__
! The above line allows for similar error checking as is done with standard
! LAMMPS files.

module LAMMPS

   use, intrinsic :: ISO_C_binding, only : C_double, C_int, C_ptr, C_char, &
      C_NULL_CHAR, C_loc, C_F_pointer, lammps_instance => C_ptr
   implicit none
   private
   public :: lammps_open, lammps_open_no_mpi, lammps_close, lammps_file, &
      lammps_command, lammps_free, lammps_extract_global, &
      lammps_extract_atom, lammps_extract_compute, lammps_extract_fix, &
      lammps_extract_variable, lammps_get_natoms, lammps_gather_atoms, &
      lammps_scatter_atoms
   public :: lammps_instance

   !! Functions supplemental to the prototypes in library.h. {{{1
   !! The function definitions (in C++) are contained in LAMMPS-wrapper.cpp.
   !! I would have written the first in Fortran, but the MPI libraries (which
   !! were written in C) have C-based functions to convert from Fortran MPI
   !! handles to C MPI handles, and there is no Fortran equivalent for those
   !! functions.
   interface
      subroutine lammps_open_wrapper (argc, argv, communicator, ptr) &
      bind (C, name='lammps_open_fortran_wrapper')
         import :: C_int, C_ptr
         integer (C_int), value :: argc
         type (C_ptr), dimension(*) :: argv
         integer, value :: communicator
         type (C_ptr) :: ptr
      end subroutine lammps_open_wrapper
      subroutine lammps_actual_error_all (ptr, file, line, str) &
      bind (C, name='lammps_error_all')
         import :: C_int, C_char, C_ptr
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*), intent(in) :: file, str
         integer (C_int), value :: line
      end subroutine lammps_actual_error_all
      function lammps_get_ntypes (ptr) result (ntypes) &
      bind (C, name='lammps_get_ntypes')
         import :: C_int, C_ptr
         type (C_ptr), value :: ptr
         integer (C_int) :: ntypes
      end function lammps_get_ntypes
      function lammps_actual_extract_compute_vectorsize (ptr, id, style) &
      result (vectorsize) bind (C, name='lammps_extract_compute_vectorsize')
         import :: C_int, C_char, C_ptr
         integer (C_int) :: vectorsize
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: id
         integer (C_int), value :: style
      end function lammps_actual_extract_compute_vectorsize
      subroutine lammps_actual_extract_compute_arraysize (ptr, id, style, &
            nrows, ncols) bind (C, name='lammps_extract_compute_arraysize')
         import :: C_int, C_char, C_ptr
         integer (C_int) :: arraysize
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: id
         integer (C_int), value :: style
         integer (C_int) :: nrows, ncols
      end subroutine lammps_actual_extract_compute_arraysize
      function lammps_actual_extract_fix_vectorsize (ptr, id, style) &
      result (vectorsize) bind (C, name='lammps_extract_fix_vectorsize')
         import :: C_int, C_char, C_ptr
         integer (C_int) :: vectorsize
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: id
         integer (C_int), value :: style
      end function lammps_actual_extract_fix_vectorsize
      subroutine lammps_actual_extract_fix_arraysize (ptr, id, style, &
            nrows, ncols) bind (C, name='lammps_extract_fix_arraysize')
         import :: C_int, C_char, C_ptr
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: id
         integer (C_int), value :: style
         integer (C_int) :: nrows, ncols
      end subroutine lammps_actual_extract_fix_arraysize
   end interface

   !! Functions/subroutines defined in library.h and library.cpp {{{1
   interface
      subroutine lammps_actual_open_no_mpi (argc, argv, ptr) &
      bind (C, name='lammps_open_no_mpi')
         import :: C_int, C_ptr
         integer (C_int), value :: argc
         type (C_ptr), dimension(*) :: argv
         type (C_ptr) :: ptr
      end subroutine lammps_actual_open_no_mpi

      subroutine lammps_close (ptr) bind (C, name='lammps_close')
         import :: C_ptr
         type (C_ptr), value :: ptr
      end subroutine lammps_close

      subroutine lammps_actual_file (ptr, str) bind (C, name='lammps_file')
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: str
      end subroutine lammps_actual_file

      function lammps_actual_command (ptr, str) result (command) &
      bind (C, name='lammps_command')
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: str
         type (C_ptr) :: command
      end function lammps_actual_command

      subroutine lammps_free (ptr) bind (C, name='lammps_free')
         import :: C_ptr
         type (C_ptr), value :: ptr
      end subroutine lammps_free

      function lammps_actual_extract_global (ptr, name) &
      bind (C, name='lammps_extract_global') result (global)
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name
         type (C_ptr) :: global
      end function lammps_actual_extract_global

      function lammps_actual_extract_atom (ptr, name) &
      bind (C, name='lammps_extract_atom') result (atom)
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name
         type (C_ptr) :: atom
      end function lammps_actual_extract_atom

      function lammps_actual_extract_compute (ptr, id, style, type) &
      result (compute) bind (C, name='lammps_extract_compute')
         import :: C_ptr, C_char, C_int
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: id
         integer (C_int), value :: style, type
         type (C_ptr) :: compute
      end function lammps_actual_extract_compute

      function lammps_actual_extract_fix (ptr, id, style, type, i, j) &
      result (fix) bind (C, name='lammps_extract_fix')
         import :: C_ptr, C_char, C_int
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: id
         integer (C_int), value :: style, type, i, j
         type (C_ptr) :: fix
      end function lammps_actual_extract_fix

      function lammps_actual_extract_variable (ptr, name, group) &
      result (variable) bind (C, name='lammps_extract_variable')
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name, group
         type (C_ptr) :: variable
      end function lammps_actual_extract_variable

      function lammps_get_natoms (ptr) result (natoms) &
      bind (C, name='lammps_get_natoms')
         import :: C_ptr, C_int
         type (C_ptr), value :: ptr
         integer (C_int) :: natoms
      end function lammps_get_natoms

      subroutine lammps_actual_gather_atoms (ptr, name, type, count, data) &
      bind (C, name='lammps_gather_atoms')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count
      end subroutine lammps_actual_gather_atoms

      subroutine lammps_actual_scatter_atoms (ptr, name, type, count, data) &
      bind (C, name='lammps_scatter_atoms')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count
      end subroutine lammps_actual_scatter_atoms
   end interface

   ! Generic functions for the wrappers below {{{1

   ! Check the dimensions of the arrays these return; they are not always
   ! easy to find.  Note that I consider returning pointers to arbitrary
   ! memory locations with no information as to array size/shape to be
   ! extremely sloppy and error-prone.  It would appear the Fortran standards
   ! committee would agree, as they chose not to allow that sort of nonsense.

   interface lammps_extract_global
      module procedure lammps_extract_global_i, lammps_extract_global_r, &
         lammps_extract_global_dp
   end interface lammps_extract_global

   interface lammps_extract_atom
      module procedure lammps_extract_atom_ia, lammps_extract_atom_ra, &
         lammps_extract_atom_dpa, lammps_extract_atom_dp2a, &
         lammps_extract_atom_r2a
   end interface lammps_extract_atom

   interface lammps_extract_compute
      module procedure lammps_extract_compute_r, lammps_extract_compute_dp, &
         lammps_extract_compute_ra, lammps_extract_compute_dpa, &
         lammps_extract_compute_r2a, lammps_extract_compute_dp2a
   end interface lammps_extract_compute

   interface lammps_extract_fix
      module procedure lammps_extract_fix_r, lammps_extract_fix_dp, &
         lammps_extract_fix_ra, lammps_extract_fix_dpa, &
         lammps_extract_fix_r2a, lammps_extract_fix_dp2a
   end interface lammps_extract_fix

   interface lammps_extract_variable
      module procedure lammps_extract_variable_i, &
         lammps_extract_variable_dp, &
         lammps_extract_variable_r, &
         lammps_extract_variable_ra, &
         lammps_extract_variable_ia, &
         lammps_extract_variable_dpa
   end interface lammps_extract_variable

   interface lammps_gather_atoms
      module procedure lammps_gather_atoms_ia, lammps_gather_atoms_dpa, &
         lammps_gather_atoms_ra
   end interface lammps_gather_atoms

   interface lammps_scatter_atoms
      module procedure lammps_scatter_atoms_ia, lammps_scatter_atoms_dpa, &
         lammps_scatter_atoms_ra
   end interface lammps_scatter_atoms

contains !! Wrapper functions local to this module {{{1

   subroutine lammps_open (command_line, communicator, ptr)
      character (len=*), intent(in) :: command_line
      integer, intent(in) :: communicator
      type (C_ptr) :: ptr
      integer (C_int) :: argc
      type (C_ptr), dimension(:), allocatable :: argv
      character (kind=C_char), dimension(len_trim(command_line)+1), target :: &
         c_command_line
      c_command_line = string2Cstring (command_line)
      call Cstring2argcargv (c_command_line, argc, argv)
      call lammps_open_wrapper (argc, argv, communicator, ptr)
      deallocate (argv)
   end subroutine lammps_open

!-----------------------------------------------------------------------------

   subroutine lammps_open_no_mpi (command_line, ptr)
      character (len=*), intent(in) :: command_line
      type (C_ptr) :: ptr
      integer (C_int) :: argc
      type (C_ptr), dimension(:), allocatable :: argv
      character (kind=C_char), dimension(len_trim(command_line)+1), target :: &
         c_command_line
      c_command_line = string2Cstring (command_line)
      call Cstring2argcargv (c_command_line, argc, argv)
      call lammps_actual_open_no_mpi (argc, argv, ptr)
      deallocate (argv)
   end subroutine lammps_open_no_mpi

!-----------------------------------------------------------------------------

   subroutine lammps_file (ptr, str)
      type (C_ptr) :: ptr
      character (len=*) :: str
      character (kind=C_char), dimension(len_trim(str)+1) :: Cstr
      Cstr = string2Cstring (str)
      call lammps_actual_file (ptr, Cstr)
   end subroutine lammps_file

!-----------------------------------------------------------------------------

   subroutine lammps_command (ptr, str)
      type (C_ptr) :: ptr
      character (len=*) :: str
      character (kind=C_char), dimension(len_trim(str)+1) :: Cstr
      type (C_ptr) :: dummy
      Cstr = string2Cstring (str)
      dummy = lammps_actual_command (ptr, Cstr)
   end subroutine lammps_command

!-----------------------------------------------------------------------------

! lammps_extract_global {{{2
   function lammps_extract_global_Cptr (ptr, name) result (global)
      type (C_ptr) :: global
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      Cname = string2Cstring (name)
      global = lammps_actual_extract_global (ptr, Cname)
   end function lammps_extract_global_Cptr
   subroutine lammps_extract_global_i (global, ptr, name)
      integer, intent(out) :: global
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      integer (C_int), pointer :: Fptr
      Cptr = lammps_extract_global_Cptr (ptr, name)
      call C_F_pointer (Cptr, Fptr)
      global = Fptr
      nullify (Fptr)
   end subroutine lammps_extract_global_i
   subroutine lammps_extract_global_dp (global, ptr, name)
      double precision, intent(out) :: global
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      real (C_double), pointer :: Fptr
      Cptr = lammps_extract_global_Cptr (ptr, name)
      call C_F_pointer (Cptr, Fptr)
      global = Fptr
      nullify (Fptr)
   end subroutine lammps_extract_global_dp
   subroutine lammps_extract_global_r (global, ptr, name)
      real :: global
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      real (C_double), pointer :: Fptr
      Cptr = lammps_extract_global_Cptr (ptr, name)
      call C_F_pointer (Cptr, Fptr)
      global = real (Fptr)
      nullify (Fptr)
   end subroutine lammps_extract_global_r

!-----------------------------------------------------------------------------

! lammps_extract_atom {{{2
   function lammps_extract_atom_Cptr (ptr, name) result (atom)
      type (C_ptr) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      Cname = string2Cstring (name)
      atom = lammps_actual_extract_atom (ptr, Cname)
   end function lammps_extract_atom_Cptr
   subroutine lammps_extract_atom_ia (atom, ptr, name)
      integer, dimension(:), allocatable, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      integer (C_int), pointer :: Fptr
      integer :: natoms
      natoms = lammps_get_natoms (ptr)
      allocate (atom(natoms))
      Cptr = lammps_extract_atom_Cptr (ptr, name)
      call C_F_pointer (Cptr, Fptr, (/natoms/))
      atom = Fptr
      nullify (Fptr)
   end subroutine lammps_extract_atom_ia
   subroutine lammps_extract_atom_dpa (atom, ptr, name)
      double precision, dimension(:), allocatable, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      real (C_double), dimension(:), pointer :: Fptr
      integer :: nelements
      if ( name == 'mass' ) then
         nelements = lammps_get_ntypes (ptr)
      else if ( name == 'x' .or. name == 'v' .or. name == 'f' ) then
         ! We should not be getting 'x' or 'v' or 'f' here!
         call lammps_error_all (ptr, FLERR, 'You cannot extract those atom&
            & data (x, v, or f) into a rank 1 array.')
         return
      else
         ! Everything else we can get is probably nlocal units long
         call lammps_extract_global_i (nelements, ptr, 'nlocal')
      end if
      allocate (atom(nelements))
      Cptr = lammps_extract_atom_Cptr (ptr, name)
      if ( name == 'mass' ) then
         call C_F_pointer (Cptr, Fptr, (/nelements + 1/))
         atom = Fptr(2:) ! LAMMPS starts numbering at 1 (C does not)
      else
         call C_F_pointer (Cptr, Fptr, (/nelements/))
         atom = Fptr
      end if
      nullify (Fptr)
   end subroutine lammps_extract_atom_dpa
   subroutine lammps_extract_atom_ra (atom, ptr, name)
      real, dimension(:), allocatable, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      double precision, dimension(:), allocatable :: d_atom
      call lammps_extract_atom_dpa (d_atom, ptr, name)
      allocate (atom(size(d_atom)))
      atom = real(d_atom)
      deallocate (d_atom)
   end subroutine lammps_extract_atom_ra
   subroutine lammps_extract_atom_dp2a (atom, ptr, name)
      double precision, dimension(:,:), allocatable, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      integer :: nelements
      if ( name /= 'x' .and. name /= 'v' .and. name /= 'f' ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract ' // name // &
            ' into a rank 2 array.')
         return
      end if
      Cptr = lammps_extract_atom_Cptr (ptr, name)
      nelements = lammps_get_natoms (ptr)
      allocate (atom(nelements,3))
      atom = Cdoublestar_to_2darray (Cptr, nelements, 3)
   end subroutine lammps_extract_atom_dp2a
   subroutine lammps_extract_atom_r2a (atom, ptr, name)
      real, dimension(:,:), allocatable, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      double precision, dimension(:,:), allocatable :: d_atom
      call lammps_extract_atom_dp2a (d_atom, ptr, name)
      if ( allocated (d_atom) ) then
         allocate (atom(size(d_atom,1), size(d_atom,2)))
      else
         return
      end if
      atom = real(d_atom)
      deallocate (d_atom)
   end subroutine lammps_extract_atom_r2a

!-----------------------------------------------------------------------------

! lammps_extract_compute {{{2
   function lammps_extract_compute_Cptr (ptr, id, style, type) result (compute)
      type (C_ptr) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      integer (kind=C_int) :: Cstyle, Ctype
      character (kind=C_char), dimension(len_trim(id)+1) :: Cid
      Cid = string2Cstring (id)
      Cstyle = style
      Ctype = type
      compute = lammps_actual_extract_compute (ptr, Cid, Cstyle, Ctype)
   end function lammps_extract_compute_Cptr
   subroutine lammps_extract_compute_dp (compute, ptr, id, style, type)
      double precision, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      type (C_ptr) :: Cptr
      real (C_double), pointer :: Fptr
      ! The only valid values of (style,type) are (0,0) for scalar 'compute'
      if ( style /= 0 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot pack per-atom/local&
            & data into a scalar.')
         return
      end if
      if ( type == 1 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a compute&
            & vector (rank 1) into a scalar.')
         return
      else if ( type == 2 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a compute&
            & array (rank 2) into a scalar.')
         return
      end if
      Cptr = lammps_extract_compute_Cptr (ptr, id, style, type)
      call C_F_pointer (Cptr, Fptr)
      compute = Fptr
      nullify (Fptr)
      ! C pointer should not be freed!
   end subroutine lammps_extract_compute_dp
   subroutine lammps_extract_compute_r (compute, ptr, id, style, type)
      real, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      double precision :: d_compute
      call lammps_extract_compute_dp (d_compute, ptr, id, style, type)
      compute = real(d_compute)
   end subroutine lammps_extract_compute_r
   subroutine lammps_extract_compute_dpa (compute, ptr, id, style, type)
      double precision, dimension(:), allocatable, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      type (C_ptr) :: Cptr
      real (C_double), dimension(:), pointer :: Fptr
      integer :: nelements
      ! Check for the correct dimensionality
      if ( type == 0 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a compute&
            & scalar (rank 0) into a rank 1 variable.')
         return
      else if ( type == 2 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a compute&
            & array (rank 2) into a rank 1 variable.')
         return
      end if
      nelements = lammps_extract_compute_vectorsize (ptr, id, style)
      allocate (compute(nelements))
      Cptr = lammps_extract_compute_Cptr (ptr, id, style, type)
      call C_F_pointer (Cptr, Fptr, (/nelements/))
      compute = Fptr
      nullify (Fptr)
      ! C pointer should not be freed
   end subroutine lammps_extract_compute_dpa
   subroutine lammps_extract_compute_ra (compute, ptr, id, style, type)
      real, dimension(:), allocatable, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      double precision, dimension(:), allocatable :: d_compute
      call lammps_extract_compute_dpa (d_compute, ptr, id, style, type)
      allocate (compute(size(d_compute)))
      compute = real(d_compute)
      deallocate (d_compute)
   end subroutine lammps_extract_compute_ra
   subroutine lammps_extract_compute_dp2a (compute, ptr, id, style, type)
      double precision, dimension(:,:), allocatable, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      type (C_ptr) :: Cptr
      real (C_double), dimension(:,:), pointer :: Fptr
      integer :: nr, nc
      ! Check for the correct dimensionality
      if ( type == 0 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a compute&
            & scalar (rank 0) into a rank 2 variable.')
         return
      else if ( type == 1 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a compute&
            & array (rank 1) into a rank 2 variable.')
         return
      end if
      call lammps_extract_compute_arraysize (ptr, id, style, nr, nc)
      allocate (compute(nr, nc))
      Cptr = lammps_extract_compute_Cptr (ptr, id, style, type)
      call C_F_pointer (Cptr, Fptr, (/nr, nc/))
      compute = Fptr
      nullify (Fptr)
      ! C pointer should not be freed
   end subroutine lammps_extract_compute_dp2a
   subroutine lammps_extract_compute_r2a (compute, ptr, id, style, type)
      real, dimension(:,:), allocatable, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      double precision, dimension(:,:), allocatable :: d_compute
      call lammps_extract_compute_dp2a (d_compute, ptr, id, style, type)
      allocate (compute(size(d_compute,1), size(d_compute,2)))
      compute = real(d_compute)
      deallocate (d_compute)
   end subroutine lammps_extract_compute_r2a

!-----------------------------------------------------------------------------

! lammps_extract_fix {{{2
   function lammps_extract_fix_Cptr (ptr, id, style, type, i, j) &
   result (fix)
      type (C_ptr) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      character (kind=C_char), dimension(len_trim(id)+1) :: Cid
      integer (kind=C_int) :: Cstyle, Ctype, Ci, Cj
      Cid = string2Cstring (id)
      Cstyle = style
      Ctype = type
      Ci = i - 1  ! This is for consistency with the values from f_ID[i],
      Cj = j - 1  ! which is different from what library.cpp uses!
      if ( (type >= 1 .and. Ci < 0) .or. &
           (type == 2 .and. (Ci < 0 .or. Cj < 0) ) ) then
         call lammps_error_all (ptr, FLERR, 'Index out of range in&
            & lammps_extract_fix')
      end if
      fix = lammps_actual_extract_fix (ptr, Cid, Cstyle, Ctype, Ci, Cj)
   end function lammps_extract_fix_Cptr
   subroutine lammps_extract_fix_dp (fix, ptr, id, style, type, i, j)
      double precision, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      type (C_ptr) :: Cptr
      real (C_double), pointer :: Fptr
      ! Check for the correct dimensionality
      if ( style /= 0 ) then
         select case (type)
         case (0)
            call lammps_error_all (ptr, FLERR, 'There is no per-atom or local&
               & scalar data available from fixes.')
         case (1)
            call lammps_error_all (ptr, FLERR, 'You cannot extract a fix''s &
               &per-atom/local vector (rank 1) into a scalar.')
         case (2)
            call lammps_error_all (ptr, FLERR, 'You cannot extract a fix''s &
               &per-atom/local array (rank 2) into a scalar.')
         case default
            call lammps_error_all (ptr, FLERR, 'Invalid extract_fix style&
               & value.')
         end select
         return
      end if
      Cptr = lammps_extract_fix_Cptr (ptr, id, style, type, i, j)
      call C_F_pointer (Cptr, Fptr)
      fix = Fptr
      nullify (Fptr)
      ! Memory is only allocated for "global" fix variables
      if ( style == 0 ) call lammps_free (Cptr)
   end subroutine lammps_extract_fix_dp
   subroutine lammps_extract_fix_r (fix, ptr, id, style, type, i, j)
      real, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      double precision :: d_fix
      call lammps_extract_fix_dp (d_fix, ptr, id, style, type, i, j)
      fix = real(d_fix)
   end subroutine lammps_extract_fix_r
   subroutine lammps_extract_fix_dpa (fix, ptr, id, style, type, i, j)
      double precision, dimension(:), allocatable, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      type (C_ptr) :: Cptr
      real (C_double), dimension(:), pointer :: Fptr
      integer :: fix_len
      ! Check for the correct dimensionality
      if ( style == 0 ) then
         call lammps_error_all (ptr, FLERR, 'You can''t extract the&
            & whole vector from global fix data')
         return
      else if ( type == 0 ) then
         call lammps_error_all (ptr, FLERR, 'You can''t extract a fix&
            & scalar into a rank 1 variable')
         return
      else if ( type == 2 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a fix&
            & array into a rank 1 variable.')
         return
      else if ( type /= 1 ) then
         call lammps_error_all (ptr, FLERR, 'Invalid type for fix extraction.')
         return
      end if
      fix_len = lammps_extract_fix_vectorsize (ptr, id, style)
      allocate (fix(fix_len))
      Cptr = lammps_extract_fix_Cptr (ptr, id, style, type, i, j)
      call C_F_pointer (Cptr, Fptr, (/fix_len/))
      fix = Fptr
      nullify (Fptr)
      ! Memory is only allocated for "global" fix variables
      if ( style == 0 ) call lammps_free (Cptr)
   end subroutine lammps_extract_fix_dpa
   subroutine lammps_extract_fix_ra (fix, ptr, id, style, type, i, j)
      real, dimension(:), allocatable, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      double precision, dimension(:), allocatable :: d_fix
      call lammps_extract_fix_dpa (d_fix, ptr, id, style, type, i, j)
      allocate (fix(size(d_fix)))
      fix = real(d_fix)
      deallocate (d_fix)
   end subroutine lammps_extract_fix_ra
   subroutine lammps_extract_fix_dp2a (fix, ptr, id, style, type, i, j)
      double precision, dimension(:,:), allocatable, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      type (C_ptr) :: Cptr
      real (C_double), dimension(:,:), pointer :: Fptr
      integer :: nr, nc
      ! Check for the correct dimensionality
      if ( style == 0 ) then
         call lammps_error_all (ptr, FLERR, 'It is not possible to extract the&
            & entire array from global fix data.')
         return
      else if ( type == 0 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a fix&
            & scalar (rank 0) into a rank 2 variable.')
         return
      else if ( type == 1 ) then
         call lammps_error_all (ptr, FLERR, 'You cannot extract a fix&
            & vector (rank 1) into a rank 2 variable.')
         return
      end if
      call lammps_extract_fix_arraysize (ptr, id, style, nr, nc)
      allocate (fix(nr, nc))
      Cptr = lammps_extract_fix_Cptr (ptr, id, style, type, i, j)
      call C_F_pointer (Cptr, Fptr, (/nr, nc/))
      fix = Fptr
      nullify (Fptr)
      ! C pointer should not be freed
   end subroutine lammps_extract_fix_dp2a
   subroutine lammps_extract_fix_r2a (fix, ptr, id, style, type, i, j)
      real, dimension(:,:), allocatable, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      double precision, dimension(:,:), allocatable :: d_fix
      call lammps_extract_fix_dp2a (d_fix, ptr, id, style, type, i, j)
      allocate (fix(size(d_fix,1), size(d_fix,2)))
      fix = real(d_fix)
      deallocate (d_fix)
   end subroutine lammps_extract_fix_r2a

!-----------------------------------------------------------------------------

! lammps_extract_variable {{{2
   function lammps_extract_variable_Cptr (ptr, name, group) result (variable)
      type (C_ptr) :: ptr, variable
      character (len=*) :: name
      character (len=*), optional :: group
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      character (kind=C_char), dimension(:), allocatable :: Cgroup
      Cname = string2Cstring (name)
      if ( present(group) ) then
         allocate (Cgroup(len_trim(group)+1))
         Cgroup = string2Cstring (group)
      else
         allocate (Cgroup(1))
         Cgroup(1) = C_NULL_CHAR
      end if
      variable = lammps_actual_extract_variable (ptr, Cname, Cgroup)
      deallocate (Cgroup)
   end function lammps_extract_variable_Cptr
   subroutine lammps_extract_variable_i (variable, ptr, name, group)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
      integer, intent(out) :: variable
      type (C_ptr) :: Cptr
      integer (C_int), pointer :: Fptr
      if ( present(group) ) then
         Cptr = lammps_extract_variable_Cptr (ptr, name, group)
      else
         Cptr = lammps_extract_variable_Cptr (ptr, name)
      end if
      call C_F_pointer (Cptr, Fptr)
      variable = Fptr
      nullify (Fptr)
      call lammps_free (Cptr)
   end subroutine lammps_extract_variable_i
   subroutine lammps_extract_variable_dp (variable, ptr, name, group)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
      double precision, intent(out) :: variable
      type (C_ptr) :: Cptr
      real (C_double), pointer :: Fptr
      if ( present(group) ) then
         Cptr = lammps_extract_variable_Cptr (ptr, name, group)
      else
         Cptr = lammps_extract_variable_Cptr (ptr, name)
      end if
      call C_F_pointer (Cptr, Fptr)
      variable = Fptr
      nullify (Fptr)
      call lammps_free (Cptr)
   end subroutine lammps_extract_variable_dp
   subroutine lammps_extract_variable_r (variable, ptr, name, group)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
      real, intent(out) :: variable
      double precision :: d_var
      if ( present (group) ) then
         call lammps_extract_variable_dp (d_var, ptr, name, group)
      else
         call lammps_extract_variable_dp (d_var, ptr, name)
      end if
      variable = real(d_var)
   end subroutine lammps_extract_variable_r

   subroutine lammps_extract_variable_ia (variable, ptr, name, group)
      integer, dimension(:), allocatable, intent(out) :: variable
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
      type (C_ptr) :: Cptr
      integer (C_int), dimension(:), pointer :: Fptr
      integer :: natoms
      nullify (Fptr)
      if ( present(group) ) then
         Cptr = lammps_extract_variable_Cptr (ptr, name, group)
      else
         Cptr = lammps_extract_variable_Cptr (ptr, name)
      end if
      natoms = lammps_get_natoms (ptr)
      allocate (variable(natoms))
      call C_F_pointer (Cptr, Fptr, (/natoms/))
      variable = Fptr
      nullify (Fptr)
      call lammps_free (Cptr)
   end subroutine lammps_extract_variable_ia
   subroutine lammps_extract_variable_dpa (variable, ptr, name, group)
      double precision, dimension(:), allocatable, intent(out) :: variable
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
      type (C_ptr) :: Cptr
      real (C_double), dimension(:), pointer :: Fptr
      integer :: natoms
      if ( present(group) ) then
         Cptr = lammps_extract_variable_Cptr (ptr, name, group)
      else
         Cptr = lammps_extract_variable_Cptr (ptr, name)
      end if
      natoms = lammps_get_natoms (ptr)
      allocate (variable(natoms))
      call C_F_pointer (Cptr, Fptr, (/natoms/))
      variable = Fptr
      nullify (Fptr)
      call lammps_free (Cptr)
   end subroutine lammps_extract_variable_dpa
   subroutine lammps_extract_variable_ra (variable, ptr, name, group)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
      real, dimension(:), allocatable, intent(out) :: variable
      double precision, dimension(:), allocatable :: d_var
      if ( present (group) ) then
         call lammps_extract_variable_dpa (d_var, ptr, name, group)
      else
         call lammps_extract_variable_dpa (d_var, ptr, name)
      end if
      allocate (variable(size(d_var)))
      variable = real(d_var)
      deallocate (d_var)
   end subroutine lammps_extract_variable_ra

!-------------------------------------------------------------------------2}}}

   subroutine lammps_gather_atoms_ia (ptr, name, count, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count
      integer, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      integer (C_int), dimension(:), pointer :: Fdata
      integer (C_int) :: natoms
      character (kind=C_char), dimension(len_trim(name)) :: Cname
      integer (C_int), parameter :: Ctype = 0
      integer (C_int) :: Ccount
      natoms = lammps_get_natoms (ptr)
      Cname = string2Cstring (name)
      if ( count /= 1 .and. count /= 3 ) then
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      else
         Ccount = count
      end if
      allocate ( Fdata(count*natoms) )
      allocate ( data(count*natoms) )
      Cdata = C_loc (Fdata(1))
      call lammps_actual_gather_atoms (ptr, Cname, Ctype, Ccount, Cdata)
      data = Fdata
      deallocate (Fdata)
   end subroutine lammps_gather_atoms_ia
   subroutine lammps_gather_atoms_dpa (ptr, name, count, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count
      double precision, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      real (C_double), dimension(:), pointer :: Fdata
      integer (C_int) :: natoms
      character (kind=C_char), dimension(len_trim(name)) :: Cname
      integer (C_int), parameter :: Ctype = 1
      integer (C_int) :: Ccount
      natoms = lammps_get_natoms (ptr)
      Cname = string2Cstring (name)
      if ( count /= 1 .and. count /= 3 ) then
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      else
         Ccount = count
      end if
      allocate ( Fdata(count*natoms) )
      allocate ( data(count*natoms) )
      Cdata = C_loc (Fdata(1))
      call lammps_actual_gather_atoms (ptr, Cname, Ctype, Ccount, Cdata)
      data = Fdata(:)
      deallocate (Fdata)
   end subroutine lammps_gather_atoms_dpa
   subroutine lammps_gather_atoms_ra (ptr, name, count, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count
      real, dimension(:), allocatable, intent(out) :: data
      double precision, dimension(:), allocatable :: d_data
      call lammps_gather_atoms_dpa (ptr, name, count, d_data)
      allocate (data(size(d_data)))
      data = d_data
      deallocate (d_data)
   end subroutine lammps_gather_atoms_ra

!-----------------------------------------------------------------------------

   subroutine lammps_scatter_atoms_ia (ptr, name, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: data
      integer (kind=C_int) :: natoms, Ccount
      integer (kind=C_int), parameter :: Ctype = 0
      character (kind=C_char), dimension(len_trim(name)) :: Cname
      integer, dimension(size(data)), target :: Fdata
      type (C_ptr) :: Cdata
      natoms = lammps_get_natoms (ptr)
      Cname = string2Cstring (name)
      Ccount = size(data) / natoms
      if ( Ccount /= 1 .and. Ccount /= 3 ) &
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      Fdata = data
      Cdata = C_loc (Fdata(1))
      call lammps_actual_scatter_atoms (ptr, Cname, Ctype, Ccount, Cdata)
   end subroutine lammps_scatter_atoms_ia
   subroutine lammps_scatter_atoms_dpa (ptr, name, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      double precision, dimension(:), intent(in) :: data
      integer (kind=C_int) :: natoms, Ccount
      integer (kind=C_int), parameter :: Ctype = 0
      character (kind=C_char), dimension(len_trim(name)) :: Cname
      double precision, dimension(size(data)), target :: Fdata
      type (C_ptr) :: Cdata
      natoms = lammps_get_natoms (ptr)
      Cname = string2Cstring (name)
      Ccount = size(data) / natoms
      if ( Ccount /= 1 .and. Ccount /= 3 ) &
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      Fdata = data
      Cdata = C_loc (Fdata(1))
      call lammps_actual_scatter_atoms (ptr, Cname, Ctype, Ccount, Cdata)
   end subroutine lammps_scatter_atoms_dpa
   subroutine lammps_scatter_atoms_ra (ptr, name, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      real, dimension(:), intent(out) :: data
      double precision, dimension(size(data)) :: d_data
      d_data = real (data, kind(d_data))
      call lammps_scatter_atoms_dpa (ptr, name, d_data)
   end subroutine lammps_scatter_atoms_ra

!-----------------------------------------------------------------------------

!   subroutine lammps_get_coords (ptr, coords)
!      type (C_ptr) :: ptr
!      double precision, dimension(:), allocatable, intent(out) :: coords
!      real (C_double), dimension(:), allocatable, target :: C_coords
!      integer :: natoms
!      natoms = lammps_get_natoms (ptr)
!      allocate (coords(3*natoms))
!      allocate (C_coords(3*natoms))
!      call lammps_actual_get_coords (ptr, C_loc(C_coords))
!      coords = C_coords
!      deallocate (C_coords)
!   end subroutine lammps_get_coords
!
!!-----------------------------------------------------------------------------
!
!   subroutine lammps_put_coords (ptr, coords)
!      type (C_ptr) :: ptr
!      double precision, dimension(:) :: coords
!      real (C_double), dimension(size(coords)) :: C_coords
!      C_coords = coords
!      call lammps_actual_put_coords (ptr, C_coords)
!   end subroutine lammps_put_coords
!
!!-----------------------------------------------------------------------------

   function lammps_extract_compute_vectorsize (ptr, id, style) &
   result (vectorsize)
      integer :: vectorsize
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style
      integer (C_int) :: Cvectorsize, Cstyle
      character (kind=C_char), dimension(len_trim(id)+1) :: Cid
      Cid = string2Cstring (id)
      Cstyle = int(style, C_int)
      Cvectorsize = lammps_actual_extract_compute_vectorsize (ptr, Cid, Cstyle)
      vectorsize = int(Cvectorsize, kind(vectorsize))
   end function lammps_extract_compute_vectorsize

!-----------------------------------------------------------------------------

   function lammps_extract_fix_vectorsize (ptr, id, style) &
   result (vectorsize)
      integer :: vectorsize
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style
      integer (C_int) :: Cvectorsize, Cstyle
      character (kind=C_char), dimension(len_trim(id)+1) :: Cid
      Cid = string2Cstring (id)
      Cstyle = int(style, C_int)
      Cvectorsize = lammps_actual_extract_fix_vectorsize (ptr, Cid, Cstyle)
      vectorsize = int(Cvectorsize, kind(vectorsize))
   end function lammps_extract_fix_vectorsize

!-----------------------------------------------------------------------------

   subroutine lammps_extract_compute_arraysize (ptr, id, style, nrows, ncols)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style
      integer, intent(out) :: nrows, ncols
      integer (C_int) :: Cstyle, Cnrows, Cncols
      character (kind=C_char), dimension(len_trim(id)+1) :: Cid
      Cid = string2Cstring (id)
      Cstyle = int (style, C_int)
      call lammps_actual_extract_compute_arraysize (ptr, Cid, Cstyle, &
         Cnrows, Cncols)
      nrows = int (Cnrows, kind(nrows))
      ncols = int (Cncols, kind(ncols))
   end subroutine lammps_extract_compute_arraysize

!-----------------------------------------------------------------------------

   subroutine lammps_extract_fix_arraysize (ptr, id, style, nrows, ncols)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style
      integer, intent(out) :: nrows, ncols
      integer (C_int) :: Cstyle, Cnrows, Cncols
      character (kind=C_char), dimension(len_trim(id)+1) :: Cid
      Cid = string2Cstring (id)
      Cstyle = int (style, kind(Cstyle))
      call lammps_actual_extract_fix_arraysize (ptr, Cid, Cstyle, &
         Cnrows, Cncols)
      nrows = int (Cnrows, kind(nrows))
      ncols = int (Cncols, kind(ncols))
   end subroutine lammps_extract_fix_arraysize

!-----------------------------------------------------------------------------

   subroutine lammps_error_all (ptr, file, line, str)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: file, str
      integer, intent(in) :: line
      character (kind=C_char), dimension(len_trim(file)+1) :: Cfile
      character (kind=C_char), dimension(len_trim(str)+1) :: Cstr
      integer (C_int) :: Cline
      Cline = int(line, kind(Cline))
      Cfile = string2Cstring (file)
      Cstr = string2Cstring (str)
      call lammps_actual_error_all (ptr, Cfile, Cline, Cstr)
   end subroutine lammps_error_all

!-----------------------------------------------------------------------------

! Locally defined helper functions {{{1

   pure function string2Cstring (string) result (C_string)
      use, intrinsic :: ISO_C_binding, only : C_char, C_NULL_CHAR
      character (len=*), intent(in) :: string
      character (len=1, kind=C_char) :: C_string (len_trim(string)+1)
      integer :: i, n
      n = len_trim (string)
      forall (i = 1:n)
         C_string(i) = string(i:i)
      end forall
      C_string(n+1) = C_NULL_CHAR
   end function string2Cstring

!-----------------------------------------------------------------------------

   subroutine Cstring2argcargv (Cstring, argc, argv)
   !! Converts a C-style string to argc and argv, that is, words in Cstring
   !! become C-style strings in argv.  IMPORTANT:  Cstring is modified by
   !! this routine!  I would make Cstring local TO this routine and accept
   !! a Fortran-style string instead, but we run into scoping and
   !! allocation problems that way.  This routine assumes the string is
   !! null-terminated, as all C-style strings must be.

      character (kind=C_char), dimension(*), target, intent(inout) :: Cstring
      integer (C_int), intent(out) :: argc
      type (C_ptr), dimension(:), allocatable, intent(out) :: argv

      integer :: StringStart, SpaceIndex, strlen, argnum

      argc = 1_C_int

      ! Find the length of the string
      strlen = 1
      do while ( Cstring(strlen) /= C_NULL_CHAR )
         strlen = strlen + 1
      end do

      ! Find the number of non-escaped spaces
      SpaceIndex = 2
      do while ( SpaceIndex < strlen )
         if ( Cstring(SpaceIndex) == ' ' .and. &
              Cstring(SpaceIndex-1) /= '\' ) then
            argc = argc + 1_C_int
            ! Find the next non-space character
            do while ( Cstring(SpaceIndex+1) == ' ')
               SpaceIndex = SpaceIndex + 1
            end do
         end if
         SpaceIndex = SpaceIndex + 1
      end do

      ! Now allocate memory for argv
      allocate (argv(argc))

      ! Now find the string starting and ending locations
      StringStart = 1
      SpaceIndex = 2
      argnum = 1
      do while ( SpaceIndex < strlen )
         if ( Cstring(SpaceIndex) == ' ' .and. &
              Cstring(SpaceIndex-1) /= '\' ) then
            ! Found a real space => split strings and store this one
            Cstring(Spaceindex) = C_NULL_CHAR ! Replaces space with NULL
            argv(argnum) = C_loc(Cstring(StringStart))
            argnum = argnum + 1
            ! Find the next non-space character
            do while ( Cstring(SpaceIndex+1) == ' ')
               SpaceIndex = SpaceIndex + 1
            end do
            StringStart = SpaceIndex + 1
         else if ( Cstring(SpaceIndex) == ' ' .and. &
              Cstring(SpaceIndex-1) == '\' ) then
            ! Escaped space => remove backslash and move rest of array
            Cstring(SpaceIndex-1:strlen-1) = Cstring(SpaceIndex:strlen)
            strlen = strlen - 1 ! Last character is still C_NULL_CHAR
         end if
         SpaceIndex = SpaceIndex + 1
      end do
      ! Now handle the last argument
      argv(argnum) = C_loc(Cstring(StringStart))

   end subroutine Cstring2argcargv

!-----------------------------------------------------------------------------

   function Cdoublestar_to_2darray (Carray, nrows, ncolumns) result (Farray)

   ! Take a C/C++ array of pointers to pointers to doubles (sort of like a
   ! two-dimensional array, and handled the same way from the programmer's
   ! perspective) into a Fortran-style array.  Note that columns in C still
   ! correspond to columns in Fortran here and the same for rows.

      type (C_ptr), intent(in) :: Carray
      integer, intent(in) :: nrows, ncolumns
      double precision, dimension(nrows, ncolumns) :: Farray
      type (C_ptr), dimension(:), pointer :: C_rows
      real (C_double), dimension(:), pointer :: F_row
      integer :: i

      ! Convert each "C row pointer" into an array of rows
      call C_F_pointer (Carray, C_rows, (/nrows/))
      do i = 1, nrows
         ! Convert each C pointer (an entire row) into a Fortran pointer
         call C_F_pointer (C_rows(i), F_row, (/ncolumns/))
         Farray (i,:) = real(F_row, kind(0.0D0))
      end do

   end function Cdoublestar_to_2darray
! 1}}}

end module LAMMPS

! vim: foldmethod=marker ts=3 sts=3 expandtab

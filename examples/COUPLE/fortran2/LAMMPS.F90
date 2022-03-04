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
!   Contributing author:  Karl D. Hammond <hammondkd@missouri.edu>
!                         University of Tennessee, Knoxville (USA), 2012
!   Updated October 2020 by the author (now at the University of Missouri).
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
!!    integer (kind=C_int) lammps_version (ptr)
!!    subroutine lammps_file (ptr, str)
!!    subroutine lammps_command (ptr, str)
!!    subroutine lammps_commands_list (ptr, cmds)
!!    subroutine lammps_commands_string (ptr, str)
!!    subroutine lammps_free (ptr)
!!    integer function lammps_extract_setting (ptr, name)
!!    subroutine lammps_extract_global (global, ptr, name)
!!    subroutine lammps_extract_box (ptr, boxlo, boxhi, xy, yz, xz, &
!!       periodicity, box_change)
!!    subroutine lammps_extract_atom (atom, ptr, name)
!!    subroutine lammps_extract_compute (compute, ptr, id, style, type)
!!    subroutine lammps_extract_fix (fix, ptr, id, style, type, i, j)
!!    subroutine lammps_extract_variable (variable, ptr, name, group)
!!    double precision function lammps_get_thermo (ptr, name)
!!    double precision function lammps_get_natoms (ptr)
!!    subroutine lammps_set_variable (ptr, name, str, [err])
!!    subroutine lammps_reset_box (ptr, boxlo, boxhi, xy, yz, xz)
!!    subroutine lammps_gather_atoms (ptr, name, count, data)
!!    subroutine lammps_gather_atoms_concat (ptr, name, count, data)
!!    subroutine lammps_gather_atoms_subset (ptr, name, count, ids, data)
!!    subroutine lammps_scatter_atoms (ptr, name, data)
!!    subroutine lammps_scatter_atoms_subset (ptr, name, ids, data)
!!    logical function lammps_config_has_package (package_name)
!!    integer function lammps_config_package_count ()
!!    logical function lammps_config_package_name (index, buffer)
!!    logical function lammps_config_has_gzip_support ()
!!    logical function lammps_config_has_png_support ()
!!    logical function lammps_config_has_jpeg_support ()
!!    logical function lammps_config_has_ffmpeg_support ()
!!    logical function lammps_config_has_exceptions ()
!!    integer function lammps_find_pair_neighlist (ptr, style, exact, nsub,
!!       request)
!!    integer function lammps_find_fix_neighlist (ptr, id, request)
!!    integer function lammps_find_compute_neighlist (ptr, id, request)
!!    integer function lammps_neighlist_num_elements (ptr, idx)
!!    subroutine lammps_neighlist_element_neighbors (ptr, idx, element, iatom,
!!       numneigh, neighbors)
!!    subroutine lammps_create_atoms (ptr, n, id, type, x, v, image,
!1       shrinkexceed)
!!
!! The following are also available if compiled with -DLAMMPS_EXCEPTIONS
!!    function lammps_has_error (ptr)
!!    function lammps_get_last_error_message (ptr, buffer)
!!
!! Note that the following function is not implemented from this interface:
!!    lammps_set_fix_external_callback

#define FLERR __FILE__,__LINE__
! The above line allows for similar error checking as is done with standard
! LAMMPS files.

! This should (?) allow this module to work with settings in lmptype.h
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#ifdef LAMMPS_SMALLBIG
#define C_smallint C_int
#define C_imageint C_int
#define C_tagint C_int
#define C_bigint C_int64_t
#endif

#ifdef LAMMPS_BIGBIG
#define C_smallint C_int
#define C_imageint C_int64_t
#define C_tagint C_int64_t
#define C_bigint C_int64_t
#endif

#ifdef LAMMPS_SMALLSMALL
#define C_smallint C_int
#define C_imageint C_int
#define C_tagint C_int
#define C_bigint C_int
#endif

module LAMMPS

   use, intrinsic :: ISO_C_binding, only : C_double, C_int, C_ptr, C_char, &
      C_NULL_CHAR, C_NULL_PTR, C_loc, C_F_pointer, C_int64_t, &
      lammps_instance => C_ptr
   implicit none
   private

   ! We inherit some ISO_C_BINDING entities for ease of use
   public :: lammps_instance, C_ptr, C_double, C_int
   ! Only the following functions may be called by the user:
   public :: lammps_open, lammps_open_no_mpi, lammps_close, &
      lammps_version, lammps_file, lammps_command, lammps_commands_list, &
      lammps_commands_string, lammps_free, lammps_extract_setting, &
      lammps_extract_global, lammps_extract_box, lammps_extract_atom, &
      lammps_extract_compute, lammps_extract_fix, lammps_extract_variable, &
      lammps_get_thermo, lammps_get_natoms, lammps_set_variable, &
      lammps_reset_box, lammps_gather_atoms, lammps_gather_atoms_concat, &
      lammps_gather_atoms_subset, lammps_scatter_atoms, &
      lammps_scatter_atoms_subset, lammps_config_has_package, &
      lammps_config_package_count, lammps_config_package_name, &
      lammps_config_has_gzip_support, lammps_config_has_png_support, &
      lammps_config_has_jpeg_support, lammps_config_has_ffmpeg_support, &
      lammps_config_has_exceptions, lammps_find_pair_neighlist, &
      lammps_find_fix_neighlist, lammps_find_compute_neighlist, &
      lammps_neighlist_num_elements, lammps_neighlist_element_neighbors, &
      lammps_create_atoms

#ifdef LAMMPS_EXCEPTIONS
   public :: lammps_has_error, lammps_get_last_error_message
#endif

   !! constants for extracting data from computes and fixes
   !! and data types

   INTEGER, PARAMETER :: LMP_STYLE_GLOBAL = 0, LMP_STYLE_ATOM = 1, &
       LMP_STYLE_LOCAL = 2, LMP_TYPE_SCALAR = 0, LMP_TYPE_VECTOR = 1, &
       LMP_TYPE_ARRAY = 2, LMP_SIZE_VECTOR = 3, LMP_SIZE_ROWS = 4, &
       LMP_SIZE_COLS = 5, LAMMPS_INT = 0, LAMMPS_INT_2D = 1, &
       LAMMPS_DOUBLE = 2, LAMMPS_DOUBLE_2D = 3, LAMMPS_INT64 = 4, &
       LAMMPS_INT64_2D = 5, LAMMPS_STRING = 6

   PUBLIC :: LMP_STYLE_GLOBAL, LMP_STYLE_ATOM, LMP_STYLE_LOCAL, &
       LMP_TYPE_SCALAR, LMP_TYPE_VECTOR, LMP_TYPE_ARRAY, &
       LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS,  LAMMPS_INT, &
       LAMMPS_INT_2D, LAMMPS_DOUBLE, LAMMPS_DOUBLE_2D, LAMMPS_INT64, &
       LAMMPS_INT64_2D, LAMMPS_STRING 

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
         character (kind=C_char), dimension(*) :: file, str
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

      function lammps_version (ptr) result (version) &
      bind (C, name='lammps_version')
         import :: C_ptr, C_int
         type (C_ptr), value :: ptr
         integer (kind=C_int) :: version
      end function lammps_version

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

      subroutine lammps_actual_commands_list (ptr, ncmd, cmds) &
      bind (C, name='lammps_commands_list')
         import :: C_ptr, C_int
         type (C_ptr), value :: ptr
         type (C_ptr), dimension(*) :: cmds
         integer (C_int), value :: ncmd
      end subroutine lammps_actual_commands_list

      subroutine lammps_actual_commands_string (ptr, str) &
      bind (C, name='lammps_commands_string')
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: str
      end subroutine lammps_actual_commands_string

      subroutine lammps_free (ptr) bind (C, name='lammps_free')
         import :: C_ptr
         type (C_ptr), value :: ptr
      end subroutine lammps_free

      function lammps_actual_extract_setting (ptr, name) result (setting) &
      bind(C, name='lammps_extract_setting')
         import :: C_ptr, C_char, C_int
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name
         integer (kind=C_int) :: setting
      end function lammps_actual_extract_setting

      function lammps_actual_extract_global (ptr, name) &
      bind (C, name='lammps_extract_global') result (global)
         import :: C_ptr, C_char
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name
         type (C_ptr) :: global
      end function lammps_actual_extract_global

      subroutine lammps_actual_extract_box (ptr, boxlo, boxhi, xy, yz, &
            xz, periodicity, box_change) bind (C, name='lammps_extract_box')
         import :: C_ptr, C_double, C_int
         type (C_ptr), value :: ptr
         real (C_double) :: boxlo(3), boxhi(3), xy, yz, xz
         integer (C_int) :: periodicity(3), box_change
      end subroutine lammps_actual_extract_box

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

      function lammps_actual_get_thermo (ptr, name) result (dval) &
      bind (C, name='lammps_get_thermo')
         import :: C_ptr, C_char, C_double
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name
         real (C_double) :: dval
      end function lammps_actual_get_thermo

      function lammps_get_natoms (ptr) result (natoms) &
      bind (C, name='lammps_get_natoms')
         import :: C_ptr, C_double
         type (C_ptr), value :: ptr
         real (C_double) :: natoms
      end function lammps_get_natoms

      function lammps_actual_set_variable (ptr, name, str) result (err) &
      bind (C, name='lammps_set_variable')
         import :: C_ptr, C_char, C_int
         type (C_ptr), value :: ptr
         character (kind=C_char), dimension(*) :: name, str
         integer (C_int) :: err
      end function lammps_actual_set_variable

      subroutine lammps_actual_reset_box (ptr, boxlo, boxhi, xy, yz, xz) &
      bind (C, name='lammps_reset_box')
         import :: C_ptr, C_double, C_int
         type (C_ptr), value :: ptr
         real (C_double) :: boxlo(3), boxhi(3), xy, yz, xz
      end subroutine lammps_actual_reset_box

      subroutine lammps_actual_gather_atoms (ptr, name, type, count, data) &
      bind (C, name='lammps_gather_atoms')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count
      end subroutine lammps_actual_gather_atoms

      subroutine lammps_actual_gather_atoms_concat (ptr, name, type, count, &
      data) bind (C, name='lammps_gather_atoms_concat')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count
      end subroutine lammps_actual_gather_atoms_concat

      subroutine lammps_actual_gather_atoms_subset (ptr, name, type, count, &
      ndata, ids, data) bind (C, name='lammps_gather_atoms_subset')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count, ndata
         integer (C_int) :: ids(:)
      end subroutine lammps_actual_gather_atoms_subset

      subroutine lammps_actual_scatter_atoms (ptr, name, type, count, data) &
      bind (C, name='lammps_scatter_atoms')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count
      end subroutine lammps_actual_scatter_atoms

      subroutine lammps_actual_scatter_atoms_subset (ptr, name, type, count, &
      ndata, ids, data) bind (C, name='lammps_scatter_atoms_subset')
         import :: C_ptr, C_int, C_char
         type (C_ptr), value :: ptr, data
         character (kind=C_char), dimension(*) :: name
         integer (C_int), value :: type, count, ndata
         integer (C_int), dimension(*) :: ids
      end subroutine lammps_actual_scatter_atoms_subset

      function lammps_actual_config_has_package (package_name) &
      result (has_it) bind (C, name='lammps_config_has_package')
         import :: C_char, C_int
         character (len=1, kind=C_char), dimension(*) :: package_name
         integer (C_int) :: has_it
      end function lammps_actual_config_has_package

      function lammps_config_package_count () result (count) &
      bind (C, name='lammps_config_package_count')
         import :: C_int
         integer (C_int) :: count
      end function lammps_config_package_count

      function lammps_actual_config_package_name (index, buffer, max_size) &
      result (num) bind (C, name='lammps_config_package_name')
         import :: C_int, C_char
         integer (C_int), value :: index, max_size
         character (len=1, kind=C_char), dimension(*) :: buffer
         integer (C_int) :: num
      end function lammps_actual_config_package_name

      function lammps_actual_config_has_gzip_support () result (C_has_it) &
      bind (C, name='lammps_config_has_gzip_support')
         import :: C_int
         integer (C_int) :: C_has_it
      end function lammps_actual_config_has_gzip_support

      function lammps_actual_config_has_png_support () result (C_has_it) &
      bind (C, name='lammps_config_has_png_support')
         import :: C_int
         integer (C_int) :: C_has_it
      end function lammps_actual_config_has_png_support

      function lammps_actual_config_has_jpeg_support () result (C_has_it) &
      bind (C, name='lammps_config_has_jpeg_support')
         import :: C_int
         integer (C_int) :: C_has_it
      end function lammps_actual_config_has_jpeg_support

      function lammps_actual_config_has_ffmpeg_support () result (C_has_it) &
      bind (C, name='lammps_config_has_ffmpeg_support')
         import :: C_int
         integer (C_int) :: C_has_it
      end function lammps_actual_config_has_ffmpeg_support

      function lammps_actual_config_has_exceptions () result (C_has_it) &
      bind (C, name='lammps_config_has_exceptions')
         import :: C_int
         integer (C_int) :: C_has_it
      end function lammps_actual_config_has_exceptions

      function lammps_actual_find_pair_neighlist (ptr, style, exact, nsub, &
            request) result (C_neighlist) &
      bind (C, name='lammps_find_pair_neighlist')
         import :: C_ptr, C_int, C_char
         integer (C_int) :: C_neighlist
         type (C_ptr), value :: ptr
         character (len=1, kind=C_char), dimension(*) :: style
         integer (C_int), value :: exact, nsub, request
      end function lammps_actual_find_pair_neighlist

      function lammps_actual_find_fix_neighlist (ptr, id, request) &
      result (C_neighlist) bind (C, name='lammps_find_fix_neighlist')
         import :: C_ptr, C_int, C_char
         integer (C_int) :: C_neighlist
         type (C_ptr), value :: ptr
         character (len=1, kind=C_char), dimension(*) :: id
         integer (C_int), value :: request
      end function lammps_actual_find_fix_neighlist

      function lammps_actual_find_compute_neighlist (ptr, id, request) &
      result (C_neighlist) bind (C, name='lammps_find_compute_neighlist')
         import :: C_ptr, C_int, C_char
         integer (C_int) :: C_neighlist
         type (C_ptr), value :: ptr
         character (len=1, kind=C_char), dimension(*) :: id
         integer (C_int), value :: request
      end function lammps_actual_find_compute_neighlist

      function lammps_actual_neighlist_num_elements (ptr, idx) &
      result (nelements) bind (C, name='lammps_neighlist_num_elements')
         import :: C_ptr, C_int
         integer (C_int) :: nelements
         type (C_ptr), value :: ptr
         integer (C_int), value :: idx
      end function lammps_actual_neighlist_num_elements

      subroutine lammps_actual_neighlist_element_neighbors (ptr, idx, &
         element, iatom, numneigh, neighbors) &
      bind (C, name='lammps_neighlist_element_neighbors')
         import :: C_ptr, C_int
         type (C_ptr), value :: ptr
         integer (C_int), value :: idx, element
         integer (C_int) :: iatom, numneigh
         type (C_ptr) :: neighbors
      end subroutine lammps_actual_neighlist_element_neighbors

      subroutine lammps_actual_create_atoms (ptr, n, id, type, x, v, image, &
      shrinkexceed) bind (C, name='lammps_create_atoms')
         import :: C_ptr, C_int64_t, C_double, C_int
         type (C_ptr), value :: ptr
         integer (C_int), value :: n
         integer (C_tagint) :: id
         integer (C_int) :: type
         real (C_double), dimension(*) :: x, v
         integer (C_imageint), dimension(*) :: image
         integer (C_int), value :: shrinkexceed
      end subroutine lammps_actual_create_atoms

#ifdef LAMMPS_EXCEPTIONS

      function lammps_actual_has_error (ptr) result (C_has_it) &
      bind (C, name='lammps_has_error')
         import :: C_int, C_ptr
         type (C_ptr), value :: ptr
         integer (C_int) :: C_has_it
      end function lammps_actual_has_error

      function lammps_actual_get_last_error_message (ptr, buffer, &
      buffer_size) result (error_type) &
      bind (C, name='lammps_get_last_error_message')
         import :: C_ptr, C_char, C_int
         type (C_ptr), value :: ptr
         integer (C_int), value :: buffer_size
         character (len=1, kind=C_char), dimension(*) :: buffer
         integer (C_int) :: error_type
      end function lammps_actual_get_last_error_message

#endif

   end interface

   ! Generic functions for the wrappers below {{{1

   interface lammps_extract_global
      module procedure lammps_extract_global_i, &
         lammps_extract_global_dp
   end interface lammps_extract_global

   interface lammps_extract_atom
      module procedure lammps_extract_atom_ia, &
         lammps_extract_atom_dpa, &
         lammps_extract_atom_dp2a
   end interface lammps_extract_atom

   interface lammps_extract_compute
      module procedure lammps_extract_compute_dp, &
         lammps_extract_compute_dpa, &
         lammps_extract_compute_dp2a
   end interface lammps_extract_compute

   interface lammps_extract_fix
      module procedure lammps_extract_fix_dp, &
         lammps_extract_fix_dpa, &
         lammps_extract_fix_dp2a
   end interface lammps_extract_fix

   interface lammps_extract_variable
      module procedure lammps_extract_variable_dp, &
         lammps_extract_variable_dpa
   end interface lammps_extract_variable

   interface lammps_gather_atoms
      module procedure lammps_gather_atoms_ia, lammps_gather_atoms_dpa
   end interface lammps_gather_atoms

   interface lammps_gather_atoms_concat
      module procedure lammps_gather_atoms_concat_ia, &
         lammps_gather_atoms_concat_dpa
   end interface lammps_gather_atoms_concat

   interface lammps_gather_atoms_subset
      module procedure lammps_gather_atoms_subset_ia, &
         lammps_gather_atoms_subset_dpa
   end interface lammps_gather_atoms_subset

   interface lammps_scatter_atoms
      module procedure lammps_scatter_atoms_ia, lammps_scatter_atoms_dpa
   end interface lammps_scatter_atoms

   interface lammps_scatter_atoms_subset
      module procedure lammps_scatter_atoms_subset_ia, &
         lammps_scatter_atoms_subset_dpa
   end interface lammps_scatter_atoms_subset

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

   subroutine lammps_commands_list (ptr, cmds)
      type (C_ptr), intent(in) :: ptr
      character (len=*), dimension(:) :: cmds
      integer (C_int) :: ncmd
!      character (kind=C_char,len=1), dimension(size(cmds)) :: C_cmds
      type (C_ptr), dimension(:), allocatable :: C_cmds
      character (len=1, kind=C_char), allocatable, target :: C_strings(:,:)
      integer :: i, max_len
      ncmd = size(cmds)
      allocate (C_cmds(ncmd))
      max_len = 0
      do i=1, size(cmds)
         if ( len(cmds(i)) > max_len ) max_len = len(cmds(i))
      end do
      allocate (C_strings(max_len + 1, ncmd))
      do i=1, size(cmds)
         C_strings(:,i) = string2Cstring(cmds(i))
         C_cmds(i) = C_loc(C_strings(1,i))
      end do
      call lammps_actual_commands_list (ptr, ncmd, C_cmds)
      deallocate (C_strings)
      deallocate (C_cmds)
   end subroutine lammps_commands_list

!-----------------------------------------------------------------------------

   subroutine lammps_commands_string (ptr, str)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: str
      character (kind=C_char), dimension(len_trim(str)+1) :: C_str
      C_str = string2Cstring (str)
      call lammps_actual_commands_string (ptr, C_str)
   end subroutine lammps_commands_string

!-----------------------------------------------------------------------------

   function lammps_extract_setting (ptr, name) result (setting)
      integer :: setting
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (kind=C_char), dimension(len_trim(name)+1) :: C_name
      C_name = string2Cstring (name)
      setting = lammps_actual_extract_setting (ptr, C_name)
   end function lammps_extract_setting

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
      integer (C_int), pointer, intent(out) :: global
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      Cptr = lammps_extract_global_Cptr (ptr, name)
      call C_F_pointer (Cptr, global)
   end subroutine lammps_extract_global_i
   subroutine lammps_extract_global_dp (global, ptr, name)
      real (C_double), pointer, intent(out) :: global
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      Cptr = lammps_extract_global_Cptr (ptr, name)
      call C_F_pointer (Cptr, global)
   end subroutine lammps_extract_global_dp

!-----------------------------------------------------------------------------

! lammps_extract_box {{{2
   subroutine lammps_extract_box (ptr, boxlo, boxhi, xy, yz, xz, &
         periodicity, box_change)
      type (C_ptr), intent(in) :: ptr
      double precision, dimension(3), intent(out) :: boxlo, boxhi
      double precision, intent(out) :: xy, yz, xz
      logical, intent(out) :: periodicity(3), box_change
      integer (kind=C_int) :: C_periodicity(3), C_box_change
      real (C_double) :: C_boxlo(3), C_boxhi(3), C_xy, C_yz, C_xz
      call lammps_actual_extract_box (ptr, C_boxlo, C_boxhi, C_xy, C_yz, &
         C_xz, C_periodicity, C_box_change)
      boxlo = C_boxlo
      boxhi = C_boxhi
      xy = C_xy
      yz = C_yz
      xz = C_xz
      periodicity = (C_periodicity == 1) 
      box_change = (C_box_change == 1)
   end subroutine

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
      integer (C_int), dimension(:), pointer, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      integer (C_int), pointer :: nelements
      call lammps_extract_global_i (nelements, ptr, 'nlocal')
      Cptr = lammps_extract_atom_Cptr (ptr, name)
      call C_F_pointer (Cptr, atom, (/nelements/))
   end subroutine lammps_extract_atom_ia
   subroutine lammps_extract_atom_dpa (atom, ptr, name)
      real (C_double), dimension(:), pointer, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      integer (C_int), pointer :: nlocal
      integer :: nelements
      real (C_double), dimension(:), pointer :: Fptr
      if ( name == 'mass' ) then
         nelements = lammps_get_ntypes (ptr) + 1
      else if ( name == 'x' .or. name == 'v' .or. name == 'f' .or. &
                name == 'mu' .or. name == 'omega' .or. name == 'torque' .or. &
                name == 'angmom' ) then
         ! We should not be getting a rank-2 array here!
         call lammps_error_all (ptr, FLERR, 'You cannot extract those atom&
            & data (' // trim(name) // ') into a rank 1 array.')
         return
      else
         ! Everything else we can get is probably nlocal units long
         call lammps_extract_global_i (nlocal, ptr, 'nlocal')
         nelements = nlocal
      end if
      Cptr = lammps_extract_atom_Cptr (ptr, name)
      call C_F_pointer (Cptr, Fptr, (/nelements/))
      if ( name == 'mass' ) then
         atom(0:) => Fptr
      else
         atom => Fptr
      end if
   end subroutine lammps_extract_atom_dpa
   subroutine lammps_extract_atom_dp2a (atom, ptr, name)
      real (C_double), dimension(:,:), pointer, intent(out) :: atom
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      type (C_ptr) :: Cptr
      type (C_ptr), pointer, dimension(:) :: Catom
      integer (C_int), pointer :: nelements
      if ( name /= 'x' .and. name /= 'v' .and. name /= 'f' .and. &
           name /= 'mu' .and. name /= 'omega' .and. name /= 'tandque' .and. &
           name /= 'angmom' ) then
         ! We should not be getting a rank-2 array here!
         call lammps_error_all (ptr, FLERR, 'You cannot extract those atom&
            & data (' // trim(name) // ') into a rank 2 array.')
         return
      end if
      Cptr = lammps_extract_atom_Cptr (ptr, name)
      call lammps_extract_global_i (nelements, ptr, 'nlocal')
      ! Catom will now be the array of void* pointers that the void** pointer
      ! pointed to.  Catom(1) is now the pointer to the first element.
      call C_F_pointer (Cptr, Catom, (/nelements/))
      ! Now get the actual array, which has its shape transposed from what we
      ! might think of it in C
      call C_F_pointer (Catom(1), atom, (/3, nelements/))
   end subroutine lammps_extract_atom_dp2a

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
      real (C_double), pointer, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      type (C_ptr) :: Cptr
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
      call C_F_pointer (Cptr, compute)
   end subroutine lammps_extract_compute_dp
   subroutine lammps_extract_compute_dpa (compute, ptr, id, style, type)
      real (C_double), dimension(:), pointer, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      type (C_ptr) :: Cptr
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
      Cptr = lammps_extract_compute_Cptr (ptr, id, style, type)
      call C_F_pointer (Cptr, compute, (/nelements/))
   end subroutine lammps_extract_compute_dpa
   subroutine lammps_extract_compute_dp2a (compute, ptr, id, style, type)
      real (C_double), dimension(:,:), pointer, intent(out) :: compute
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type
      type (C_ptr) :: Cptr
      type (C_ptr), pointer, dimension(:) :: Ccompute
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
      Cptr = lammps_extract_compute_Cptr (ptr, id, style, type)
      call C_F_pointer (Cptr, Ccompute, (/nr/))
      ! Note that the matrix is transposed, from Fortran's perspective
      call C_F_pointer (Ccompute(1), compute, (/nc, nr/))
   end subroutine lammps_extract_compute_dp2a

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
      real (C_double), intent(out) :: fix
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
            call lammps_error_all (ptr, FLERR, 'Invalid extract_fix style/&
               &type combination.')
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
   subroutine lammps_extract_fix_dpa (fix, ptr, id, style, type, i, j)
      real (C_double), dimension(:), pointer, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      type (C_ptr) :: Cptr
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
      call C_F_pointer (Cptr, fix, (/fix_len/))
      ! Memory is only allocated for "global" fix variables, which we should
      ! never get here, so no need to call lammps_free!
   end subroutine lammps_extract_fix_dpa
   subroutine lammps_extract_fix_dp2a (fix, ptr, id, style, type, i, j)
      real (C_double), dimension(:,:), pointer, intent(out) :: fix
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: style, type, i, j
      type (C_ptr) :: Cptr
      type (C_ptr), pointer, dimension(:) :: Cfix
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
      ! Extract pointer to first element as Cfix(1)
      call C_F_pointer (Cptr, Cfix, (/nr/))
      ! Now extract the array, which is transposed
      call C_F_pointer (Cfix(1), fix, (/nc, nr/))
   end subroutine lammps_extract_fix_dp2a

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
   subroutine lammps_extract_variable_dp (variable, ptr, name, group)
      real (C_double), intent(out) :: variable
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: group
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
   subroutine lammps_extract_variable_dpa (variable, ptr, name, group)
      real (C_double), dimension(:), allocatable, intent(out) :: variable
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

!-------------------------------------------------------------------------2}}}

   function lammps_get_thermo (ptr, name) result (dval)
      double precision :: dval
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      Cname = string2Cstring (name)
      dval = lammps_actual_get_thermo (ptr, Cname)
   end function lammps_get_thermo

   subroutine lammps_set_variable (ptr, name, str, err)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name, str
      integer, optional :: err
      integer (C_int) :: Cerr
      character (kind=C_char) :: Cname(len_trim(name)+1), Cstr(len_trim(str)+1)
      Cname = string2Cstring (name)
      Cstr = string2Cstring (str)
      Cerr = lammps_actual_set_variable (ptr, Cname, Cstr)
      if ( present(err) ) err = Cerr
   end subroutine lammps_set_variable

   subroutine lammps_reset_box (ptr, boxlo, boxhi, xy, yz, xz)
      type (C_ptr), intent(in) :: ptr
      double precision, dimension(3), intent(in) :: boxlo, boxhi
      double precision, intent(in) :: xy, yz, xz
      real (C_double) :: C_boxlo(3), C_boxhi(3), C_xy, C_yz, C_xz
      C_boxlo = boxlo
      C_boxhi = boxhi
      C_xy = xy
      C_xz = xz
      C_yz = yz
      call lammps_actual_reset_box (ptr, C_boxlo, C_boxhi, C_xy, C_xz, C_yz)
   end subroutine lammps_reset_box

! lammps_gather_atoms {{{2
   subroutine lammps_gather_atoms_ia (ptr, name, count, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count
      integer, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      integer (C_int), dimension(:), pointer :: Fdata
      integer (C_int) :: natoms
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), parameter :: Ctype = 0_C_int
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
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), parameter :: Ctype = 1_C_int
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

! lammps_gather_atoms_concat {{{2
   subroutine lammps_gather_atoms_concat_ia (ptr, name, count, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count
      integer, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      integer (C_int), dimension(:), pointer :: Fdata
      integer (C_int) :: natoms
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), parameter :: Ctype = 0_C_int
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
      call lammps_actual_gather_atoms_concat (ptr, Cname, Ctype, Ccount, Cdata)
      data = Fdata
      deallocate (Fdata)
   end subroutine lammps_gather_atoms_concat_ia
   subroutine lammps_gather_atoms_concat_dpa (ptr, name, count, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count
      double precision, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      real (C_double), dimension(:), pointer :: Fdata
      integer (C_int) :: natoms
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), parameter :: Ctype = 1_C_int
      integer (C_int) :: Ccount
      natoms = lammps_get_natoms (ptr)
      Cname = string2Cstring (name)
      if ( count /= 1 .and. count /= 3 ) then
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms_concat&
            & requires count to be either 1 or 3')
      else
         Ccount = count
      end if
      allocate ( Fdata(count*natoms) )
      allocate ( data(count*natoms) )
      Cdata = C_loc (Fdata(1))
      call lammps_actual_gather_atoms_concat (ptr, Cname, Ctype, Ccount, Cdata)
      data = Fdata(:)
      deallocate (Fdata)
   end subroutine lammps_gather_atoms_concat_dpa

!-----------------------------------------------------------------------------

! lammps_gather_atoms_subset {{{2
   subroutine lammps_gather_atoms_subset_ia (ptr,name,count,ids,data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count, ids(:)
      integer, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      integer (C_int), dimension(:), pointer :: Fdata
      integer (C_int) :: ndata, Cids(size(ids))
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), parameter :: Ctype = 0_C_int
      integer (C_int) :: Ccount
      ndata = size(ids)
      Cname = string2Cstring (name)
      Cids = ids
      if ( count /= 1 .and. count /= 3 ) then
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms_subset&
            & requires count to be either 1 or 3')
      else
         Ccount = count
      end if
      allocate ( Fdata(count*ndata) )
      allocate ( data(count*ndata) )
      Cdata = C_loc (Fdata(1))
      call lammps_actual_gather_atoms_subset (ptr, Cname, Ctype, Ccount, &
         ndata, Cids, Cdata)
      data = Fdata
      deallocate (Fdata)
   end subroutine lammps_gather_atoms_subset_ia
   subroutine lammps_gather_atoms_subset_dpa (ptr,name,count,ids,data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, intent(in) :: count, ids(:)
      double precision, dimension(:), allocatable, intent(out) :: data
      type (C_ptr) :: Cdata
      real (C_double), dimension(:), pointer :: Fdata
      integer (C_int) :: ndata, Cids(size(ids))
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), parameter :: Ctype = 1_C_int
      integer (C_int) :: Ccount
      ndata = size(ids)
      Cname = string2Cstring (name)
      Cids = ids
      if ( count /= 1 .and. count /= 3 ) then
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      else
         Ccount = count
      end if
      allocate ( Fdata(count*ndata) )
      allocate ( data(count*ndata) )
      Cdata = C_loc (Fdata(1))
      call lammps_actual_gather_atoms_subset (ptr, Cname, Ctype, Ccount, &
         ndata, Cids, Cdata)
      data = Fdata
      deallocate (Fdata)
   end subroutine lammps_gather_atoms_subset_dpa

!-----------------------------------------------------------------------------

! lammps_scatter_atoms {{{2
   subroutine lammps_scatter_atoms_ia (ptr, name, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: data
      integer (kind=C_int) :: natoms, Ccount
      integer (kind=C_int), parameter :: Ctype = 0_C_int
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), dimension(size(data)), target :: Fdata
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
      integer (kind=C_int), parameter :: Ctype = 1_C_int
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      real (C_double), dimension(size(data)), target :: Fdata
      type (C_ptr) :: Cdata
      natoms = lammps_get_natoms (ptr)
      Cname = string2Cstring (name)
      Ccount = size(data) / natoms
      if ( Ccount /= 1 .and. Ccount /= 3 ) &
         call lammps_error_all (ptr, FLERR, 'lammps_scatter_atoms requires&
            & count to be either 1 or 3')
      Fdata = data
      Cdata = C_loc (Fdata(1))
      call lammps_actual_scatter_atoms (ptr, Cname, Ctype, Ccount, Cdata)
   end subroutine lammps_scatter_atoms_dpa

!-----------------------------------------------------------------------------

! lammps_scatter_atoms_subset {{{2
   subroutine lammps_scatter_atoms_subset_ia (ptr, name, ids, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: data, ids
      integer (kind=C_int) :: ndata, Ccount, Cids(size(ids))
      integer (kind=C_int), parameter :: Ctype = 0_C_int
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      integer (C_int), dimension(size(data)), target :: Fdata
      type (C_ptr) :: Cdata
      ndata = size(ids)
      Cname = string2Cstring (name)
      Ccount = size(data) / ndata
      if ( Ccount /= 1 .and. Ccount /= 3 ) &
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      Fdata = data
      Cdata = C_loc (Fdata(1))
      Cids = ids
      call lammps_actual_scatter_atoms_subset (ptr, Cname, Ctype, Ccount, &
         ndata, Cids, Cdata)
   end subroutine lammps_scatter_atoms_subset_ia
   subroutine lammps_scatter_atoms_subset_dpa (ptr, name, ids, data)
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: name
      double precision, dimension(:), intent(in) :: data, ids
      integer (kind=C_int) :: ndata, Ccount, Cids(size(ids))
      integer (kind=C_int), parameter :: Ctype = 1_C_int
      character (kind=C_char), dimension(len_trim(name)+1) :: Cname
      real (C_double), dimension(size(data)), target :: Fdata
      type (C_ptr) :: Cdata
      ndata = size(ids)
      Cname = string2Cstring (name)
      Ccount = size(data) / ndata
      if ( Ccount /= 1 .and. Ccount /= 3 ) &
         call lammps_error_all (ptr, FLERR, 'lammps_gather_atoms requires&
            & count to be either 1 or 3')
      Fdata = data
      Cdata = C_loc (Fdata(1))
      call lammps_actual_scatter_atoms_subset (ptr, Cname, Ctype, Ccount, &
         ndata, Cids, Cdata)
   end subroutine lammps_scatter_atoms_subset_dpa

!-----------------------------------------------------------------------------

   subroutine lammps_create_atoms (ptr, id, type, x, v, image, shrinkexceed)
      type (C_ptr), intent(in) :: ptr
      integer (kind=C_tagint), dimension(:), optional :: id
      integer, dimension(:) :: type
      double precision, dimension(:,:) :: x
      double precision, dimension(:,:), optional :: v
      integer (kind=C_imageint), dimension(:), target, optional :: image
      logical, optional :: shrinkexceed
      real (C_double), dimension(size(x)) :: C_x, C_v
      integer (C_int) :: C_shrinkexceed, n
      integer (kind=C_tagint) :: C_id
      integer (C_imageint), dimension(size(x)/3) :: C_image
      integer (C_int) :: C_type
      if (shrinkexceed) then
         C_shrinkexceed = 1_C_int
      else
         C_shrinkexceed = 0_C_int
      end if
      C_x = reshape(transpose(x), shape(C_x))
      if ( present(v) ) then
         C_v = reshape(transpose(v), shape(C_v))
      else
         C_v = 0.0_C_double
      end if
      if ( present(image) ) then
         C_image = image
      else
         C_image = int(0,kind=C_imageint)
      end if
      n = size(type)
      call lammps_actual_create_atoms (ptr, n, C_id, C_type, C_x, C_v, &
         C_image, C_shrinkexceed)
   end subroutine lammps_create_atoms

!-----------------------------------------------------------------------------

   function lammps_config_has_package (package_name) result (has_it)
      character (len=*), intent(in) :: package_name
      character (len=1, kind=C_char), dimension(len_trim(package_name)+1) :: &
         C_package_name
      logical :: has_it
      integer (C_int) :: C_has_it
      C_package_name = string2Cstring (package_name)
      C_has_it = lammps_actual_config_has_package (C_package_name)
      has_it = (C_has_it == 1)
   end function lammps_config_has_package

!-----------------------------------------------------------------------------

   function lammps_config_package_name (index, buffer) result (installed)
      character (len=*), intent(inout) :: buffer
      integer, intent(in) :: index
      logical :: installed
      integer (kind=C_int) :: C_installed, C_index, max_size
      character (len=1, kind=C_char), dimension(len_trim(buffer)+1) :: C_buffer
      C_buffer = string2Cstring (buffer)
      max_size = len(buffer)
      C_index = index
      C_installed = lammps_actual_config_package_name (C_index, C_buffer, &
         max_size)
      installed = (C_installed == 1_C_int)
      buffer = Cstring2string (C_buffer)
   end function lammps_config_package_name

!-----------------------------------------------------------------------------

   logical function lammps_config_has_gzip_support () result (has_it)
      integer (C_int) :: C_has_it
      C_has_it = lammps_actual_config_has_gzip_support ()
      has_it = (C_has_it == 1_C_int)
   end function lammps_config_has_gzip_support

!-----------------------------------------------------------------------------

   logical function lammps_config_has_png_support () result (has_it)
      integer (C_int) :: C_has_it
      C_has_it = lammps_actual_config_has_png_support ()
      has_it = (C_has_it == 1_C_int)
   end function lammps_config_has_png_support

!-----------------------------------------------------------------------------

   logical function lammps_config_has_jpeg_support () result (has_it)
      integer (C_int) :: C_has_it
      C_has_it = lammps_actual_config_has_jpeg_support ()
      has_it = (C_has_it == 1_C_int)
   end function lammps_config_has_jpeg_support

!-----------------------------------------------------------------------------

   logical function lammps_config_has_ffmpeg_support () result (has_it)
      integer (C_int) :: C_has_it
      C_has_it = lammps_actual_config_has_ffmpeg_support ()
      has_it = (C_has_it == 1_C_int)
   end function lammps_config_has_ffmpeg_support

!-----------------------------------------------------------------------------

   logical function lammps_config_has_exceptions () result (has_it)
      integer (C_int) :: C_has_it
      C_has_it = lammps_actual_config_has_exceptions ()
      has_it = (C_has_it == 1_C_int)
   end function lammps_config_has_exceptions

!-----------------------------------------------------------------------------

   function lammps_find_pair_neighlist (ptr, style, exact, nsub, request) &
   result (neighlist)
      integer :: neighlist
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: style
      logical, intent(in) :: exact
      integer, intent(in) :: nsub, request
      integer (C_int) :: C_exact, C_nsub, C_neighlist, C_request
      character (kind=C_char), dimension(len_trim(style)+1) :: C_style
      if (exact) then
         C_exact = 1_C_int
      else
         C_exact = 0_C_int
      end if
      C_nsub = nsub
      C_request = request
      C_style = string2Cstring (style)
      C_neighlist = lammps_actual_find_pair_neighlist (ptr, C_style, C_exact, &
         C_nsub, C_request)
      neighlist = C_neighlist
   end function lammps_find_pair_neighlist

!-----------------------------------------------------------------------------

   function lammps_find_fix_neighlist (ptr, id, request) result (neighlist)
      integer :: neighlist
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: request
      integer (C_int) :: C_request, C_neighlist
      character (kind=C_char), dimension(len_trim(id)+1) :: C_id
      C_id = string2Cstring(id)
      C_request = request
      C_neighlist = lammps_actual_find_fix_neighlist (ptr, C_id, C_request)
      neighlist = C_neighlist
   end function lammps_find_fix_neighlist

!-----------------------------------------------------------------------------

   function lammps_find_compute_neighlist (ptr, id, request) result (neighlist)
      integer :: neighlist
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(in) :: id
      integer, intent(in) :: request
      integer (C_int) :: C_request
      character (kind=C_char), dimension(len_trim(id)+1) :: C_id
      C_id = string2Cstring(id)
      C_request = request
      neighlist = lammps_actual_find_compute_neighlist (ptr, C_id, C_request)
   end function lammps_find_compute_neighlist

!-----------------------------------------------------------------------------

   function lammps_neighlist_num_elements (ptr, idx) result (nelements)
      integer :: nelements
      type (C_ptr), intent(in) :: ptr
      integer, intent(in) :: idx
      integer (C_int) :: C_idx
      C_idx = idx
      nelements = lammps_actual_neighlist_num_elements (ptr, C_idx)
   end function lammps_neighlist_num_elements

!-----------------------------------------------------------------------------

   subroutine lammps_neighlist_element_neighbors (ptr, idx, element, iatom, &
         numneigh, neighbors)
      type (C_ptr), intent(in) :: ptr
      integer, intent(in) :: idx, element
      integer, intent(out) :: iatom, numneigh
      integer (C_int), dimension(:), pointer, intent(out) :: neighbors
      integer (C_int) :: C_idx, C_element, C_iatom, C_numneigh
      type (C_ptr) :: C_neighbors
      C_idx = idx
      C_element = element
      call lammps_actual_neighlist_element_neighbors (ptr, C_idx, C_element, &
         C_iatom, C_numneigh, C_neighbors)
      iatom = C_iatom
      numneigh = C_numneigh
      call C_F_pointer (C_neighbors, neighbors, [numneigh])
   end subroutine lammps_neighlist_element_neighbors

!-----------------------------------------------------------------------------

! These are only defined if -DLAMMPS_EXCEPTIONS was issued
#ifdef LAMMPS_EXCEPTIONS
   logical function lammps_has_error (ptr) result (has_it)
      type (C_ptr), intent(in) :: ptr
      integer (kind=C_int) :: C_has_it
      C_has_it = lammps_actual_has_error (ptr)
      has_it = (C_has_it == 1_C_int)
   end function lammps_has_error

!-----------------------------------------------------------------------------

   function lammps_get_last_error_message (ptr, buffer) result (error_type)
      integer (C_int) :: error_type
      type (C_ptr), intent(in) :: ptr
      character (len=*), intent(out) :: buffer
      integer (C_int) :: buffer_size
      character (len=1, kind=C_char), dimension(len(buffer)+1) :: C_buffer
      buffer_size = len(buffer)
      error_type = lammps_actual_get_last_error_message (ptr, C_buffer, &
         buffer_size)
      buffer = Cstring2string (C_buffer)
   end function lammps_get_last_error_message
#endif
!-------------------------------------------------------------------------2}}}

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

   pure function Cstrlen (Cstring) result (n)
      character (len=1, kind=C_char), dimension(:), intent(in) :: Cstring
      integer :: n, i
      n = size(Cstring)
      do i = 1, size(Cstring)
         if ( Cstring(i) == C_NULL_CHAR ) then
            n = i - 1
            return
         end if
      end do
   end function Cstrlen

!-----------------------------------------------------------------------------

   pure function Cstring2string (Cstring) result (string)
      !use, intrinsic :: ISO_C_binding, only : C_char, C_NULL_CHAR
      character (len=1, kind=C_char), intent(in) :: Cstring (:)
      character (len=Cstrlen(Cstring)) :: string
      integer :: i, n
      n = Cstrlen(Cstring)
      string = ''
      forall (i = 1:n)
         string(i:i) = Cstring(i)
      end forall
   end function Cstring2string

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

! 1}}}

end module LAMMPS

! vim: foldmethod=marker tabstop=3 softtabstop=3 shiftwidth=3 expandtab

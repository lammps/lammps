FUNCTION f_lammps_with_args() BIND(C, name="f_lammps_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_with_args
  CHARACTER(len=12), DIMENSION(12), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps', '-log', 'none', &
      '-echo','screen','-nocite','-var','zpos','1.5','-var','x','2']

  lmp = lammps(args)
  f_lammps_with_args = lmp%handle
END FUNCTION f_lammps_with_args

SUBROUTINE f_lammps_close() BIND(C, name="f_lammps_close")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_atom() BIND(C)
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp, big_input, cont_input, pair_input
   IMPLICIT NONE

   CALL lmp%commands_list(big_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(pair_input)
END SUBROUTINE f_lammps_setup_extract_atom

FUNCTION f_lammps_extract_atom_mass() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_atom_mass
   REAL(c_double), DIMENSION(:), POINTER :: mass => NULL()

   mass = lmp%extract_atom('mass')
   f_lammps_extract_atom_mass = mass(1)
END FUNCTION f_lammps_extract_atom_mass

FUNCTION f_lammps_extract_atom_tag_int(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   INTEGER(c_int) :: f_lammps_extract_atom_tag_int
   INTEGER(c_int), DIMENSION(:), POINTER :: tag => NULL()

   tag = lmp%extract_atom('id')
   f_lammps_extract_atom_tag_int = tag(i)
END FUNCTION f_lammps_extract_atom_tag_int

FUNCTION f_lammps_extract_atom_tag_int64(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int64_t
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int64_t), INTENT(IN), VALUE :: i
   INTEGER(c_int64_t) :: f_lammps_extract_atom_tag_int64
   INTEGER(c_int64_t), DIMENSION(:), POINTER :: tag => NULL()

   tag = lmp%extract_atom('id')
   f_lammps_extract_atom_tag_int64 = tag(i)
END FUNCTION f_lammps_extract_atom_tag_int64

FUNCTION f_lammps_extract_atom_type(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   INTEGER(c_int) :: f_lammps_extract_atom_type
   INTEGER(c_int), DIMENSION(:), POINTER :: atype => NULL()

   atype = lmp%extract_atom('type')
   f_lammps_extract_atom_type = atype(i)
END FUNCTION f_lammps_extract_atom_type

FUNCTION f_lammps_extract_atom_mask(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   INTEGER(c_int) :: f_lammps_extract_atom_mask
   INTEGER(c_int), DIMENSION(:), POINTER :: mask => NULL()

   mask = lmp%extract_atom('mask')
   f_lammps_extract_atom_mask = mask(i)
END FUNCTION f_lammps_extract_atom_mask

SUBROUTINE f_lammps_extract_atom_x(i, x) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double), DIMENSION(3) :: x
   REAL(c_double), DIMENSION(:,:), POINTER :: xptr => NULL()

   xptr = lmp%extract_atom('x')
   x = xptr(:,i)
END SUBROUTINE f_lammps_extract_atom_x

SUBROUTINE f_lammps_extract_atom_v(i, v) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double), DIMENSION(3) :: v
   REAL(c_double), DIMENSION(:,:), POINTER :: vptr => NULL()

   vptr = lmp%extract_atom('v')
   v = vptr(:,i)
END SUBROUTINE f_lammps_extract_atom_v

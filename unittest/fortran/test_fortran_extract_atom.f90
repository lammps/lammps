MODULE keepatom
  USE liblammps
  TYPE(LAMMPS) :: lmp
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: demo_input = &
      [ CHARACTER(len=40) ::                                &
      'region       box block 0 $x 0 3 0 4',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(LEN=40), DIMENSION(2), PARAMETER :: cont_input = &
      [ CHARACTER(len=40) ::                                &
      'create_atoms 1 single &',                            &
      ' 0.2 0.1 0.1' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: pair_input = &
      [ CHARACTER(LEN=40) ::                                &
      'pair_style lj/cut 2.5',                              &
      'pair_coeff 1 1 1.0 1.0',                             &
      'mass 1 2.0' ]
END MODULE keepatom

FUNCTION f_lammps_with_args() BIND(C, name="f_lammps_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepatom, ONLY: lmp
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
  USE keepatom, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_atom () BIND(C)
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp, demo_input, cont_input, pair_input
   IMPLICIT NONE

   CALL lmp%commands_list(demo_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(pair_input)
!   CALL lmp%command('run 0')
END SUBROUTINE f_lammps_setup_extract_atom

FUNCTION f_lammps_extract_atom_mass () BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   REAL(C_double) :: f_lammps_extract_atom_mass
   REAL(C_double), DIMENSION(:), POINTER :: mass => NULL()

   mass = lmp%extract_atom('mass')
   f_lammps_extract_atom_mass = mass(1)
END FUNCTION f_lammps_extract_atom_mass

FUNCTION f_lammps_extract_atom_tag_int (i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double, C_int
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int), INTENT(IN), VALUE :: i
   INTEGER(C_int) :: f_lammps_extract_atom_tag_int
   INTEGER(C_int), DIMENSION(:), POINTER :: tag => NULL()

   tag = lmp%extract_atom('id')
   f_lammps_extract_atom_tag_int = tag(i)
END FUNCTION f_lammps_extract_atom_tag_int

FUNCTION f_lammps_extract_atom_tag_int64 (i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double, C_int64_t
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int64_t), INTENT(IN), VALUE :: i
   INTEGER(C_int64_t) :: f_lammps_extract_atom_tag_int64
   INTEGER(C_int64_t), DIMENSION(:), POINTER :: tag => NULL()

   tag = lmp%extract_atom('id')
   f_lammps_extract_atom_tag_int64 = tag(i)
END FUNCTION f_lammps_extract_atom_tag_int64

FUNCTION f_lammps_extract_atom_type(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int), INTENT(IN), VALUE :: i
   INTEGER(C_int) :: f_lammps_extract_atom_type
   INTEGER(C_int), DIMENSION(:), POINTER :: atype => NULL()

   atype = lmp%extract_atom('type')
   f_lammps_extract_atom_type = atype(i)
END FUNCTION f_lammps_extract_atom_type

FUNCTION f_lammps_extract_atom_mask(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int), INTENT(IN), VALUE :: i
   INTEGER(C_int) :: f_lammps_extract_atom_mask
   INTEGER(C_int), DIMENSION(:), POINTER :: mask => NULL()

   mask = lmp%extract_atom('mask')
   f_lammps_extract_atom_mask = mask(i)
END FUNCTION f_lammps_extract_atom_mask

SUBROUTINE f_lammps_extract_atom_x (i, x) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double, C_int
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int), INTENT(IN), VALUE :: i
   REAL(C_double), DIMENSION(3) :: x
   REAL(C_double), DIMENSION(:,:), POINTER :: xptr => NULL()

   xptr = lmp%extract_atom('x')
   x = xptr(:,i)
END SUBROUTINE f_lammps_extract_atom_x

SUBROUTINE f_lammps_extract_atom_v (i, v) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double, C_int
   USE LIBLAMMPS
   USE keepatom, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int), INTENT(IN), VALUE :: i
   REAL(C_double), DIMENSION(3) :: v
   REAL(C_double), DIMENSION(:,:), POINTER :: vptr => NULL()

   vptr = lmp%extract_atom('v')
   v = vptr(:,i)
END SUBROUTINE f_lammps_extract_atom_v

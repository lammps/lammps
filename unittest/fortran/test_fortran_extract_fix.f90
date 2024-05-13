FUNCTION f_lammps_with_args() BIND(C)
  USE ISO_C_BINDING, ONLY: C_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE
  TYPE(C_ptr) :: f_lammps_with_args
  CHARACTER(len=12), DIMENSION(12), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps', '-log', 'none', &
      '-echo','screen','-nocite','-var','zpos','1.5','-var','x','2']

  lmp = lammps(args)
  f_lammps_with_args = lmp%handle
END FUNCTION f_lammps_with_args

SUBROUTINE f_lammps_close() BIND(C)
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = C_NULL_PTR
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_fix() BIND(C)
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp, big_input, cont_input, pair_input, more_input
   IMPLICIT NONE

   CALL lmp%commands_list(big_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(more_input)
   CALL lmp%commands_list(pair_input)
   CALL lmp%command("fix state all store/state 0 z") ! per-atom vector
   CALL lmp%command("fix move all move linear 0 0 0") ! for per-atom array
   CALL lmp%command("fix recenter all recenter 1 1 1") ! global scalar, vector
   CALL lmp%command("variable natoms equal count(all)")
   CALL lmp%command("variable ts equal step")
   CALL lmp%command("fix vec all vector 1 v_natoms v_ts") ! global array
   CALL lmp%command("run 1") ! must be 1, otherwise move/recenter won't happen
END SUBROUTINE f_lammps_setup_extract_fix

FUNCTION f_lammps_extract_fix_global_scalar() BIND(C) RESULT(scalar)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: scalar

   scalar = lmp%extract_fix("recenter", lmp%style%global, lmp%type%scalar)
END FUNCTION f_lammps_extract_fix_global_scalar

FUNCTION f_lammps_extract_fix_global_vector(i) BIND(C) RESULT(element)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double) :: element

   element = lmp%extract_fix("recenter", lmp%style%global, lmp%type%vector, i)
END FUNCTION f_lammps_extract_fix_global_vector

FUNCTION f_lammps_extract_fix_global_array(i,j) BIND(C) RESULT(element)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i, j
   REAL(c_double) :: element

   element = lmp%extract_fix("vec", lmp%style%global, lmp%type%array, i, j)
END FUNCTION f_lammps_extract_fix_global_array

FUNCTION f_lammps_extract_fix_peratom_vector(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double) :: f_lammps_extract_fix_peratom_vector
   REAL(c_double), DIMENSION(:), POINTER :: vector

   vector = lmp%extract_fix("state", lmp%style%atom, lmp%type%vector, -1, -1)
   f_lammps_extract_fix_peratom_vector = vector(i)
END FUNCTION f_lammps_extract_fix_peratom_vector

FUNCTION f_lammps_extract_fix_peratom_array(i,j) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i, j
   REAL(c_double) :: f_lammps_extract_fix_peratom_array
   REAL(c_double), DIMENSION(:,:), POINTER :: array

   array = lmp%extract_fix("move", lmp%style%atom, lmp%type%array, -1, -1)
   f_lammps_extract_fix_peratom_array = array(i,j)
END FUNCTION f_lammps_extract_fix_peratom_array

FUNCTION f_lammps_with_args() BIND(C, name="f_lammps_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepstuff,      ONLY: lmp
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

SUBROUTINE f_lammps_box_setup() BIND(C)
   USE liblammps
   USE keepstuff, ONLY : lmp, demo_input
   IMPLICIT NONE

   CALL lmp%commands_list(demo_input)
END SUBROUTINE f_lammps_box_setup

SUBROUTINE f_lammps_delete_everything() BIND(C)
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE

   CALL lmp%command("delete_atoms group all");
END SUBROUTINE f_lammps_delete_everything

FUNCTION f_lammps_extract_box_xlo() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_box_xlo
   REAL(c_double) :: boxdim(3)

   CALL lmp%extract_box(boxlo=boxdim)
   f_lammps_extract_box_xlo = boxdim(1)
END FUNCTION f_lammps_extract_box_xlo

FUNCTION f_lammps_extract_box_xhi() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_box_xhi
   REAL(c_double) :: boxdim(3)

   CALL lmp%extract_box(boxhi=boxdim)
   f_lammps_extract_box_xhi = boxdim(1)
END FUNCTION f_lammps_extract_box_xhi

FUNCTION f_lammps_extract_box_ylo() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_box_ylo
   REAL(c_double) :: boxdim(3)

   CALL lmp%extract_box(boxlo=boxdim)
   f_lammps_extract_box_ylo = boxdim(2)
END FUNCTION f_lammps_extract_box_ylo

FUNCTION f_lammps_extract_box_yhi() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_box_yhi
   REAL(c_double) :: boxdim(3)

   CALL lmp%extract_box(boxhi=boxdim)
   f_lammps_extract_box_yhi = boxdim(2)
END FUNCTION f_lammps_extract_box_yhi

FUNCTION f_lammps_extract_box_zlo() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_box_zlo
   REAL(c_double) :: boxdim(3)

   CALL lmp%extract_box(boxlo=boxdim)
   f_lammps_extract_box_zlo = boxdim(2)
END FUNCTION f_lammps_extract_box_zlo

FUNCTION f_lammps_extract_box_zhi() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_box_zhi
   REAL(c_double) :: boxdim(3)

   CALL lmp%extract_box(boxhi=boxdim)
   f_lammps_extract_box_zhi = boxdim(2)
END FUNCTION f_lammps_extract_box_zhi

SUBROUTINE f_lammps_reset_box_2x() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: newlo(3), newhi(3), xy, yz, xz

   xy = 0.0_c_double
   yz = 0.0_c_double
   xz = 0.0_c_double
   newlo = [-1.0_c_double, -1.0_c_double, -1.0_c_double]
   newhi = [3.0_c_double, 3.0_c_double, 3.0_c_double]
   CALL lmp%reset_box(newlo, newhi, xy, yz, xz)
END SUBROUTINE f_lammps_reset_box_2x

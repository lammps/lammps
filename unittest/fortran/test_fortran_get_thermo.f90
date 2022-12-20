FUNCTION f_lammps_with_args() BIND(C)
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

SUBROUTINE f_lammps_close() BIND(C)
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_get_thermo_setup() BIND(C)
   USE liblammps
   USE keepstuff, ONLY : lmp, big_input, cont_input, pair_input
   IMPLICIT NONE

   CALL lmp%commands_list(big_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(pair_input)
END SUBROUTINE f_lammps_get_thermo_setup

FUNCTION f_lammps_get_thermo_natoms() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_natoms

   f_lammps_get_thermo_natoms = lmp%get_thermo('atoms')
END FUNCTION f_lammps_get_thermo_natoms

FUNCTION f_lammps_get_thermo_dt() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_dt

   f_lammps_get_thermo_dt = lmp%get_thermo('dt')
END FUNCTION f_lammps_get_thermo_dt

FUNCTION f_lammps_get_thermo_vol() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_vol

   f_lammps_get_thermo_vol = lmp%get_thermo('vol')
END FUNCTION f_lammps_get_thermo_vol

FUNCTION f_lammps_get_thermo_lx() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_lx

   f_lammps_get_thermo_lx = lmp%get_thermo('lx')
END FUNCTION f_lammps_get_thermo_lx

FUNCTION f_lammps_get_thermo_ly() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_ly

   f_lammps_get_thermo_ly = lmp%get_thermo('ly')
END FUNCTION f_lammps_get_thermo_ly

FUNCTION f_lammps_get_thermo_lz() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_lz

   f_lammps_get_thermo_lz = lmp%get_thermo('lz')
END FUNCTION f_lammps_get_thermo_lz

FUNCTION f_lammps_get_thermo_xlo() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_xlo

   f_lammps_get_thermo_xlo = lmp%get_thermo('xlo')
END FUNCTION f_lammps_get_thermo_xlo

FUNCTION f_lammps_get_thermo_xhi() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_xhi

   f_lammps_get_thermo_xhi = lmp%get_thermo('xhi')
END FUNCTION f_lammps_get_thermo_xhi

FUNCTION f_lammps_get_thermo_ylo() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_ylo

   f_lammps_get_thermo_ylo = lmp%get_thermo('ylo')
END FUNCTION f_lammps_get_thermo_ylo

FUNCTION f_lammps_get_thermo_yhi() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_yhi

   f_lammps_get_thermo_yhi = lmp%get_thermo('yhi')
END FUNCTION f_lammps_get_thermo_yhi

FUNCTION f_lammps_get_thermo_zlo() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_zlo

   f_lammps_get_thermo_zlo = lmp%get_thermo('zlo')
END FUNCTION f_lammps_get_thermo_zlo

FUNCTION f_lammps_get_thermo_zhi() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_get_thermo_zhi

   f_lammps_get_thermo_zhi = lmp%get_thermo('zhi')
END FUNCTION f_lammps_get_thermo_zhi

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

SUBROUTINE f_lammps_last_thermo_setup() BIND(C)
  USE liblammps
  USE keepstuff, ONLY : lmp, big_input, cont_input, pair_input
  IMPLICIT NONE

  CALL lmp%commands_list(big_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('thermo 10')
  CALL lmp%command('run 15 post no')
END SUBROUTINE f_lammps_last_thermo_setup

FUNCTION f_lammps_last_thermo_step() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_int64_t
  USE liblammps
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_last_thermo_step
  INTEGER :: size_bigint
  INTEGER(c_int), POINTER :: ival
  INTEGER(c_int64_t), POINTER :: bval

  size_bigint = lmp%extract_setting('bigint')
  IF (size_bigint == 4) THEN
      ival = lmp%last_thermo('step',1)
      f_lammps_last_thermo_step = INT(ival)
  ELSE
      bval = lmp%last_thermo('step',1)
      f_lammps_last_thermo_step = INT(bval)
  END IF
END FUNCTION f_lammps_last_thermo_step

FUNCTION f_lammps_last_thermo_num() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  USE liblammps
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_last_thermo_num
  INTEGER(c_int), POINTER :: ival

  ival = lmp%last_thermo('num',1)
  f_lammps_last_thermo_num = ival
END FUNCTION f_lammps_last_thermo_num

FUNCTION f_lammps_last_thermo_type(idx) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  USE liblammps
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), VALUE :: idx
  INTEGER(c_int) :: f_lammps_last_thermo_type
  INTEGER(c_int), POINTER :: ival

  ival = lmp%last_thermo('type',idx)
  f_lammps_last_thermo_type = ival
END FUNCTION f_lammps_last_thermo_type

FUNCTION f_lammps_last_thermo_string(idx) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_ptr, c_null_ptr
  USE liblammps
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  INTEGER(c_int), VALUE :: idx
  TYPE(c_ptr) :: f_lammps_last_thermo_string
  CHARACTER(LEN=12) :: buffer

  buffer = lmp%last_thermo('keyword',idx)
  IF (LEN_TRIM(buffer) > 0) THEN
      f_lammps_last_thermo_string = f2c_string(buffer)
  ELSE
      f_lammps_last_thermo_string = c_null_ptr
  END IF
END FUNCTION f_lammps_last_thermo_string

FUNCTION f_lammps_last_thermo_int(idx) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_int64_t
  USE liblammps
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), VALUE :: idx
  INTEGER(c_int), POINTER :: ival
  INTEGER(c_int64_t), POINTER :: bval
  INTEGER(c_int) :: f_lammps_last_thermo_int
  INTEGER :: size_bigint

  size_bigint = lmp%extract_setting('bigint')
  IF (size_bigint == 4) THEN
      ival = lmp%last_thermo('data',idx)
      f_lammps_last_thermo_int = INT(ival)
  ELSE
      bval = lmp%last_thermo('data',idx)
      f_lammps_last_thermo_int = INT(bval)
  END IF
END FUNCTION f_lammps_last_thermo_int

FUNCTION f_lammps_last_thermo_double(idx) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  USE liblammps
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), VALUE :: idx
  REAL(c_double), POINTER :: dval
  REAL(c_double) :: f_lammps_last_thermo_double

  dval = lmp%last_thermo('data',idx)
  f_lammps_last_thermo_double = dval
END FUNCTION f_lammps_last_thermo_double

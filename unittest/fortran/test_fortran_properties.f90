SUBROUTINE f_lammps_memory_usage(meminfo) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double), DIMENSION(3), INTENT(OUT) :: meminfo

   CALL lmp%memory_usage(meminfo)
END SUBROUTINE f_lammps_memory_usage

FUNCTION f_lammps_get_mpi_comm() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE liblammps
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int) :: f_lammps_get_mpi_comm

   f_lammps_get_mpi_comm = lmp%get_mpi_comm()
END FUNCTION f_lammps_get_mpi_comm

FUNCTION f_lammps_extract_setting(Cstr) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_char, c_null_char
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int) :: f_lammps_extract_setting
   CHARACTER(KIND=c_char, LEN=1), DIMENSION(*), INTENT(IN) :: Cstr
   INTEGER :: strlen, i
   CHARACTER(LEN=:), ALLOCATABLE :: Fstr

   i = 1
   DO WHILE (Cstr(i) /= c_null_char)
      i = i + 1
   END DO
   strlen = i
   ALLOCATE(CHARACTER(LEN=strlen) :: Fstr)
   DO i = 1, strlen
      Fstr(i:i) = Cstr(i)
   END DO
   f_lammps_extract_setting = lmp%extract_setting(Fstr)
   DEALLOCATE(Fstr)
END FUNCTION f_lammps_extract_setting

FUNCTION f_lammps_has_error() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int) :: f_lammps_has_error

   IF (lmp%has_error()) THEN
      f_lammps_has_error = 1_c_int
   ELSE
      f_lammps_has_error = 0_c_int
   END IF
END FUNCTION f_lammps_has_error

FUNCTION f_lammps_get_last_error_message(errmesg, errlen) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_char, c_ptr, C_F_POINTER, &
      c_null_char
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int) :: f_lammps_get_last_error_message
   CHARACTER(KIND=c_char), DIMENSION(*) :: errmesg
   INTEGER(c_int), VALUE, INTENT(IN) :: errlen
   CHARACTER(LEN=:), ALLOCATABLE :: buffer
   INTEGER :: status, i

   ! copy error message to buffer
   ALLOCATE(CHARACTER(errlen) :: buffer)
   CALL lmp%get_last_error_message(buffer, status)
   f_lammps_get_last_error_message = status
   ! and copy to C style string
   errmesg(1:errlen) = c_null_char
   DO i=1, errlen
      errmesg(i) = buffer(i:i)
      IF (buffer(i:i) == c_null_char) EXIT
   END DO
   errmesg(errlen) = c_null_char
   DEALLOCATE(buffer)
END FUNCTION f_lammps_get_last_error_message

FUNCTION f_lammps_get_image_flags_int(ix, iy, iz) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: ix, iy, iz
  INTEGER(c_int) :: f_lammps_get_image_flags_int

  f_lammps_get_image_flags_int = lmp%encode_image_flags(ix, iy, iz)
END FUNCTION f_lammps_get_image_flags_int

FUNCTION f_lammps_get_image_flags_bigint(ix, iy, iz) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: ix, iy, iz
  INTEGER(c_int64_t) :: f_lammps_get_image_flags_bigint

  f_lammps_get_image_flags_bigint = lmp%encode_image_flags(ix, iy, iz)
END FUNCTION f_lammps_get_image_flags_bigint

SUBROUTINE f_lammps_decode_image_flags(image, flag) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: image
  INTEGER(c_int), INTENT(OUT) :: flag(3)

  CALL lmp%decode_image_flags(image, flag)
END SUBROUTINE f_lammps_decode_image_flags

SUBROUTINE f_lammps_decode_image_flags_bigbig(image, flag) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int64_t), INTENT(IN), VALUE :: image
  INTEGER(c_int), INTENT(OUT) :: flag(3)

  CALL lmp%decode_image_flags(image, flag)
END SUBROUTINE f_lammps_decode_image_flags_bigbig

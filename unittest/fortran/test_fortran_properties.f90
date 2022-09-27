MODULE keepprops
  USE liblammps
  IMPLICIT NONE
  TYPE(LAMMPS) :: lmp
  CHARACTER(len=40), DIMENSION(3), PARAMETER :: demo_input = &
      [ CHARACTER(len=40) ::                                &
      'region       box block 0 $x 0 2 0 2',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(len=40), DIMENSION(2), PARAMETER :: cont_input = &
      [ CHARACTER(len=40) ::                                &
      'create_atoms 1 single &',                            &
      ' 0.2 0.1 0.1' ]
END MODULE keepprops

FUNCTION f_lammps_version () BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int
   USE liblammps
   USE keepprops, ONLY : lmp
   IMPLICIT NONE
   INTEGER (C_int) :: f_lammps_version

   f_lammps_version = lmp%version()
END FUNCTION f_lammps_version

SUBROUTINE f_lammps_memory_usage (meminfo) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double
   USE liblammps
   USE keepprops, ONLY : lmp
   IMPLICIT NONE
   REAL (C_double), DIMENSION(3), INTENT(OUT) :: meminfo

   CALL lmp%memory_usage(meminfo)
END SUBROUTINE f_lammps_memory_usage

FUNCTION f_lammps_get_mpi_comm () BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int
   USE liblammps
   USE keepprops, ONLY : lmp
   IMPLICIT NONE
   INTEGER (C_int) :: f_lammps_get_mpi_comm

   f_lammps_get_mpi_comm = lmp%get_mpi_comm()
END FUNCTION f_lammps_get_mpi_comm

FUNCTION f_lammps_extract_setting (Cstr) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int, C_char
   USE keepprops, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER (C_int) :: f_lammps_extract_setting
   CHARACTER (KIND=C_char, LEN=1), DIMENSION(*), INTENT(IN) :: Cstr
   INTEGER :: strlen, i
   CHARACTER (LEN=:), ALLOCATABLE :: Fstr

   i = 1
   DO WHILE (Cstr(i) /= ACHAR(0))
      i = i + 1
   END DO
   strlen = i
   allocate ( CHARACTER(LEN=strlen) :: Fstr)
   FORALL (i=1:strlen)
      Fstr(i:i) = Cstr(i)
   END FORALL
   f_lammps_extract_setting = lmp%extract_setting(Fstr)
   deallocate (Fstr)
END FUNCTION f_lammps_extract_setting

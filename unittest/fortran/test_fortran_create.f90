MODULE MYMPI
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lmp_comm_split

  INTERFACE
      FUNCTION lmp_comm_split(color, key) BIND(C,name='create_mpi_comm_split')
        IMPORT :: c_int
        IMPLICIT NONE
        INTEGER(c_int), VALUE, INTENT(IN) :: color, key
        INTEGER(c_int)                    :: lmp_comm_split
      END FUNCTION lmp_comm_split
  END INTERFACE
END MODULE MYMPI

FUNCTION f_lammps_no_mpi_no_args() BIND(C, name="f_lammps_no_mpi_no_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepstuff,      ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_no_mpi_no_args

  lmp = lammps()
  f_lammps_no_mpi_no_args = lmp%handle
END FUNCTION f_lammps_no_mpi_no_args

FUNCTION f_lammps_no_mpi_with_args() BIND(C, name="f_lammps_no_mpi_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepstuff,      ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_no_mpi_with_args

  CHARACTER(len=12), DIMENSION(4), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps', '-log', 'none', '-nocite' ]

  lmp = lammps(args)
  f_lammps_no_mpi_with_args = lmp%handle
END FUNCTION f_lammps_no_mpi_with_args

FUNCTION f_lammps_open_no_args() BIND(C, name="f_lammps_open_no_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE MYMPI,           ONLY: lmp_comm_split
  USE liblammps
  USE keepstuff,      ONLY: lmp,mycomm
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_no_args
  INTEGER     :: color, key

  color = 1
  key = 1
  mycomm = lmp_comm_split(color, key)
  lmp = lammps(comm=mycomm)
  f_lammps_open_no_args = lmp%handle
END FUNCTION f_lammps_open_no_args

FUNCTION f_lammps_open_with_args() BIND(C, name="f_lammps_open_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE MYMPI,         ONLY: lmp_comm_split
  USE liblammps
  USE keepstuff,      ONLY: lmp,mycomm
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_with_args
  INTEGER     :: color, key

  CHARACTER(len=12), DIMENSION(4), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps', '-log', 'none', '-nocite' ]

  color = 2
  key = 1
  mycomm = lmp_comm_split(color, key)
  lmp = lammps(args,mycomm)
  f_lammps_open_with_args = lmp%handle
END FUNCTION f_lammps_open_with_args

SUBROUTINE f_lammps_close() BIND(C, name="f_lammps_close")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

FUNCTION f_lammps_get_comm() BIND(C, name="f_lammps_get_comm")
  USE ISO_C_BINDING, ONLY: c_int
  USE liblammps
  USE keepstuff, ONLY: mycomm
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_get_comm

  f_lammps_get_comm = mycomm
END FUNCTION f_lammps_get_comm

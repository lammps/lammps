MODULE keepcreate
  USE liblammps
  TYPE(LAMMPS) :: lmp
  INTEGER :: mycomm
END MODULE keepcreate

FUNCTION f_lammps_no_mpi_no_args() BIND(C, name="f_lammps_no_mpi_no_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepcreate,      ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_no_mpi_no_args

  lmp = lammps()
  f_lammps_no_mpi_no_args = lmp%handle
END FUNCTION f_lammps_no_mpi_no_args

FUNCTION f_lammps_no_mpi_with_args() BIND(C, name="f_lammps_no_mpi_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepcreate,      ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_no_mpi_with_args

  CHARACTER(len=12), DIMENSION(4), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps', '-log', 'none', '-nocite' ]

  lmp = lammps(args)
  f_lammps_no_mpi_with_args = lmp%handle
END FUNCTION f_lammps_no_mpi_with_args

FUNCTION f_lammps_open_no_args() BIND(C, name="f_lammps_open_no_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE MPI,           ONLY: MPI_COMM_WORLD, mpi_comm_split
  USE liblammps
  USE keepcreate,      ONLY: lmp,mycomm
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_no_args
  INTEGER     :: color, key, ierr

  color = 1
  key = 1
  CALL mpi_comm_split(MPI_COMM_WORLD, color, key, mycomm, ierr)
  lmp = lammps(comm=mycomm)
  f_lammps_open_no_args = lmp%handle
END FUNCTION f_lammps_open_no_args

FUNCTION f_lammps_open_with_args() BIND(C, name="f_lammps_open_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE MPI,           ONLY: MPI_COMM_WORLD, mpi_comm_split
  USE liblammps
  USE keepcreate,      ONLY: lmp,mycomm
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_with_args
  INTEGER     :: color, key, ierr

  CHARACTER(len=12), DIMENSION(4), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps', '-log', 'none', '-nocite' ]

  color = 2
  key = 1
  CALL mpi_comm_split(MPI_COMM_WORLD, color, key, mycomm, ierr)
  lmp = lammps(args,mycomm)
  f_lammps_open_with_args = lmp%handle
END FUNCTION f_lammps_open_with_args

SUBROUTINE f_lammps_close() BIND(C, name="f_lammps_close")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepcreate, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

FUNCTION f_lammps_get_comm() BIND(C, name="f_lammps_get_comm")
  USE liblammps
  USE keepcreate, ONLY: mycomm
  IMPLICIT NONE
  INTEGER :: f_lammps_get_comm

  f_lammps_get_comm = mycomm
END FUNCTION f_lammps_get_comm

MODULE keepdata
  USE liblammps
  TYPE(LAMMPS) :: lmp
  INTEGER :: mycomm
END MODULE keepdata


FUNCTION f_lammps_open_no_args() BIND(C, name="f_lammps_open_no_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepdata,      ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_no_args
  
  lmp = lammps()
  f_lammps_open_no_args = lmp%handle
END FUNCTION f_lammps_open_no_args

FUNCTION f_lammps_open_with_args() BIND(C, name="f_lammps_open_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE MPI,           ONLY: MPI_COMM_WORLD, mpi_comm_split
  USE liblammps
  USE keepdata,      ONLY: lmp,mycomm
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_with_args
  INTEGER     :: color, key, ierr
  
  CHARACTER(len=*), DIMENSION(4), PARAMETER :: args = &
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
  USE keepdata, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close
  
FUNCTION f_lammps_get_comm() BIND(C, name="f_lammps_get_comm")
  USE liblammps
  USE keepdata, ONLY: mycomm
  IMPLICIT NONE
  INTEGER :: f_lammps_get_comm

  f_lammps_get_comm = mycomm
END FUNCTION f_lammps_get_comm
  


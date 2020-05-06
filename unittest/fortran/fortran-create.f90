MODULE keepdata
  USE liblammps
  TYPE(LAMMPS) :: lmp
END MODULE keepdata


FUNCTION f_lammps_open_noargs() BIND(C, name="f_lammps_open_noargs")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepdata,      ONLY: lmp
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_open_noargs
  
  CHARACTER(len=*), DIMENSION(1), PARAMETER :: args = &
      [ CHARACTER(len=12) :: 'liblammps' ]

  lmp = lammps()
  f_lammps_open_noargs = lmp%handle
END FUNCTION f_lammps_open_noargs


SUBROUTINE f_lammps_close() BIND(C, name="f_lammps_close")
  USE liblammps
  USE keepdata, ONLY: lmp
  IMPLICIT NONE

  call lmp%close()
END SUBROUTINE f_lammps_close
  


MODULE keepcmds
  USE liblammps
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
END MODULE keepcmds

FUNCTION f_lammps_with_args() BIND(C, name="f_lammps_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepcmds,      ONLY: lmp
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
  USE keepcmds, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

FUNCTION f_lammps_get_natoms() BIND(C, name="f_lammps_get_natoms")
  USE ISO_C_BINDING, ONLY: c_null_ptr, c_double
  USE liblammps
  USE keepcmds, ONLY: lmp
  IMPLICIT NONE
  REAL(c_double) :: f_lammps_get_natoms

  f_lammps_get_natoms = lmp%get_natoms()
END FUNCTION f_lammps_get_natoms

SUBROUTINE f_lammps_file() BIND(C, name="f_lammps_file")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepcmds, ONLY: lmp, demo_input, cont_input
  IMPLICIT NONE
  INTEGER :: i
  CHARACTER(len=*), PARAMETER :: demo_file = 'in.test', cont_file = 'in.cont'

  OPEN(10, file=demo_file, status='replace')
  WRITE(10, fmt='(A)') (demo_input(i),i=1,SIZE(demo_input))
  CLOSE(10)
  OPEN(11, file=cont_file, status='replace')
  WRITE(11, fmt='(A)') (cont_input(i),i=1,SIZE(cont_input))
  CLOSE(11)
  CALL lmp%file(demo_file)
  CALL lmp%file(cont_file)
  OPEN(12, file=demo_file, status='old')
  CLOSE(12, status='delete')
  OPEN(13, file=cont_file, status='old')
  CLOSE(13, status='delete')
END SUBROUTINE f_lammps_file

SUBROUTINE f_lammps_command() BIND(C, name="f_lammps_command")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepcmds, ONLY: lmp, demo_input
  IMPLICIT NONE
  INTEGER :: i

  DO i=1,SIZE(demo_input)
      call lmp%command(demo_input(i))
  END DO
END SUBROUTINE f_lammps_command

SUBROUTINE f_lammps_commands_list() BIND(C, name="f_lammps_commands_list")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepcmds, ONLY: lmp, demo_input, cont_input
  IMPLICIT NONE

  CALL lmp%commands_list(demo_input)
  CALL lmp%commands_list(cont_input)
END SUBROUTINE f_lammps_commands_list

SUBROUTINE f_lammps_commands_string() BIND(C, name="f_lammps_commands_string")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepcmds, ONLY: lmp, demo_input, cont_input
  IMPLICIT NONE
  INTEGER :: i
  CHARACTER(len=512) :: cmds

  cmds = ''
  DO i=1,SIZE(demo_input)
      cmds = TRIM(cmds) // TRIM(demo_input(i)) // NEW_LINE('A')
  END DO
  DO i=1,SIZE(cont_input)
      cmds = TRIM(cmds) // TRIM(cont_input(i)) // NEW_LINE('A')
  END DO

  CALL lmp%commands_string(cmds)
END SUBROUTINE f_lammps_commands_string

FUNCTION f_lammps_with_args() BIND(C)
  USE ISO_C_BINDING, ONLY: c_ptr
  USE LIBLAMMPS
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

SUBROUTINE f_lammps_setup_create_atoms() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, big_input, cont_input, more_input
  IMPLICIT NONE

  !CALL lmp%command('atom_modify map array')
  CALL lmp%commands_list(big_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(more_input)
END SUBROUTINE f_lammps_setup_create_atoms

SUBROUTINE f_lammps_create_three_atoms() BIND(C)
  USE ISO_C_BINDING, ONLY: c_double, c_int, c_int64_t
  USE keepstuff, ONLY : lmp, big_input, cont_input, more_input
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(3) :: new_ids, new_images, new_types
  INTEGER(c_int64_t), DIMENSION(3) :: new_big_ids, new_big_images
  REAL(c_double), DIMENSION(9) :: new_x, new_v
  LOGICAL :: wrap
  INTEGER(c_int) :: tagint_size

  new_ids = [4, 6, 5]
  new_big_ids = [4, 6, 5]
  new_images = [0, 0, 1]
  new_big_images = [0, 0, 1]
  new_types = [1, 1, 1]
  new_x = [ 1.0_c_double, 1.8_c_double, 2.718281828_c_double, &
            0.6_c_double, 0.8_c_double, 2.2_c_double, &
            1.8_c_double, 0.1_c_double, 1.8_c_double ]
  new_v = [ 0.0_c_double, 1.0_c_double, -1.0_c_double, &
            0.1_c_double, 0.2_c_double, -0.2_c_double, &
            1.0_c_double, -1.0_c_double, 3.0_c_double ]
  wrap = .FALSE.
  tagint_size = lmp%extract_setting('tagint')
  IF ( tagint_size == 4_c_int ) THEN
    CALL lmp%create_atoms(new_ids, new_types, new_x, new_v, new_images, wrap)
  ELSE
    CALL lmp%create_atoms(new_big_ids, new_types, new_x, new_v, &
      new_big_images, wrap)
  END IF
END SUBROUTINE f_lammps_create_three_atoms

! vim: ts=2 sts=2 sw=2 et

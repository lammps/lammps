FUNCTION f_lammps_with_args() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr
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
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_null_ptr
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
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int, c_int64_t
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(3) :: new_ids, new_images, new_types
  INTEGER(c_int64_t), DIMENSION(3) :: new_big_ids, new_big_images
  REAL(c_double), DIMENSION(9) :: new_x, new_v
  LOGICAL :: wrap
  INTEGER(c_int) :: tagint_size

  new_ids = [4, 6, 5]
  new_big_ids = [4, 6, 5]
  tagint_size = lmp%extract_setting('tagint')
  IF ( tagint_size == 4_c_int ) THEN
    new_images(1) = lmp%encode_image_flags(1, -1, 3)
    new_images(2) = lmp%encode_image_flags(-2, 0, 0)
    new_images(3) = lmp%encode_image_flags(-2, -2, 1)
  ELSE
    new_big_images(1) = lmp%encode_image_flags(1, -1, 3)
    new_big_images(2) = lmp%encode_image_flags(-2, 0, 0)
    new_big_images(3) = lmp%encode_image_flags(-2, -2, 1)
  END IF
  new_types = [1, 1, 1]
  new_x = [ 1.0_c_double, 1.8_c_double, 2.718281828_c_double, &
            0.6_c_double, 0.8_c_double, 2.2_c_double, &
            1.8_c_double, 0.1_c_double, 1.8_c_double ]
  new_v = [ 0.0_c_double, 1.0_c_double, -1.0_c_double, &
            0.1_c_double, 0.2_c_double, -0.2_c_double, &
            1.0_c_double, -1.0_c_double, 3.0_c_double ]
  wrap = .FALSE.
  IF ( tagint_size == 4_c_int ) THEN
    CALL lmp%create_atoms(new_ids, new_types, new_x, new_v, new_images, wrap)
  ELSE
    CALL lmp%create_atoms(new_big_ids, new_types, new_x, new_v, &
      new_big_images, wrap)
  END IF
END SUBROUTINE f_lammps_create_three_atoms

SUBROUTINE f_lammps_create_two_more() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(2) :: new_types
  REAL(c_double), DIMENSION(6) :: new_x

  new_types = [1_c_int, 1_c_int]
  new_x = [0.1_c_double, 1.9_c_double, 3.8_c_double, &
           1.2_c_double, 2.1_c_double, 1.25_c_double]
  CALL lmp%create_atoms(type=new_types, x=new_x)
END SUBROUTINE f_lammps_create_two_more

SUBROUTINE f_lammps_create_two_more_small() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(2) :: new_types
  REAL(c_double), DIMENSION(6) :: new_x
  INTEGER(c_int), DIMENSION(2) :: new_id, new_image

  new_types = [1_c_int, 1_c_int]
  new_x = [0.1_c_double, 1.9_c_double, 3.8_c_double, &
           1.2_c_double, 2.1_c_double, 1.25_c_double]
  new_id = [8_c_int, 7_c_int]
  new_image(1) = lmp%encode_image_flags(1,0,0)
  new_image(2) = lmp%encode_image_flags(-1,0,0)
  CALL lmp%create_atoms(id=new_id, image=new_image, type=new_types, x=new_x)
END SUBROUTINE f_lammps_create_two_more_small

SUBROUTINE f_lammps_create_two_more_big() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int, c_int64_t
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(2) :: new_types
  REAL(c_double), DIMENSION(6) :: new_x
  INTEGER(c_int64_t), DIMENSION(2) :: new_id, new_image

  new_types = [1_c_int, 1_c_int]
  new_x = [0.1_c_double, 1.9_c_double, 3.8_c_double, &
           1.2_c_double, 2.1_c_double, 1.25_c_double]
  new_id = [8_c_int64_t, 7_c_int64_t]
  new_image(1) = lmp%encode_image_flags(1,0,0)
  new_image(2) = lmp%encode_image_flags(-1,0,0)
  CALL lmp%create_atoms(id=new_id, image=new_image, type=new_types, x=new_x)
END SUBROUTINE f_lammps_create_two_more_big

SUBROUTINE f_lammps_create_two_more_small2() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(2) :: new_types
  REAL(c_double), DIMENSION(6) :: new_x
  INTEGER(c_int), DIMENSION(2) :: new_id

  new_types = [1_c_int, 1_c_int]
  new_x = [0.1_c_double, 1.9_c_double, 3.8_c_double, &
           1.2_c_double, 2.1_c_double, 1.25_c_double]
  new_id = [8_c_int, 7_c_int]
  CALL lmp%create_atoms(id=new_id, type=new_types, x=new_x)
END SUBROUTINE f_lammps_create_two_more_small2

SUBROUTINE f_lammps_create_two_more_big2() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_double, c_int, c_int64_t
  USE keepstuff, ONLY : lmp
  USE LIBLAMMPS
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(2) :: new_types
  REAL(c_double), DIMENSION(6) :: new_x
  INTEGER(c_int64_t), DIMENSION(2) :: new_id

  new_types = [1_c_int, 1_c_int]
  new_x = [0.1_c_double, 1.9_c_double, 3.8_c_double, &
           1.2_c_double, 2.1_c_double, 1.25_c_double]
  new_id = [8_c_int64_t, 7_c_int64_t]
  CALL lmp%create_atoms(id=new_id, type=new_types, x=new_x)
END SUBROUTINE f_lammps_create_two_more_big2

! vim: ts=2 sts=2 sw=2 et

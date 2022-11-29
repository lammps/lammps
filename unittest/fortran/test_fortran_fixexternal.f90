MODULE ext_stuff
  USE, INTRINSIC :: ISO_Fortran_ENV, ONLY : error_unit
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int, c_int64_t, c_loc
  USE LIBLAMMPS
  IMPLICIT NONE

  REAL(c_double), SAVE :: direction = 1.0_c_double

CONTAINS

  SUBROUTINE f_lammps_reverse_direction() BIND(C)
    direction = -direction
  END SUBROUTINE f_lammps_reverse_direction

  SUBROUTINE f_callback_ss(instance, timestep, id, x, f)
    CLASS(*), INTENT(INOUT) :: instance
    INTEGER(c_int) :: timestep
    INTEGER(c_int), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f

    WHERE (id == 1)
      f(1,:) = 1.0_c_double
      f(2,:) = -1.0_c_double
      f(3,:) = 1.25_c_double
    ELSEWHERE
      f(1,:) = -1.0_c_double
      f(2,:) = +1.0_c_double
      f(3,:) = -1.25_c_double
    END WHERE
  END SUBROUTINE f_callback_ss

  SUBROUTINE f_callback_sb(instance, timestep, id, x, f)
    CLASS(*), INTENT(INOUT) :: instance
    INTEGER(c_int64_t) :: timestep
    INTEGER(c_int), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f

    WHERE (id == 1_c_int)
      f(1,:) = 1.0_c_double
      f(2,:) = -1.0_c_double
      f(3,:) = 1.25_c_double
    ELSEWHERE
      f(1,:) = -1.0_c_double
      f(2,:) = +1.0_c_double
      f(3,:) = -1.25_c_double
    END WHERE
  END SUBROUTINE f_callback_sb

  SUBROUTINE f_callback_bb(instance, timestep, id, x, f)
    CLASS(*), INTENT(INOUT) :: instance
    INTEGER(c_int64_t) :: timestep
    INTEGER(c_int64_t), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f

    WHERE (id == 1_c_int64_t)
      f(1,:) = 1.0_c_double
      f(2,:) = -1.0_c_double
      f(3,:) = 1.25_c_double
    ELSEWHERE
      f(1,:) = -1.0_c_double
      f(2,:) = +1.0_c_double
      f(3,:) = -1.25_c_double
    END WHERE
  END SUBROUTINE f_callback_bb

  SUBROUTINE f_callback2_ss(entity, timestep, id, x, f)
    CLASS(*), INTENT(INOUT), target :: entity
    INTEGER(c_int) :: timestep
    INTEGER(c_int), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f

    SELECT TYPE (entity)
      TYPE IS (REAL(c_double))
        WHERE (id == 1_c_int)
          f(1,:) = SIGN(1.0_c_double, entity) * 2.0_c_double
          f(2,:) = SIGN(1.0_c_double, entity) * (-2.0_c_double)
          f(3,:) = SIGN(1.0_c_double, entity) * 2.5_c_double
        ELSEWHERE
          f(1,:) = SIGN(1.0_c_double, entity) * (-2.0_c_double)
          f(2,:) = SIGN(1.0_c_double, entity) * 2.0_c_double
          f(3,:) = SIGN(1.0_c_double, entity) * (-2.5_c_double)
        END WHERE
      CLASS DEFAULT
        WRITE(error_unit,'(A)') 'ERROR: Failed to resolve "entity" in&
          & f_callback2_ss'
        STOP 1
    END SELECT
  END SUBROUTINE f_callback2_ss

  SUBROUTINE f_callback2_sb(entity, timestep, id, x, f)
    CLASS(*), INTENT(INOUT), target :: entity
    INTEGER(c_int64_t) :: timestep
    INTEGER(c_int), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f

    SELECT TYPE (entity)
      TYPE IS (REAL(c_double))
        WHERE (id == 1_c_int)
          f(1,:) = SIGN(1.0_c_double, entity) * 2.0_c_double
          f(2,:) = SIGN(1.0_c_double, entity) * (-2.0_c_double)
          f(3,:) = SIGN(1.0_c_double, entity) * 2.5_c_double
        ELSEWHERE
          f(1,:) = SIGN(1.0_c_double, entity) * (-2.0_c_double)
          f(2,:) = SIGN(1.0_c_double, entity) * 2.0_c_double
          f(3,:) = SIGN(1.0_c_double, entity) * (-2.5_c_double)
        END WHERE
      CLASS DEFAULT
        WRITE(error_unit,'(A)') 'ERROR: Failed to resolve "entity" in&
          & f_callback2_sb'
        STOP 1
    END SELECT
  END SUBROUTINE f_callback2_sb

  SUBROUTINE f_callback2_bb(entity, timestep, id, x, f)
    CLASS(*), INTENT(INOUT), target :: entity
    INTEGER(c_int64_t) :: timestep
    INTEGER(c_int64_t), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f

    SELECT TYPE (entity)
      TYPE IS (REAL(c_double))
        WHERE (id == 1_c_int64_t)
          f(1,:) = SIGN(1.0_c_double, entity) * 2.0_c_double
          f(2,:) = SIGN(1.0_c_double, entity) * (-2.0_c_double)
          f(3,:) = SIGN(1.0_c_double, entity) * 2.5_c_double
        ELSEWHERE
          f(1,:) = SIGN(1.0_c_double, entity) * (-2.0_c_double)
          f(2,:) = SIGN(1.0_c_double, entity) * 2.0_c_double
          f(3,:) = SIGN(1.0_c_double, entity) * (-2.5_c_double)
        END WHERE
      CLASS DEFAULT
        WRITE(error_unit,'(A)') 'ERROR: Failed to resolve "entity" in&
          & f_callback2_sb'
        STOP 1
    END SELECT
  END SUBROUTINE f_callback2_bb
END MODULE ext_stuff

FUNCTION f_lammps_with_args() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr
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

SUBROUTINE f_lammps_setup_fix_external() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, demo_input, cont_input, pair_input
  IMPLICIT NONE

  CALL lmp%commands_list(demo_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('neigh_modify exclude group all all')
  CALL lmp%command('fix ext1 all external pf/callback 1 1')
  CALL lmp%command('fix ext2 all external pf/callback 1 1')
END SUBROUTINE f_lammps_setup_fix_external

SUBROUTINE f_lammps_set_fix_external_callbacks() BIND(C)
  USE ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  USE ext_stuff
  IMPLICIT NONE
  INTEGER :: size_bigint, size_tagint, nlocal

  nlocal = lmp%extract_setting('nlocal')

  size_bigint = lmp%extract_setting('bigint')
  size_tagint = lmp%extract_setting('tagint')
  IF (size_bigint == 4_c_int .AND. size_tagint == 4_c_int) THEN
    CALL lmp%set_fix_external_callback('ext1', f_callback_ss)
    CALL lmp%set_fix_external_callback('ext2', f_callback2_ss, direction)
  ELSE IF (size_bigint == 8_c_int .AND. size_tagint == 8_c_int) THEN
    CALL lmp%set_fix_external_callback('ext1', f_callback_bb)
    CALL lmp%set_fix_external_callback('ext2', f_callback2_bb, direction)
  ELSE
    CALL lmp%set_fix_external_callback('ext1', f_callback_sb)
    CALL lmp%set_fix_external_callback('ext2', f_callback2_sb, direction)
  END IF
END SUBROUTINE f_lammps_set_fix_external_callbacks

SUBROUTINE f_lammps_get_force (i, ptr) BIND(C)
  USE ISO_C_BINDING, ONLY : c_int, c_double, c_ptr, C_F_POINTER
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  TYPE(c_ptr), INTENT(IN), VALUE :: ptr
  REAL(c_double), DIMENSION(:,:), POINTER :: force => NULL()
  REAL(c_double), DIMENSION(:), POINTER :: f => NULL()

  CALL C_F_POINTER(ptr, f, [3])
  force = lmp%extract_atom('f')
  f = force(:,i)
END SUBROUTINE f_lammps_get_force

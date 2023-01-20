MODULE ext_stuff
  USE, INTRINSIC :: ISO_Fortran_ENV, ONLY : error_unit
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int, c_int64_t, c_loc
  USE LIBLAMMPS
  IMPLICIT NONE

  INTEGER, PARAMETER :: vec_length = 8
  REAL(c_double), SAVE :: direction = 1.0_c_double
  REAL(c_double), DIMENSION(:,:), POINTER, SAVE :: f3 => NULL(), f4 => NULL()

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
    REAL(c_double), DIMENSION(SIZE(id)) :: e
    REAL(c_double), DIMENSION(6,SIZE(id)) :: v

    WHERE (id == 1)
      f(1,:) = 1.0_c_double
      f(2,:) = -1.0_c_double
      f(3,:) = 1.25_c_double
      e = 1.0_c_double
      v(1,:) = 1.0_c_double
      v(2,:) = 2.0_c_double
      v(3,:) = -1.0_c_double
      v(4,:) = -2.0_c_double
      v(5,:) = 3.0_c_double
      v(6,:) = -3.0_c_double
    ELSEWHERE
      f(1,:) = -1.0_c_double
      f(2,:) = +1.0_c_double
      f(3,:) = -1.25_c_double
      e = 10.0_c_double
      v(1,:) = 10.0_c_double
      v(2,:) = 20.0_c_double
      v(3,:) = -10.0_c_double
      v(4,:) = -20.0_c_double
      v(5,:) = 30.0_c_double
      v(6,:) = -30.0_c_double
    END WHERE
    SELECT TYPE (instance)
      CLASS IS (lammps)
        CALL instance%fix_external_set_energy_peratom('ext1', e)
        CALL instance%fix_external_set_virial_peratom('ext1', v)
      CLASS DEFAULT
        WRITE(error_unit,*) 'UMM...this should never happen.'
        STOP 1
    END SELECT
  END SUBROUTINE f_callback_ss

  SUBROUTINE f_callback_sb(instance, timestep, id, x, f)
    CLASS(*), INTENT(INOUT) :: instance
    INTEGER(c_int64_t) :: timestep
    INTEGER(c_int), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f
    REAL(c_double), DIMENSION(SIZE(id)) :: e
    REAL(c_double), DIMENSION(6,SIZE(id)) :: v

    WHERE (id == 1_c_int)
      f(1,:) = 1.0_c_double
      f(2,:) = -1.0_c_double
      f(3,:) = 1.25_c_double
      e = 1.0_c_double
      v(1,:) = 1.0_c_double
      v(2,:) = 2.0_c_double
      v(3,:) = -1.0_c_double
      v(4,:) = -2.0_c_double
      v(5,:) = 3.0_c_double
      v(6,:) = -3.0_c_double
    ELSEWHERE
      f(1,:) = -1.0_c_double
      f(2,:) = +1.0_c_double
      f(3,:) = -1.25_c_double
      e = 10.0_c_double
      v(1,:) = 10.0_c_double
      v(2,:) = 20.0_c_double
      v(3,:) = -10.0_c_double
      v(4,:) = -20.0_c_double
      v(5,:) = 30.0_c_double
      v(6,:) = -30.0_c_double
    END WHERE
    SELECT TYPE (instance)
      CLASS IS (lammps)
        CALL instance%fix_external_set_energy_peratom('ext1', e)
        CALL instance%fix_external_set_virial_peratom('ext1', v)
      CLASS DEFAULT
        WRITE(error_unit,*) 'UMM...this should never happen.'
        STOP 1
    END SELECT
  END SUBROUTINE f_callback_sb

  SUBROUTINE f_callback_bb(instance, timestep, id, x, f)
    CLASS(*), INTENT(INOUT) :: instance
    INTEGER(c_int64_t) :: timestep
    INTEGER(c_int64_t), DIMENSION(:), INTENT(IN) :: id
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: x
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: f
    REAL(c_double), DIMENSION(SIZE(id)) :: e
    REAL(c_double), DIMENSION(6,SIZE(id)) :: v

    WHERE (id == 1_c_int64_t)
      f(1,:) = 1.0_c_double
      f(2,:) = -1.0_c_double
      f(3,:) = 1.25_c_double
      e = 1.0_c_double
      v(1,:) = 1.0_c_double
      v(2,:) = 2.0_c_double
      v(3,:) = -1.0_c_double
      v(4,:) = -2.0_c_double
      v(5,:) = 3.0_c_double
      v(6,:) = -3.0_c_double
    ELSEWHERE
      f(1,:) = -1.0_c_double
      f(2,:) = +1.0_c_double
      f(3,:) = -1.25_c_double
      e = 10.0_c_double
      v(1,:) = 10.0_c_double
      v(2,:) = 20.0_c_double
      v(3,:) = -10.0_c_double
      v(4,:) = -20.0_c_double
      v(5,:) = 30.0_c_double
      v(6,:) = -30.0_c_double
    END WHERE
    SELECT TYPE (instance)
      CLASS IS (lammps)
        CALL instance%fix_external_set_energy_peratom('ext1', e)
        CALL instance%fix_external_set_virial_peratom('ext1', v)
      CLASS DEFAULT
        WRITE(error_unit,*) 'UMM...this should never happen.'
        STOP 1
    END SELECT
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
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_fix_external_callback() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, demo_input, cont_input, pair_input
  USE ext_stuff, ONLY : vec_length
  IMPLICIT NONE

  CALL lmp%commands_list(demo_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('neigh_modify exclude group all all')
  CALL lmp%command('fix ext1 all external pf/callback 1 1')
  CALL lmp%command('fix ext2 all external pf/callback 1 1')
  CALL lmp%fix_external_set_vector_length('ext2', vec_length)
END SUBROUTINE f_lammps_setup_fix_external_callback

SUBROUTINE f_lammps_setup_fix_external_array() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, demo_input, cont_input, pair_input
  USE ext_stuff, ONLY : f3, f4
  IMPLICIT NONE

  CALL lmp%commands_list(demo_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('neigh_modify exclude group all all')
  CALL lmp%command('fix ext3 all external pf/array 1')
  CALL lmp%command('fix ext4 all external pf/array 1')
  CALL lmp%command('thermo_style custom step pxx pe etotal')
  CALL lmp%command('thermo_modify norm no')
  CALL lmp%command('thermo 100')
  f3 = lmp%fix_external_get_force('ext3')
  f4 = lmp%fix_external_get_force('ext4')
END SUBROUTINE f_lammps_setup_fix_external_array

SUBROUTINE f_lammps_set_fix_external_callbacks() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  USE ext_stuff
  IMPLICIT NONE
  INTEGER :: size_bigint, size_tagint, nlocal

  nlocal = lmp%extract_setting('nlocal')

  size_bigint = lmp%extract_setting('bigint')
  size_tagint = lmp%extract_setting('tagint')
  IF (size_bigint == 4_c_int .AND. size_tagint == 4_c_int) THEN
    CALL lmp%set_fix_external_callback('ext1', f_callback_ss, lmp)
    CALL lmp%set_fix_external_callback('ext2', f_callback2_ss, direction)
  ELSE IF (size_bigint == 8_c_int .AND. size_tagint == 8_c_int) THEN
    CALL lmp%set_fix_external_callback('ext1', f_callback_bb, lmp)
    CALL lmp%set_fix_external_callback('ext2', f_callback2_bb, direction)
  ELSE
    CALL lmp%set_fix_external_callback('ext1', f_callback_sb, lmp)
    CALL lmp%set_fix_external_callback('ext2', f_callback2_sb, direction)
  END IF
END SUBROUTINE f_lammps_set_fix_external_callbacks

SUBROUTINE f_lammps_get_force (i, ptr) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double, c_ptr, C_F_POINTER
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

SUBROUTINE f_lammps_find_forces() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  USE ext_stuff, ONLY : f3, f4
  IMPLICIT NONE
  INTEGER(c_int) :: size_tagint
  INTEGER(c_int), DIMENSION(:), POINTER :: id
  INTEGER(c_int64_t), DIMENSION(:), POINTER :: tag

  f3(:,:) = 0.0_c_double
  f4(:,:) = 0.0_c_double
  size_tagint = lmp%extract_setting('tagint')
  IF (size_tagint == 4_c_int) THEN
    id = lmp%extract_atom('id')
    WHERE (id == 1_c_int)
      f3(1,:) = 4.0_c_double
      f3(2,:) = -4.0_c_double
      f3(3,:) = 6.0_c_double
      f4(1,:) = 10.0_c_double
      f4(2,:) = -10.0_c_double
      f4(3,:) = 12.0_c_double
    ELSEWHERE
      f3(1,:) = 5.0_c_double
      f3(2,:) = -5.0_c_double
      f3(3,:) = 7.0_c_double
      f4(1,:) = 11.0_c_double
      f4(2,:) = -11.0_c_double
      f4(3,:) = 13.0_c_double
    END WHERE
  ELSE
    tag = lmp%extract_atom('id')
    WHERE (tag == 1_c_int64_t)
      f3(1,:) = 4.0_c_double
      f3(2,:) = -4.0_c_double
      f3(3,:) = 6.0_c_double
      f4(1,:) = 10.0_c_double
      f4(2,:) = -10.0_c_double
      f4(3,:) = 12.0_c_double
    ELSEWHERE
      f3(1,:) = 5.0_c_double
      f3(2,:) = -5.0_c_double
      f3(3,:) = 7.0_c_double
      f4(1,:) = 11.0_c_double
      f4(2,:) = -11.0_c_double
      f4(3,:) = 13.0_c_double
    END WHERE
  END IF
END SUBROUTINE f_lammps_find_forces

SUBROUTINE f_lammps_add_energy() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE

  CALL lmp%fix_external_set_energy_global('ext3', -20.2_c_double);
END SUBROUTINE f_lammps_add_energy

SUBROUTINE f_lammps_set_virial() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE

  CALL lmp%fix_external_set_virial_global('ext4', [1.0_c_double, &
    2.0_c_double, 2.5_c_double, -1.0_c_double, -2.25_c_double, -3.02_c_double])
END SUBROUTINE f_lammps_set_virial

FUNCTION f_lammps_find_peratom_energy(i) RESULT(energy) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: energy
  REAL(c_double), DIMENSION(:), POINTER :: e

  e = lmp%extract_compute('peratom', lmp%style%atom, lmp%type%vector)
  energy = e(i)
END FUNCTION f_lammps_Find_peratom_energy

SUBROUTINE f_lammps_find_peratom_virial(v, i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  REAL(c_double), DIMENSION(6) :: v
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double), DIMENSION(:,:), POINTER :: virial

  virial = lmp%extract_compute('vperatom', lmp%style%atom, lmp%type%array)
  v = virial(:,i)
END SUBROUTINE f_lammps_find_peratom_virial

SUBROUTINE f_lammps_fixexternal_set_vector() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  USE ext_stuff, ONLY : vec_length
  IMPLICIT NONE
  REAL(c_double), DIMENSION(vec_length) :: v
  INTEGER :: i
  DO i = 1, vec_length
    v(i) = REAL(i, c_double)
    CALL lmp%fix_external_set_vector('ext2', i, v(i))
  END DO
END SUBROUTINE f_lammps_fixexternal_set_vector

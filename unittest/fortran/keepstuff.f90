MODULE keepstuff
  USE LIBLAMMPS
  IMPLICIT NONE
  TYPE(LAMMPS), SAVE :: lmp
  INTEGER, SAVE :: mycomm
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: demo_input = &
      [ CHARACTER(LEN=40) ::                                 &
      'region       box block 0 $x 0 2 0 2',                 &
      'create_box 1 box',                                    &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: big_input = &
      [ CHARACTER(LEN=40) ::                                &
      'region       box block 0 $x 0 3 0 4',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(LEN=40), DIMENSION(2), PARAMETER :: cont_input = &
      [ CHARACTER(LEN=40) ::                                 &
      'create_atoms 1 single &',                             &
      ' 0.2 0.1 0.1' ]
  CHARACTER(LEN=60), DIMENSION(18), PARAMETER :: prop_input =         &
      [ CHARACTER(LEN=60) ::  'fix 1 all nve', 'mass 1 3.0',          &
      'fix 2 all property/atom mol q rmass ghost yes',                &
      'fix 3 all property/atom i_one i2_two 2 d_three d2_four 2',     &
      'set group all mass 2.0', 'set atom 1 charge -1',               &
      'set atom 2 charge  1', 'set atom 1 mol 2', 'set atom 2 mol 1', &
      'set atom 1 i_one -3', 'set atom 2 i_one  3',                   &
      'set atom 1 d_three -1.3', 'set atom 2 d_three  3.5',           &
      'set atom 1 i_two[1] -3', 'set atom 2 i_two[2]  3',             &
      'set atom * d_four[1] -1.3', 'set atom * d_four[2]  3.5',       &
      'run 0 post no' ]
  CHARACTER(LEN=40), DIMENSION(1), PARAMETER :: more_input = &
      [ CHARACTER(LEN=40) :: 'create_atoms 1 single 0.5 0.5 0.5' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: pair_input = &
      [ CHARACTER(LEN=40) ::                                 &
      'pair_style lj/cut 2.5',                               &
      'pair_coeff 1 1 1.0 1.0',                              &
      'mass 1 2.0' ]

  INTERFACE
    FUNCTION c_strlen(str) BIND(C,name='strlen')
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr, c_size_t
      IMPLICIT NONE
      TYPE(c_ptr), INTENT(IN), VALUE :: str
      INTEGER(c_size_t) :: c_strlen
    END FUNCTION c_strlen
  END INTERFACE

CONTAINS

  FUNCTION f2c_string(f_string) RESULT(ptr)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_char, c_ptr, c_size_t, &
      c_null_char, C_F_POINTER
    CHARACTER(LEN=*), INTENT(IN)           :: f_string
    CHARACTER(LEN=1, KIND=c_char), POINTER :: c_string(:)
    TYPE(c_ptr) :: ptr
    INTEGER(c_size_t) :: i, n

    INTERFACE
      FUNCTION lammps_malloc(size) BIND(C, name='malloc')
        USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr, c_size_t
        IMPLICIT NONE
        INTEGER(c_size_t), VALUE :: size
        TYPE(c_ptr) :: lammps_malloc
      END FUNCTION lammps_malloc
    END INTERFACE

    n = LEN_TRIM(f_string)
    ptr = lammps_malloc(n+1)
    CALL C_F_POINTER(ptr, c_string, [1])
    DO i=1, n
        c_string(i) = f_string(i:i)
    END DO
    c_string(n+1) = c_null_char
  END FUNCTION f2c_string

END MODULE keepstuff

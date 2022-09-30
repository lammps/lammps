MODULE keepvar
  USE liblammps
  IMPLICIT NONE
  TYPE(LAMMPS) :: lmp
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: demo_input = &
      [ CHARACTER(LEN=40) ::                                &
      'region       box block 0 $x 0 3 0 4',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: cont_input = &
      [ CHARACTER(LEN=40) ::                                &
      'create_atoms 1 single &',                            &
      ' 0.2 0.1 0.1',                                       &
      'create_atoms 1 single 0.5 0.5 0.5' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: pair_input = &
      [ CHARACTER(LEN=40) ::                                &
      'pair_style lj/cut 2.5',                              &
      'pair_coeff 1 1 1.0 1.0',                             &
      'mass 1 2.0' ]
  CHARACTER(LEN=60), DIMENSION(4), PARAMETER :: py_input = &
      [ CHARACTER(LEN=60) :: &
      'python square_it input 1 v_lp return v_square here """', &
      'def square_it(N) :', &
      '  return N*N', &
      '"""' ]

CONTAINS

  FUNCTION absolute_path(filename)
    CHARACTER(LEN=:), ALLOCATABLE :: absolute_path
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=:), ALLOCATABLE :: test_input_directory

print *, 'GOT HERE! filename is ', filename
    test_input_directory = lmp%extract_variable('input_dir')
print *, '          test_input_directory is ', test_input_directory
    absolute_path = test_input_directory // '/' // TRIM(filename)
  END FUNCTION absolute_path

END MODULE keepvar

FUNCTION f_lammps_with_C_args(argc, argv) BIND(C)
  USE ISO_C_BINDING, ONLY: c_ptr, c_char, c_int, c_size_t, c_f_pointer
  USE liblammps
  USE keepvar, ONLY: lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: argc
  TYPE(c_ptr), VALUE :: argv
  TYPE(c_ptr), DIMENSION(:), POINTER :: Fargv
  INTEGER, PARAMETER :: ARG_LENGTH = 80
  TYPE(c_ptr) :: f_lammps_with_C_args
  CHARACTER(LEN=ARG_LENGTH), DIMENSION(argc) :: args
  CHARACTER(LEN=1,KIND=c_char), DIMENSION(:), POINTER :: Cstr
  INTEGER :: i, length, j

  INTERFACE
    FUNCTION c_strlen (str) BIND(C,name='strlen')
      IMPORT :: c_ptr, c_size_t
      IMPLICIT NONE
      TYPE(c_ptr), INTENT(IN), VALUE :: str
      INTEGER(c_size_t) :: c_strlen
    END FUNCTION c_strlen
  END INTERFACE    

  CALL C_F_POINTER(argv, Fargv, [argc])
  DO i = 1, argc
    args(i) = ''
    length = c_strlen(Fargv(i))
    CALL C_F_POINTER(Fargv(i), Cstr, [length])
    FORALL (j = 1:length)
      args(i)(j:j) = Cstr(j)
    END FORALL
  END DO
  
  lmp = lammps(args)
  f_lammps_with_C_args = lmp%handle
END FUNCTION f_lammps_with_C_args

SUBROUTINE f_lammps_close() BIND(C)
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepvar, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_variable () BIND(C)
   USE LIBLAMMPS
   USE keepvar, ONLY : lmp, demo_input, cont_input, pair_input, absolute_path
   IMPLICIT NONE

   CALL lmp%commands_list(demo_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(pair_input)
   CALL lmp%command('variable idx index "hello" "goodbye"')
   CALL lmp%command('variable lp loop 10')
   CALL lmp%command('variable lp_pad loop 10 pad')
   !CALL lmp%command('variable wld world "group1" "group2" "group3"')
   CALL lmp%command('variable wld world "group1"')
   CALL lmp%command('variable uni universe "universe1" "universeA"')
   CALL lmp%command('variable ulp uloop 2')
   CALL lmp%command('variable str index "this is a string"')
   CALL lmp%command('variable fmt format lp %.6G')
   CALL lmp%command('variable fmt_pad format lp %0.6g')
   CALL lmp%command('variable shell getenv SHELL')
!   CALL lmp%command('variable greet file ' // absolute_path('greetings.txt'))
!   CALL lmp%command('variable atfile atomfile ' // absolute_path('atomdata.txt')
   IF ( lmp%config_has_package('PYTHON') ) THEN
      CALL lmp%command('variable py python square_it')
   END IF
   CALL lmp%command('variable time timer')
   CALL lmp%command('variable int internal 4')
   CALL lmp%command("variable nat equal count(all)")
   CALL lmp%command("variable ts equal step")
END SUBROUTINE f_lammps_setup_extract_variable

FUNCTION f_lammps_extract_variable_index_1 () BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int
  USE LIBLAMMPS
  USE keepvar, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_index_1
  CHARACTER(LEN=80) :: str

  str = lmp%extract_variable("idx")
  IF ( trim(str) == 'hello' ) THEN
     f_lammps_extract_variable_index_1 = 1_c_int
  ELSE
     f_lammps_extract_variable_index_1 = 0_c_int
  END IF
END FUNCTION f_lammps_extract_variable_index_1

FUNCTION f_lammps_extract_variable_index_2 () BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int
  USE LIBLAMMPS
  USE keepvar, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_index_2
  CHARACTER(LEN=80) :: str

  str = lmp%extract_variable("idx")
  IF ( trim(str) == 'goodbye' ) THEN
     f_lammps_extract_variable_index_2 = 1_c_int
  ELSE
     f_lammps_extract_variable_index_2 = 0_c_int
  END IF
END FUNCTION f_lammps_extract_variable_index_2

FUNCTION f_lammps_extract_variable_loop () BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_int, C_double
  USE LIBLAMMPS
  USE keepvar, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_loop
  CHARACTER(LEN=80) :: loop

  loop = lmp%extract_variable('lp')
  READ(loop,*) f_lammps_extract_variable_loop
END FUNCTION f_lammps_extract_variable_loop

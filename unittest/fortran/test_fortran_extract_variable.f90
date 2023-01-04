MODULE keepvar
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr, c_size_t, c_char
  USE liblammps
  IMPLICIT NONE

  INTERFACE
    FUNCTION c_path_join(a, b) BIND(C)
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: a, b
      TYPE(c_ptr) :: c_path_join
    END FUNCTION c_path_join

    SUBROUTINE c_free(ptr) BIND(C,name='free')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: ptr
    END SUBROUTINE c_free
  END INTERFACE

CONTAINS

  FUNCTION absolute_path(filename)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr, c_size_t, c_char, C_F_POINTER
    USE keepstuff, ONLY : lmp, f2c_string, c_strlen
    CHARACTER(LEN=:), ALLOCATABLE :: absolute_path
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=256) :: test_input_directory
    TYPE(c_ptr) :: c_test_input_directory, c_absolute_path, c_filename
    CHARACTER(LEN=1,KIND=c_char), DIMENSION(:), POINTER :: F_absolute_path
    INTEGER(c_size_t) :: i, length

    test_input_directory = lmp%extract_variable('input_dir')
    c_test_input_directory = f2c_string(test_input_directory)
    c_filename = f2c_string(filename)
    c_absolute_path = c_path_join(c_test_input_directory, c_filename)
    length = c_strlen(c_absolute_path)
    CALL C_F_POINTER(c_absolute_path, F_absolute_path, [length])
    ALLOCATE(CHARACTER(LEN=length) :: absolute_path)
    DO i = 1, length
      absolute_path(i:i) = F_absolute_path(i)
    END DO
    CALL c_free(c_filename)
    CALL c_free(c_test_input_directory)
    CALL c_free(c_absolute_path)
  END FUNCTION absolute_path

END MODULE keepvar

FUNCTION f_lammps_with_C_args(argc, argv) BIND(C)
  USE ISO_C_BINDING, ONLY: c_ptr, c_char, c_int, c_size_t, C_F_POINTER
  USE liblammps
  USE keepstuff, ONLY: lmp, c_strlen
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: argc
  TYPE(c_ptr), VALUE :: argv
  TYPE(c_ptr), DIMENSION(:), POINTER :: Fargv
  INTEGER, PARAMETER :: ARG_LENGTH = 256
  TYPE(c_ptr) :: f_lammps_with_C_args
  CHARACTER(LEN=ARG_LENGTH), DIMENSION(argc) :: args
  CHARACTER(LEN=1,KIND=c_char), DIMENSION(:), POINTER :: Cstr
  INTEGER(c_size_t):: i, length, j

  CALL C_F_POINTER(argv, Fargv, [argc])
  DO i = 1, argc
    args(i) = ''
    length = c_strlen(Fargv(i))
    CALL C_F_POINTER(Fargv(i), Cstr, [length])
    DO j = 1, length
      args(i)(j:j) = Cstr(j)
    END DO
  END DO

  lmp = lammps(args)
  f_lammps_with_C_args = lmp%handle
END FUNCTION f_lammps_with_C_args

SUBROUTINE f_lammps_close() BIND(C)
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_variable() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, big_input, cont_input, more_input, pair_input
  USE keepvar, ONLY : absolute_path
  IMPLICIT NONE

  ! Had to do this one as one string because lammps_commands_list and
  ! lammps_commands_string do not play well with triple quotes
  CHARACTER(LEN=256), PARAMETER :: py_input = &
      'python square_it input 1 v_lp return v_py format ff here """' &
        // NEW_LINE(' ') // 'def square_it(N) :' &
        // NEW_LINE(' ') // '  return N*N' &
        // NEW_LINE(' ') // '"""'

  CALL lmp%command('atom_modify map array')
  CALL lmp%commands_list(big_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(more_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('variable idx index "hello" "goodbye"')
  CALL lmp%command('variable lp loop 10')
  CALL lmp%command('variable lp_pad loop 10 pad')
  CALL lmp%command('variable wld world "group1"')
  CALL lmp%command('variable uni universe "universe1" "universeA"')
  CALL lmp%command('variable ulp uloop 2')
  CALL lmp%command('variable str string "this is a string"')
  CALL lmp%command('variable ex equal exp(v_lp)')
  CALL lmp%command('variable fmt format ex %.6G')
  CALL lmp%command('variable fmt_pad format ex %08.6g')
  CALL lmp%command('variable username getenv FORTRAN_USER')
  CALL lmp%command('variable greeting file ' // absolute_path('greetings.txt'))
  CALL lmp%command('variable atfile atomfile ' &
    // absolute_path('atomdata.txt'))
  IF (lmp%config_has_package('PYTHON')) THEN
    CALL lmp%command(py_input)
    CALL lmp%command('variable py python square_it')
  END IF
  CALL lmp%command('variable time timer')
  CALL lmp%command('variable int internal 4')
  CALL lmp%command('variable at_z atom z')
  CALL lmp%command("compute COM all com") ! defines a global vector
  CALL lmp%command("variable center vector c_COM")
  ! make sure COM is computable...
  CALL lmp%command("thermo_style custom step pe c_COM[1] v_center[1]")
  CALL lmp%command("run 0") ! so c_COM and v_center have values
END SUBROUTINE f_lammps_setup_extract_variable

FUNCTION f_lammps_extract_variable_index_1() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_index_1
  CHARACTER(LEN=256) :: str

  str = lmp%extract_variable("idx")
  IF (trim(str) == 'hello') THEN
     f_lammps_extract_variable_index_1 = 1_c_int
  ELSE
     f_lammps_extract_variable_index_1 = 0_c_int
  END IF
END FUNCTION f_lammps_extract_variable_index_1

FUNCTION f_lammps_extract_variable_index_2() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_index_2
  CHARACTER(LEN=256) :: str

  str = lmp%extract_variable("idx")
  IF (trim(str) == 'goodbye') THEN
     f_lammps_extract_variable_index_2 = 1_c_int
  ELSE
     f_lammps_extract_variable_index_2 = 0_c_int
  END IF
END FUNCTION f_lammps_extract_variable_index_2

FUNCTION f_lammps_extract_variable_loop() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_loop
  CHARACTER(LEN=256) :: loop

  loop = lmp%extract_variable('lp')
  READ(loop,*) f_lammps_extract_variable_loop
END FUNCTION f_lammps_extract_variable_loop

FUNCTION f_lammps_extract_variable_loop_pad() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_loop_pad
  CHARACTER(LEN=20) :: loop

  loop = lmp%extract_variable('lp_pad')
  f_lammps_extract_variable_loop_pad = f2c_string(loop)
END FUNCTION f_lammps_extract_variable_loop_pad

FUNCTION f_lammps_extract_variable_world() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_world
  CHARACTER(LEN=20) :: world

  world = lmp%extract_variable('wld')
  f_lammps_extract_variable_world = f2c_string(world)
END FUNCTION f_lammps_extract_variable_world

FUNCTION f_lammps_extract_variable_universe() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_universe
  CHARACTER(LEN=20) :: universe

  universe = lmp%extract_variable('uni')
  f_lammps_extract_variable_universe = f2c_string(universe)
END FUNCTION f_lammps_extract_variable_universe

FUNCTION f_lammps_extract_variable_uloop() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: f_lammps_extract_variable_uloop
  CHARACTER(LEN=256) :: uloop

  uloop = lmp%extract_variable('ulp')
  READ(uloop,*) f_lammps_extract_variable_uloop
END FUNCTION f_lammps_extract_variable_uloop

FUNCTION f_lammps_extract_variable_string() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_string
  CHARACTER(LEN=256) :: string

  string = lmp%extract_variable('str')
  f_lammps_extract_variable_string = f2c_string(string)
END FUNCTION f_lammps_extract_variable_string

FUNCTION f_lammps_extract_variable_format() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_format
  CHARACTER(LEN=20) :: form

  form = lmp%extract_variable('fmt')
  f_lammps_extract_variable_format = f2c_string(form)
END FUNCTION f_lammps_extract_variable_format

FUNCTION f_lammps_extract_variable_format_pad() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_format_pad
  CHARACTER(LEN=20) :: form

  form = lmp%extract_variable('fmt_pad')
  f_lammps_extract_variable_format_pad = f2c_string(form)
END FUNCTION f_lammps_extract_variable_format_pad

FUNCTION f_lammps_extract_variable_getenv() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_getenv
  CHARACTER(LEN=40) :: string

  string = lmp%extract_variable('username')
  f_lammps_extract_variable_getenv = f2c_string(string)
END FUNCTION f_lammps_extract_variable_getenv

FUNCTION f_lammps_extract_variable_file() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_ptr
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  TYPE(c_ptr) :: f_lammps_extract_variable_file
  CHARACTER(LEN=40) :: string

  string = lmp%extract_variable('greeting')
  f_lammps_extract_variable_file = f2c_string(string)
END FUNCTION f_lammps_extract_variable_file

FUNCTION f_lammps_extract_variable_atomfile(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: f_lammps_extract_variable_atomfile
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: atom_data

  atom_data = lmp%extract_variable('atfile')
  f_lammps_extract_variable_atomfile = atom_data(i)
END FUNCTION f_lammps_extract_variable_atomfile

FUNCTION f_lammps_extract_variable_python() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  REAL(c_double) :: f_lammps_extract_variable_python

  f_lammps_extract_variable_python = lmp%extract_variable('py')
END FUNCTION f_lammps_extract_variable_python

FUNCTION f_lammps_extract_variable_timer() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  REAL(c_double) :: f_lammps_extract_variable_timer

  f_lammps_extract_variable_timer = lmp%extract_variable('time')
END FUNCTION f_lammps_extract_variable_timer

FUNCTION f_lammps_extract_variable_internal() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  REAL(c_double) :: f_lammps_extract_variable_internal

  f_lammps_extract_variable_internal = lmp%extract_variable('int')
END FUNCTION f_lammps_extract_variable_internal

FUNCTION f_lammps_extract_variable_equal() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  REAL(c_double) :: f_lammps_extract_variable_equal

  f_lammps_extract_variable_equal = lmp%extract_variable('ex')
END FUNCTION f_lammps_extract_variable_equal

FUNCTION f_lammps_extract_variable_atom(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: f_lammps_extract_variable_atom
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: atom

  atom = lmp%extract_variable('at_z') ! z-coordinates
  f_lammps_extract_variable_atom = atom(i)
END FUNCTION f_lammps_extract_variable_atom

FUNCTION f_lammps_extract_variable_vector(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: f_lammps_extract_variable_vector
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: vector

  vector = lmp%extract_variable('center') ! z-coordinates
  f_lammps_extract_variable_vector = vector(i)
END FUNCTION f_lammps_extract_variable_vector

SUBROUTINE f_lammps_set_variable_string() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, f2c_string
  IMPLICIT NONE
  CHARACTER(LEN=40) :: string

  string = "this is the new string"
  CALL lmp%set_variable('str', string)
END SUBROUTINE f_lammps_set_variable_string

! vim: sts=2 ts=2 sw=2 et

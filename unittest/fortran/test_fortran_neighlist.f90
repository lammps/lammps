FUNCTION f_lammps_with_args() BIND(C, name="f_lammps_with_args")
  USE ISO_C_BINDING, ONLY: c_ptr
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

SUBROUTINE f_lammps_close() BIND(C, name="f_lammps_close")
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_neigh_tests() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, big_input, cont_input, pair_input
  IMPLICIT NONE

  CALL lmp%command('atom_modify map array')
  CALL lmp%commands_list(big_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('compute c all rdf 100')
  ! We create one of the fixes that requests a neighbor list, none of which
  ! is part of LAMMPS without additional packages; as such, we only do this
  ! if REPLICA is included
  IF (lmp%config_has_package('REPLICA')) THEN
    CALL lmp%command('fix f all hyper/global 1.0 0.3 0.8 300.0')
    CALL lmp%command('compute event all event/displace 1.0')
    CALL lmp%command('hyper 0 100 f event') ! using "run 0" here segfaults (?)
  ELSE
    CALL lmp%command('run 0 post no') ! otherwise neighlists won't be requested
  END IF
END SUBROUTINE f_lammps_setup_neigh_tests

FUNCTION f_lammps_pair_neighlist_test() BIND(C) RESULT(nlist_id)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: nlist_id

  nlist_id = lmp%find_pair_neighlist('lj/cut',.TRUE., 0, 0)
END FUNCTION f_lammps_pair_neighlist_test

FUNCTION f_lammps_fix_neighlist_test() BIND(C) RESULT(nlist_id)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: nlist_id

  nlist_id = lmp%find_fix_neighlist('f',0)
END FUNCTION f_lammps_fix_neighlist_test

FUNCTION f_lammps_compute_neighlist_test() BIND(C) RESULT(nlist_id)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: nlist_id

  nlist_id = lmp%find_compute_neighlist('c',0)
END FUNCTION f_lammps_compute_neighlist_test

FUNCTION f_lammps_neighlist_num_elements(id) BIND(C) RESULT(nelements)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: id
  INTEGER(c_int) :: nelements

  nelements = lmp%neighlist_num_elements(id)
END FUNCTION f_lammps_neighlist_num_elements

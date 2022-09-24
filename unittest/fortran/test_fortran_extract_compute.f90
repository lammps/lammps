MODULE keepcompute
  USE liblammps
  TYPE(LAMMPS) :: lmp
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: demo_input = &
      [ CHARACTER(len=40) ::                                &
      'region       box block 0 $x 0 3 0 4',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(LEN=40), DIMENSION(2), PARAMETER :: cont_input = &
      [ CHARACTER(len=40) ::                                &
      'create_atoms 1 single &',                            &
      ' 0.2 0.1 0.1' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: pair_input = &
      [ CHARACTER(LEN=40) ::                                &
      'pair_style lj/cut 2.5',                              &
      'pair_coeff 1 1 1.0 1.0',                             &
      'mass 1 2.0' ]
END MODULE keepcompute

FUNCTION f_lammps_with_args() BIND(C)
  USE ISO_C_BINDING, ONLY: c_ptr
  USE liblammps
  USE keepcompute, ONLY: lmp
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
  USE keepcompute, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_compute () BIND(C)
   USE LIBLAMMPS
   USE keepcompute, ONLY : lmp, demo_input, cont_input, pair_input
   IMPLICIT NONE

   CALL lmp%commands_list(demo_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(pair_input)
   CALL lmp%command("compute peratompe all pe/atom") ! per-atom vector
   call lmp%command("compute stress all stress/atom thermo_temp") ! per-atom array
   CALL lmp%command("compute COM all com") ! global vector
   CALL lmp%command("compute totalpe all reduce sum c_peratompe") ! global scalar
   CALL lmp%command("compute 1 all rdf 100") ! global array
   CALL lmp%command("compute pairdist all pair/local dist") ! local vector
   CALL lmp%command("compute pairlocal all pair/local dist dx dy dz") ! local array
END SUBROUTINE f_lammps_setup_extract_compute

FUNCTION f_lammps_extract_compute_peratom_vector (i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_double, C_int
   USE LIBLAMMPS
   USE keepcompute, ONLY : lmp
   IMPLICIT NONE
   INTEGER(C_int), INTENT(IN), VALUE :: i
   REAL(C_double) :: f_lammps_extract_compute_peratom_vector
   REAL(C_double), DIMENSION(:), POINTER :: vector => NULL()

   vector = lmp%extract_compute('peratompe', lmp%style%atom, lmp%type%vector)
   f_lammps_extract_compute_peratom_vector = vector(i)
END FUNCTION f_lammps_extract_compute_peratom_vector

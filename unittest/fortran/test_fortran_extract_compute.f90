FUNCTION f_lammps_with_args() BIND(C)
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

SUBROUTINE f_lammps_close() BIND(C)
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_extract_compute() BIND(C)
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp, big_input, cont_input, more_input, pair_input
   IMPLICIT NONE

   CALL lmp%commands_list(big_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(more_input)
   CALL lmp%commands_list(pair_input)
   CALL lmp%command("compute peratompe all pe/atom") ! per-atom vector
   call lmp%command("compute stress all stress/atom thermo_temp") ! per-atom array
   CALL lmp%command("compute totalpe all reduce sum c_peratompe") ! global scalar
   CALL lmp%command("compute COM all com") ! global vector
   CALL lmp%command("compute RDF all rdf 100") ! global array
   CALL lmp%command("compute pairdist all pair/local dist") ! local vector
   CALL lmp%command("compute pairlocal all pair/local dist dx dy dz") ! local array
   CALL lmp%command("thermo_style custom step pe c_totalpe c_COM[1]")
   CALL lmp%command("run 0") ! must be here, otherwise will SEGFAULT
END SUBROUTINE f_lammps_setup_extract_compute

FUNCTION f_lammps_extract_compute_peratom_vector(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double) :: f_lammps_extract_compute_peratom_vector
   REAL(c_double), DIMENSION(:), POINTER :: vector => NULL()

   vector = lmp%extract_compute('peratompe', lmp%style%atom, lmp%type%vector)
   f_lammps_extract_compute_peratom_vector = vector(i)
END FUNCTION f_lammps_extract_compute_peratom_vector

FUNCTION f_lammps_extract_compute_peratom_array(i,j) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i, j
   REAL(c_double) :: f_lammps_extract_compute_peratom_array
   REAL(c_double), DIMENSION(:,:), POINTER :: array => NULL()

   array = lmp%extract_compute('stress', lmp%style%atom, lmp%type%array)
   f_lammps_extract_compute_peratom_array = array(i,j)
END FUNCTION f_lammps_extract_compute_peratom_array

FUNCTION f_lammps_extract_compute_global_scalar() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   REAL(c_double) :: f_lammps_extract_compute_global_scalar
   REAL(c_double), POINTER :: scalar

   scalar = lmp%extract_compute('totalpe', lmp%style%global, lmp%type%scalar)
   f_lammps_extract_compute_global_scalar = scalar
END FUNCTION f_lammps_extract_compute_global_scalar

FUNCTION f_lammps_extract_compute_global_vector(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double) :: f_lammps_extract_compute_global_vector
   REAL(c_double), DIMENSION(:), POINTER :: vector

   vector = lmp%extract_compute('COM', lmp%style%global, lmp%type%vector)
   f_lammps_extract_compute_global_vector = vector(i)
END FUNCTION f_lammps_extract_compute_global_vector

FUNCTION f_lammps_extract_compute_global_array(i,j) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i, j
   REAL(c_double) :: f_lammps_extract_compute_global_array
   REAL(c_double), DIMENSION(:,:), POINTER :: array

   array = lmp%extract_compute('RDF', lmp%style%global, lmp%type%array)
   f_lammps_extract_compute_global_array = array(i,j)
END FUNCTION f_lammps_extract_compute_global_array

FUNCTION f_lammps_extract_compute_local_vector(i) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i
   REAL(c_double) :: f_lammps_extract_compute_local_vector
   REAL(c_double), DIMENSION(:), POINTER :: vector

   vector = lmp%extract_compute('pairdist', lmp%style%local, lmp%type%vector)
   f_lammps_extract_compute_local_vector = vector(i)
END FUNCTION f_lammps_extract_compute_local_vector

FUNCTION f_lammps_extract_compute_local_array(i, j) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double, c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int), INTENT(IN), VALUE :: i, j
   REAL(c_double) :: f_lammps_extract_compute_local_array
   REAL(c_double), DIMENSION(:,:), POINTER :: array

   array = lmp%extract_compute('pairlocal', lmp%style%local, lmp%type%array)
   f_lammps_extract_compute_local_array = array(i,j)
END FUNCTION f_lammps_extract_compute_local_array

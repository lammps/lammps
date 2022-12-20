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

SUBROUTINE f_lammps_setup_extract_global() BIND(C)
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp, big_input, cont_input, pair_input
   IMPLICIT NONE

   CALL lmp%commands_list(big_input)
   CALL lmp%commands_list(cont_input)
   CALL lmp%commands_list(pair_input)
   CALL lmp%command('run 0')
END SUBROUTINE f_lammps_setup_extract_global

SUBROUTINE f_lammps_setup_full_extract_global() BIND(C)
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTERFACE
      SUBROUTINE f_lammps_setup_extract_global() BIND(C)
      END SUBROUTINE f_lammps_setup_extract_global
   END INTERFACE

   CALL lmp%command('atom_style full')
   CALL f_lammps_setup_extract_global
   CALL lmp%command('bond_style zero')
   CALL lmp%command('angle_style zero')
   CALL lmp%command('dihedral_style zero')
   CALL lmp%command('run 0')
END SUBROUTINE f_lammps_setup_full_extract_global

FUNCTION f_lammps_extract_global_units() BIND(C) RESULT(success)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE LIBLAMMPS
   USE keepstuff, ONLY : lmp
   IMPLICIT NONE
   INTEGER(c_int) :: success
   CHARACTER(LEN=16) :: units

   ! passing strings from Fortran to C is icky, so we do the test here and
   ! report the result instead
   units = lmp%extract_global('units')
   IF (TRIM(units) == 'lj') THEN
      success = 1_c_int
   ELSE
      success = 0_c_int
   END IF
END FUNCTION f_lammps_extract_global_units

FUNCTION f_lammps_extract_global_ntimestep() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: ntimestep
   INTEGER(c_int) :: f_lammps_extract_global_ntimestep

   ntimestep = lmp%extract_global("ntimestep")
   f_lammps_extract_global_ntimestep = ntimestep
END FUNCTION f_lammps_extract_global_ntimestep

FUNCTION f_lammps_extract_global_ntimestep_big() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int64_t), POINTER :: ntimestep
   INTEGER(c_int64_t) :: f_lammps_extract_global_ntimestep_big

   ntimestep = lmp%extract_global("ntimestep")
   f_lammps_extract_global_ntimestep_big = ntimestep
END FUNCTION f_lammps_extract_global_ntimestep_big

FUNCTION f_lammps_extract_global_dt() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double), POINTER :: dt
   REAL(c_double) :: f_lammps_extract_global_dt

   dt = lmp%extract_global("dt")
   f_lammps_extract_global_dt = dt
END FUNCTION f_lammps_extract_global_dt

SUBROUTINE f_lammps_extract_global_boxlo(C_boxlo) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double), DIMENSION(3) :: C_boxlo
   REAL(c_double), DIMENSION(:), POINTER :: boxlo

   boxlo = lmp%extract_global("boxlo")
   C_boxlo = boxlo
END SUBROUTINE f_lammps_extract_global_boxlo

SUBROUTINE f_lammps_extract_global_boxhi(C_boxhi) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double), DIMENSION(3) :: C_boxhi
   REAL(c_double), DIMENSION(:), POINTER :: boxhi

   boxhi = lmp%extract_global("boxhi")
   C_boxhi = boxhi
END SUBROUTINE f_lammps_extract_global_boxhi

FUNCTION f_lammps_extract_global_boxxlo() BIND(C) RESULT(C_boxxlo)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_boxxlo
   REAL(c_double), POINTER :: boxxlo

   boxxlo = lmp%extract_global("boxxlo")
   C_boxxlo = boxxlo
END FUNCTION f_lammps_extract_global_boxxlo

FUNCTION f_lammps_extract_global_boxxhi() BIND(C) RESULT(C_boxxhi)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_boxxhi
   REAL(c_double), POINTER :: boxxhi

   boxxhi = lmp%extract_global("boxxhi")
   C_boxxhi = boxxhi
END FUNCTION f_lammps_extract_global_boxxhi

FUNCTION f_lammps_extract_global_boxylo() BIND(C) RESULT(C_boxylo)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_boxylo
   REAL(c_double), POINTER :: boxylo

   boxylo = lmp%extract_global("boxylo")
   C_boxylo = boxylo
END FUNCTION f_lammps_extract_global_boxylo

FUNCTION f_lammps_extract_global_boxyhi() BIND(C) RESULT(C_boxyhi)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_boxyhi
   REAL(c_double), POINTER :: boxyhi

   boxyhi = lmp%extract_global("boxyhi")
   C_boxyhi = boxyhi
END FUNCTION f_lammps_extract_global_boxyhi

FUNCTION f_lammps_extract_global_boxzlo() BIND(C) RESULT(C_boxzlo)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_boxzlo
   REAL(c_double), POINTER :: boxzlo

   boxzlo = lmp%extract_global("boxzlo")
   C_boxzlo = boxzlo
END FUNCTION f_lammps_extract_global_boxzlo

FUNCTION f_lammps_extract_global_boxzhi() BIND(C) RESULT(C_boxzhi)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_boxzhi
   REAL(c_double), POINTER :: boxzhi

   boxzhi = lmp%extract_global("boxzhi")
   C_boxzhi = boxzhi
END FUNCTION f_lammps_extract_global_boxzhi

SUBROUTINE f_lammps_extract_global_periodicity(C_periodicity) BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), DIMENSION(3) :: C_periodicity
   INTEGER(c_int), DIMENSION(:), POINTER :: periodicity

   periodicity = lmp%extract_global("periodicity")
   C_periodicity = periodicity
END SUBROUTINE f_lammps_extract_global_periodicity

FUNCTION f_lammps_extract_global_triclinic() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: triclinic
   INTEGER(c_int) :: f_lammps_extract_global_triclinic

   triclinic = lmp%extract_global("triclinic")
   f_lammps_extract_global_triclinic = triclinic
END FUNCTION f_lammps_extract_global_triclinic

FUNCTION f_lammps_extract_global_xy() BIND(C) RESULT(C_xy)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_xy
   REAL(c_double), POINTER :: xy

   xy = lmp%extract_global("xy")
   C_xy = xy
END FUNCTION f_lammps_extract_global_xy

FUNCTION f_lammps_extract_global_xz() BIND(C) RESULT(C_xz)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_xz
   REAL(c_double), POINTER :: xz

   xz = lmp%extract_global("xz")
   C_xz = xz
END FUNCTION f_lammps_extract_global_xz

FUNCTION f_lammps_extract_global_yz() BIND(C) RESULT(C_yz)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_yz
   REAL(c_double), POINTER :: yz

   yz = lmp%extract_global("yz")
   C_yz = yz
END FUNCTION f_lammps_extract_global_yz

FUNCTION f_lammps_extract_global_natoms() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: natoms
   INTEGER(c_int) :: f_lammps_extract_global_natoms

   natoms = lmp%extract_global("natoms")
   f_lammps_extract_global_natoms = natoms
END FUNCTION f_lammps_extract_global_natoms

FUNCTION f_lammps_extract_global_natoms_big() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int64_t), POINTER :: natoms
   INTEGER(c_int64_t) :: f_lammps_extract_global_natoms_big

   natoms = lmp%extract_global("natoms")
   f_lammps_extract_global_natoms_big = natoms
END FUNCTION f_lammps_extract_global_natoms_big

FUNCTION f_lammps_extract_global_nbonds() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: nbonds
   INTEGER(c_int) :: f_lammps_extract_global_nbonds

   nbonds = lmp%extract_global("nbonds")
   f_lammps_extract_global_nbonds = nbonds
END FUNCTION f_lammps_extract_global_nbonds

FUNCTION f_lammps_extract_global_nbonds_big() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int64_t), POINTER :: nbonds
   INTEGER(c_int64_t) :: f_lammps_extract_global_nbonds_big

   nbonds = lmp%extract_global("nbonds")
   f_lammps_extract_global_nbonds_big = nbonds
END FUNCTION f_lammps_extract_global_nbonds_big

FUNCTION f_lammps_extract_global_nangles() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: nangles
   INTEGER(c_int) :: f_lammps_extract_global_nangles

   nangles = lmp%extract_global("nangles")
   f_lammps_extract_global_nangles = nangles
END FUNCTION f_lammps_extract_global_nangles

FUNCTION f_lammps_extract_global_nangles_big() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int64_t), POINTER :: nangles
   INTEGER(c_int64_t) :: f_lammps_extract_global_nangles_big

   nangles = lmp%extract_global("nangles")
   f_lammps_extract_global_nangles_big = nangles
END FUNCTION f_lammps_extract_global_nangles_big

FUNCTION f_lammps_extract_global_ndihedrals() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: ndihedrals
   INTEGER(c_int) :: f_lammps_extract_global_ndihedrals

   ndihedrals = lmp%extract_global("ndihedrals")
   f_lammps_extract_global_ndihedrals = ndihedrals
END FUNCTION f_lammps_extract_global_ndihedrals

FUNCTION f_lammps_extract_global_ndihedrals_big() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int64_t), POINTER :: ndihedrals
   INTEGER(c_int64_t) :: f_lammps_extract_global_ndihedrals_big

   ndihedrals = lmp%extract_global("ndihedrals")
   f_lammps_extract_global_ndihedrals_big = ndihedrals
END FUNCTION f_lammps_extract_global_ndihedrals_big

FUNCTION f_lammps_extract_global_nimpropers() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: nimpropers
   INTEGER(c_int) :: f_lammps_extract_global_nimpropers

   nimpropers = lmp%extract_global("nimpropers")
   f_lammps_extract_global_nimpropers = nimpropers
END FUNCTION f_lammps_extract_global_nimpropers

FUNCTION f_lammps_extract_global_nimpropers_big() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int64_t
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int64_t), POINTER :: nimpropers
   INTEGER(c_int64_t) :: f_lammps_extract_global_nimpropers_big

   nimpropers = lmp%extract_global("nimpropers")
   f_lammps_extract_global_nimpropers_big = nimpropers
END FUNCTION f_lammps_extract_global_nimpropers_big


FUNCTION f_lammps_extract_global_ntypes() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: ntypes
   INTEGER(c_int) :: f_lammps_extract_global_ntypes

   ntypes = lmp%extract_global("ntypes")
   f_lammps_extract_global_ntypes = ntypes
END FUNCTION f_lammps_extract_global_ntypes

FUNCTION f_lammps_extract_global_nlocal() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: nlocal
   INTEGER(c_int) :: f_lammps_extract_global_nlocal

   nlocal = lmp%extract_global("nlocal")
   f_lammps_extract_global_nlocal = nlocal
END FUNCTION f_lammps_extract_global_nlocal

FUNCTION f_lammps_extract_global_nghost() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: nghost
   INTEGER(c_int) :: f_lammps_extract_global_nghost

   nghost = lmp%extract_global("nghost")
   f_lammps_extract_global_nghost = nghost
END FUNCTION f_lammps_extract_global_nghost

FUNCTION f_lammps_extract_global_nmax() BIND(C)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   INTEGER(c_int), POINTER :: nmax
   INTEGER(c_int) :: f_lammps_extract_global_nmax

   nmax = lmp%extract_global("nmax")
   f_lammps_extract_global_nmax = nmax
END FUNCTION f_lammps_extract_global_nmax

FUNCTION f_lammps_extract_global_boltz() BIND(C) RESULT(C_k_B)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_k_B
   REAL(c_double), POINTER :: k_B

   k_B = lmp%extract_global("boltz")
   C_k_B = k_B
END FUNCTION f_lammps_extract_global_boltz

FUNCTION f_lammps_extract_global_hplanck() BIND(C) RESULT(C_h)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: C_h
   REAL(c_double), POINTER :: h

   h = lmp%extract_global("boltz")
   C_h = h
END FUNCTION f_lammps_extract_global_hplanck

FUNCTION f_lammps_extract_global_angstrom() BIND(C) RESULT(Angstrom)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: Angstrom
   REAL(c_double), POINTER :: A

   A = lmp%extract_global("angstrom")
   Angstrom = A
END FUNCTION f_lammps_extract_global_angstrom

FUNCTION f_lammps_extract_global_femtosecond() BIND(C) RESULT(fs)
   USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_double
   USE keepstuff, ONLY : lmp
   USE LIBLAMMPS
   IMPLICIT NONE
   REAL(c_double) :: fs
   REAL(c_double), POINTER :: femtosecond

   femtosecond = lmp%extract_global("femtosecond")
   fs = femtosecond
END FUNCTION f_lammps_extract_global_femtosecond

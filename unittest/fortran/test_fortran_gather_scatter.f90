FUNCTION f_lammps_with_args() BIND(C)
  USE ISO_C_BINDING, ONLY: c_ptr
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
  USE ISO_C_BINDING, ONLY: c_null_ptr
  USE liblammps
  USE keepstuff, ONLY: lmp
  IMPLICIT NONE

  CALL lmp%close()
  lmp%handle = c_null_ptr
END SUBROUTINE f_lammps_close

SUBROUTINE f_lammps_setup_gather_scatter() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, big_input, cont_input, more_input, pair_input
  IMPLICIT NONE

  CALL lmp%command('atom_modify map array')
  CALL lmp%commands_list(big_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(more_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('mass 1 1.0')
  CALL lmp%command("compute pe all pe/atom")
  CALL lmp%command("fix dummy all ave/atom 1 1 1 c_pe")
END SUBROUTINE f_lammps_setup_gather_scatter

FUNCTION f_lammps_gather_atoms_mask(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  INTEGER(c_int) :: f_lammps_gather_atoms_mask
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: mask

  CALL lmp%gather_atoms('mask', 1_c_int, mask)
  f_lammps_gather_atoms_mask = mask(i)
END FUNCTION f_lammps_gather_atoms_mask

FUNCTION f_lammps_gather_atoms_position(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: f_lammps_gather_atoms_position
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: positions

  CALL lmp%gather_atoms('x', 3_c_int, positions)
  f_lammps_gather_atoms_position = positions(i)
END FUNCTION f_lammps_gather_atoms_position

FUNCTION f_lammps_gather_atoms_concat_mask(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  INTEGER(c_int) :: f_lammps_gather_atoms_concat_mask
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: mask, tag
  INTEGER :: j

  CALL lmp%gather_atoms_concat('mask', 1_c_int, mask)
  CALL lmp%gather_atoms_concat('id', 1_c_int, tag)
  f_lammps_gather_atoms_concat_mask = -1
  DO j = 1, SIZE(tag)
    IF (tag(j) == i) THEN
      f_lammps_gather_atoms_concat_mask = mask(j)
      EXIT
    END IF
  END DO
END FUNCTION f_lammps_gather_atoms_concat_mask

FUNCTION f_lammps_gather_atoms_concat_position(xyz, id) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: id, xyz
  REAL(c_double) :: f_lammps_gather_atoms_concat_position
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: positions
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: tag
  INTEGER :: j

  CALL lmp%gather_atoms_concat('x', 3_c_int, positions)
  CALL lmp%gather_atoms_concat('id', 1_c_int, tag)
  f_lammps_gather_atoms_concat_position = -1.0_c_double
  DO j = 1, SIZE(tag)
    IF (tag(j) == id) THEN
       f_lammps_gather_atoms_concat_position = positions((j-1)*3 + xyz)
    END IF
  END DO
END FUNCTION f_lammps_gather_atoms_concat_position

FUNCTION f_lammps_gather_atoms_subset_mask(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  INTEGER(c_int) :: f_lammps_gather_atoms_subset_mask
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: mask
  INTEGER :: j
  INTEGER(c_int), DIMENSION(2), PARAMETER :: tag = [3,2]

  CALL lmp%gather_atoms_subset('mask', 1_c_int, tag, mask)
  f_lammps_gather_atoms_subset_mask = -1
  DO j = 1, SIZE(tag)
    IF (tag(j) == i) THEN
      f_lammps_gather_atoms_subset_mask = mask(j)
      EXIT
    END IF
  END DO
END FUNCTION f_lammps_gather_atoms_subset_mask

FUNCTION f_lammps_gather_atoms_subset_position(xyz,id) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: id, xyz
  REAL(c_double) :: f_lammps_gather_atoms_subset_position
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: positions
  INTEGER(c_int), DIMENSION(2), PARAMETER :: tag = [3,2]
  INTEGER :: j

  CALL lmp%gather_atoms_subset('x', 3_c_int, tag, positions)
  f_lammps_gather_atoms_subset_position = -1.0_c_double
  DO j = 1, SIZE(tag)
    IF (tag(j) == id) THEN
      f_lammps_gather_atoms_subset_position = positions((j-1)*3 + xyz)
      EXIT
    END IF
  END DO
END FUNCTION f_lammps_gather_atoms_subset_position

SUBROUTINE f_lammps_scatter_atoms_masks() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: masks
  INTEGER(c_int) :: swap

  CALL lmp%gather_atoms('mask', 1_c_int, masks)

  ! swap masks of atoms 1 and 3
  swap=masks(1)
  masks(1) = masks(3)
  masks(3) = swap

  CALL lmp%scatter_atoms('mask', masks) ! push the swap back to LAMMPS
END SUBROUTINE f_lammps_scatter_atoms_masks

SUBROUTINE f_lammps_scatter_atoms_positions() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: tags
  REAL(c_double), DIMENSION(:), ALLOCATABLE, TARGET :: xvec
  REAL(c_double), DIMENSION(:,:), POINTER :: x
  REAL(c_double) :: swap(3)

  CALL lmp%gather_atoms('id',1_c_int,tags)
  CALL lmp%gather_atoms('x',3_c_int,xvec)
  x(1:3,1:SIZE(xvec)/3) => xvec

  ! swap positions of atoms 1 and 3
  swap=x(:,1)
  x(:,1) = x(:,3)
  x(:,3) = swap

  CALL lmp%scatter_atoms('x', xvec) ! push the swap back to LAMMPS
END SUBROUTINE f_lammps_scatter_atoms_positions

SUBROUTINE f_lammps_scatter_atoms_subset_mask() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: all_masks
  INTEGER(c_int), DIMENSION(2), PARAMETER :: tags = [3,1]
  INTEGER(c_int), DIMENSION(2) :: masks

  CALL lmp%gather_atoms('mask', 1_c_int, all_masks)

  ! swap masks of atoms 1 and 3 in the new array (because 'tags' is reversed)
  masks(1) = all_masks(1)
  masks(2) = all_masks(3)

  CALL lmp%scatter_atoms_subset('mask', tags, masks) ! push the swap to LAMMPS
END SUBROUTINE f_lammps_scatter_atoms_subset_mask

SUBROUTINE f_lammps_setup_gather_bonds() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp, cont_input, more_input, pair_input
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE f_lammps_setup_gather_scatter() BIND(C)
    END SUBROUTINE f_lammps_setup_gather_scatter
  END INTERFACE

  CALL lmp%command('atom_modify map array')
  CALL lmp%command('atom_style full')
  CALL lmp%command('region simbox block 0 4 0 5 0 4')
  CALL lmp%command('create_box 1 simbox bond/types 1 extra/bond/per/atom 2')
  CALL lmp%command('create_atoms 1 single 1.0 1.0 ${zpos}')
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(more_input)
  CALL lmp%commands_list(pair_input)
  CALL lmp%command('bond_style zero')
  CALL lmp%command('bond_coeff *')
  CALL lmp%command('create_bonds many all all 1 0.0 1.5')
  CALL lmp%command('run 0')
END SUBROUTINE f_lammps_setup_gather_bonds

FUNCTION f_lammps_test_gather_bonds_small() BIND(C) RESULT(success)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: success
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE, TARGET :: bonds
  INTEGER(c_int), DIMENSION(:,:), POINTER :: bonds_array

  CALL lmp%gather_bonds(bonds)
  bonds_array(1:3,1:SIZE(bonds)/3) => bonds
  IF ( ALL(bonds_array(:,1) == [INTEGER(c_int) :: 1,1,3]) &
      .AND. ALL(bonds_array(:,2) == [INTEGER(c_int) :: 1,2,3])) THEN
    success = 1_c_int
  ELSE
    success = 0_c_int
  END IF
END FUNCTION f_lammps_test_gather_bonds_small

FUNCTION f_lammps_test_gather_bonds_big() BIND(C) RESULT(success)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int) :: success
  INTEGER(c_int64_t), DIMENSION(:), ALLOCATABLE, TARGET :: bonds
  INTEGER(c_int64_t), DIMENSION(:,:), POINTER :: bonds_array

  CALL lmp%gather_bonds(bonds)
  bonds_array(1:3,1:SIZE(bonds)/3) => bonds
  IF ( ALL(bonds_array(:,1) == [INTEGER(c_int64_t) :: 1,1,3]) &
      .AND. ALL(bonds_array(:,2) == [INTEGER(c_int64_t) :: 1,2,3])) THEN
    success = 1_c_int
  ELSE
    success = 0_c_int
  END IF
END FUNCTION f_lammps_test_gather_bonds_big

FUNCTION f_lammps_gather_pe_atom(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: f_lammps_gather_pe_atom
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: pe_atom

  CALL lmp%gather('c_pe', 1_c_int, pe_atom)
  f_lammps_gather_pe_atom = pe_atom(i)
END FUNCTION f_lammps_gather_pe_atom

FUNCTION f_lammps_gather_pe_atom_concat(i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN), VALUE :: i
  REAL(c_double) :: f_lammps_gather_pe_atom_concat
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: pe_atom
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE :: tag
  INTEGER :: j

  CALL lmp%gather_concat('id', 1_c_int, tag)
  CALL lmp%gather_concat('c_pe', 1_c_int, pe_atom)
  DO j = 1, SIZE(tag)
    IF (tag(j) == i) THEN
       f_lammps_gather_pe_atom_concat = pe_atom(j)
       EXIT
    END IF
  END DO
  f_lammps_gather_pe_atom_concat = pe_atom(i)
END FUNCTION f_lammps_gather_pe_atom_concat

SUBROUTINE f_lammps_gather_pe_atom_subset(ids, pe) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), INTENT(IN) :: ids(2)
  REAL(c_double), INTENT(OUT) :: pe(2)
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: pe_atom
  INTEGER(c_int) :: natoms

  natoms = NINT(lmp%get_natoms(), c_int)
  CALL lmp%gather_subset('c_pe', 1, ids, pe_atom)
  pe(1:2) = pe_atom(1:2)
END SUBROUTINE f_lammps_gather_pe_atom_subset

SUBROUTINE f_lammps_scatter_compute() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: pe_atom
  REAL(c_double) :: swap

  CALL lmp%gather('c_pe', 1_c_int, pe_atom)

  ! swap the computed energy of atoms 1 and 3
  swap = pe_atom(1)
  pe_atom(1) = pe_atom(3)
  pe_atom(3) = swap

  CALL lmp%scatter('c_pe', pe_atom) ! push the swap back to LAMMPS
END SUBROUTINE f_lammps_scatter_compute

SUBROUTINE f_lammps_scatter_subset_compute() BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER(c_int), PARAMETER :: ids(2) = [3,1]
  REAL(c_double), DIMENSION(:), ALLOCATABLE :: pe_atom
  REAL(c_double) :: swap

  CALL lmp%gather_subset('c_pe', 1_c_int, ids, pe_atom)

  ! swap the computed energy of atoms 1 and 3
  swap = pe_atom(1)
  pe_atom(1) = pe_atom(2)
  pe_atom(2) = swap

  CALL lmp%scatter_subset('c_pe', ids, pe_atom) ! push the swap back to LAMMPS
END SUBROUTINE f_lammps_scatter_subset_compute

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
  USE keepstuff, ONLY : lmp, big_input, cont_input, more_input
  IMPLICIT NONE

  CALL lmp%command('atom_modify map array')
  CALL lmp%commands_list(big_input)
  CALL lmp%commands_list(cont_input)
  CALL lmp%commands_list(more_input)
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
  INTEGER(c_int), DIMENSION(*), PARAMETER :: tag = [3,2]

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
  INTEGER(c_int), DIMENSION(*), PARAMETER :: tag = [3,2]
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
  INTEGER(c_int), DIMENSION(*), PARAMETER :: tags = [3,1]
  INTEGER(c_int), DIMENSION(2) :: masks

  CALL lmp%gather_atoms('mask', 1_c_int, all_masks)

  ! swap masks of atoms 1 and 3 in the new array (because 'tags' is reversed)
  masks(1) = all_masks(1)
  masks(2) = all_masks(3)

  CALL lmp%scatter_atoms_subset('mask', tags, masks) ! push the swap to LAMMPS
END SUBROUTINE f_lammps_scatter_atoms_subset_mask

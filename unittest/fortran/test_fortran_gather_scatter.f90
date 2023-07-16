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

SUBROUTINE f_lammps_setup_gather_topology() BIND(C)
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE

  CALL lmp%command('include ${input_dir}/in.fourmol')
  CALL lmp%command('run 0 post no')
END SUBROUTINE f_lammps_setup_gather_topology

FUNCTION f_lammps_test_gather_bonds_small() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, nbonds, size_bigint
  INTEGER(c_int) :: count
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE, TARGET :: bonds
  INTEGER(c_int), DIMENSION(:,:), POINTER :: bonds_array
  INTEGER(c_int), POINTER :: nbonds_small
  INTEGER(c_int64_t), POINTER :: nbonds_big

  size_bigint = lmp%extract_setting('bigint')
  IF (size_bigint == 4) THEN
      nbonds_small = lmp%extract_global('nbonds')
      nbonds = nbonds_small
  ELSE
      nbonds_big = lmp%extract_global('nbonds')
      nbonds = INT(nbonds_big)
  END IF

  CALL lmp%gather_bonds(bonds)
  bonds_array(1:3,1:SIZE(bonds)/3) => bonds
  count = 0
  DO i=1, nbonds
      count = count + check_bond(i, 5, 1, 2, bonds_array)
      count = count + check_bond(i, 3, 1, 3, bonds_array)
      count = count + check_bond(i, 2, 3, 4, bonds_array)
      count = count + check_bond(i, 2, 3, 5, bonds_array)
      count = count + check_bond(i, 1, 3, 6, bonds_array)
      count = count + check_bond(i, 3, 6, 8, bonds_array)
      count = count + check_bond(i, 4, 6, 7, bonds_array)
      count = count + check_bond(i, 5, 8, 9, bonds_array)
      count = count + check_bond(i, 5, 27, 28, bonds_array)
      count = count + check_bond(i, 5, 27, 29, bonds_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_bond(idx, batom1, batom2, btype, barray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, batom1, batom2, btype
    INTEGER(c_int), DIMENSION(:,:) :: barray
    check_bond = 0
    IF ((barray(1,idx) == batom1) .AND. (barray(2,idx) == batom2)) THEN
        IF (barray(3,idx) == btype) check_bond = 1
    END IF
    IF ((barray(1,idx) == batom2) .AND. (barray(2,idx) == batom1)) THEN
        IF (barray(3,idx) == btype) check_bond = 1
    END IF
  END FUNCTION check_bond

END FUNCTION f_lammps_test_gather_bonds_small

FUNCTION f_lammps_test_gather_bonds_big() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, nbonds
  INTEGER(c_int) :: count
  INTEGER(c_int64_t), DIMENSION(:), ALLOCATABLE, TARGET :: bonds
  INTEGER(c_int64_t), DIMENSION(:,:), POINTER :: bonds_array
  INTEGER(c_int64_t), POINTER :: nbonds_big

  nbonds_big = lmp%extract_global('nbonds')
  nbonds = INT(nbonds_big)
  CALL lmp%gather_bonds(bonds)
  bonds_array(1:3,1:SIZE(bonds)/3) => bonds
  count = 0
  DO i=1, nbonds
      count = count + check_bond(i, 5, 1, 2, bonds_array)
      count = count + check_bond(i, 3, 1, 3, bonds_array)
      count = count + check_bond(i, 2, 3, 4, bonds_array)
      count = count + check_bond(i, 2, 3, 5, bonds_array)
      count = count + check_bond(i, 1, 3, 6, bonds_array)
      count = count + check_bond(i, 3, 6, 8, bonds_array)
      count = count + check_bond(i, 4, 6, 7, bonds_array)
      count = count + check_bond(i, 5, 8, 9, bonds_array)
      count = count + check_bond(i, 5, 27, 28, bonds_array)
      count = count + check_bond(i, 5, 27, 29, bonds_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_bond(idx, batom1, batom2, btype, barray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, batom1, batom2, btype
    INTEGER(c_int64_t), DIMENSION(:,:) :: barray
    check_bond = 0
    IF ((barray(1,idx) == batom1) .AND. (barray(2,idx) == batom2)) THEN
        IF (barray(3,idx) == btype) check_bond = 1
    END IF
    IF ((barray(1,idx) == batom2) .AND. (barray(2,idx) == batom1)) THEN
        IF (barray(3,idx) == btype) check_bond = 1
    END IF
  END FUNCTION check_bond

END FUNCTION f_lammps_test_gather_bonds_big

FUNCTION f_lammps_test_gather_angles_small() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, nangles, size_bigint
  INTEGER(c_int) :: count
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE, TARGET :: angles
  INTEGER(c_int), DIMENSION(:,:), POINTER :: angles_array
  INTEGER(c_int), POINTER :: nangles_small
  INTEGER(c_int64_t), POINTER :: nangles_big

  size_bigint = lmp%extract_setting('bigint')
  IF (size_bigint == 4) THEN
      nangles_small = lmp%extract_global('nangles')
      nangles = nangles_small
  ELSE
      nangles_big = lmp%extract_global('nangles')
      nangles = INT(nangles_big)
  END IF

  CALL lmp%gather_angles(angles)
  angles_array(1:4,1:SIZE(angles)/4) => angles
  count = 0
  DO i=1, nangles
      count = count + check_angle(i, 4, 2, 1, 3, angles_array)
      count = count + check_angle(i, 4, 1, 3, 5, angles_array)
      count = count + check_angle(i, 4, 1, 3, 4, angles_array)
      count = count + check_angle(i, 4, 13, 12, 15, angles_array)
      count = count + check_angle(i, 4, 13, 12, 14, angles_array)
      count = count + check_angle(i, 2, 5, 3, 6, angles_array)
      count = count + check_angle(i, 2, 4, 3, 6, angles_array)
      count = count + check_angle(i, 3, 3, 6, 7, angles_array)
      count = count + check_angle(i, 3, 3, 6, 8, angles_array)
      count = count + check_angle(i, 1, 22, 21, 23, angles_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_angle(idx, aatom1, aatom2, aatom3, atype, aarray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, aatom1, aatom2, aatom3, atype
    INTEGER(c_int), DIMENSION(:,:) :: aarray
    check_angle = 0
    IF ((aarray(1,idx) == aatom1) .AND. (aarray(2,idx) == aatom2) .AND. (aarray(3,idx) == aatom3)) THEN
        IF (aarray(4,idx) == atype) check_angle = 1
    END IF
    IF ((aarray(1,idx) == aatom3) .AND. (aarray(2,idx) == aatom2) .AND. (aarray(3,idx) == aatom1)) THEN
        IF (aarray(4,idx) == atype) check_angle = 1
    END IF
  END FUNCTION check_angle

END FUNCTION f_lammps_test_gather_angles_small

FUNCTION f_lammps_test_gather_angles_big() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, nangles
  INTEGER(c_int) :: count
  INTEGER(c_int64_t), DIMENSION(:), ALLOCATABLE, TARGET :: angles
  INTEGER(c_int64_t), DIMENSION(:,:), POINTER :: angles_array
  INTEGER(c_int64_t), POINTER :: nangles_big

  nangles_big = lmp%extract_global('nangles')
  nangles = INT(nangles_big)
  CALL lmp%gather_angles(angles)
  angles_array(1:4,1:SIZE(angles)/4) => angles
  count = 0
  DO i=1, nangles
      count = count + check_angle(i, 4, 2, 1, 3, angles_array)
      count = count + check_angle(i, 4, 1, 3, 5, angles_array)
      count = count + check_angle(i, 4, 1, 3, 4, angles_array)
      count = count + check_angle(i, 4, 13, 12, 15, angles_array)
      count = count + check_angle(i, 4, 13, 12, 14, angles_array)
      count = count + check_angle(i, 2, 5, 3, 6, angles_array)
      count = count + check_angle(i, 2, 4, 3, 6, angles_array)
      count = count + check_angle(i, 3, 3, 6, 7, angles_array)
      count = count + check_angle(i, 3, 3, 6, 8, angles_array)
      count = count + check_angle(i, 1, 22, 21, 23, angles_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_angle(idx, aatom1, aatom2, aatom3, atype, aarray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, aatom1, aatom2, aatom3, atype
    INTEGER(c_int64_t), DIMENSION(:,:) :: aarray
    check_angle = 0
    IF ((aarray(1,idx) == aatom1) .AND. (aarray(2,idx) == aatom2) .AND. (aarray(3,idx) == aatom3)) THEN
        IF (aarray(4,idx) == atype) check_angle = 1
    END IF
    IF ((aarray(1,idx) == aatom3) .AND. (aarray(2,idx) == aatom2) .AND. (aarray(3,idx) == aatom1)) THEN
        IF (aarray(4,idx) == atype) check_angle = 1
    END IF
  END FUNCTION check_angle

END FUNCTION f_lammps_test_gather_angles_big

FUNCTION f_lammps_test_gather_dihedrals_small() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, ndihedrals, size_bigint
  INTEGER(c_int) :: count
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE, TARGET :: dihedrals
  INTEGER(c_int), DIMENSION(:,:), POINTER :: dihedrals_array
  INTEGER(c_int), POINTER :: ndihedrals_small
  INTEGER(c_int64_t), POINTER :: ndihedrals_big

  size_bigint = lmp%extract_setting('bigint')
  IF (size_bigint == 4) THEN
      ndihedrals_small = lmp%extract_global('ndihedrals')
      ndihedrals = ndihedrals_small
  ELSE
      ndihedrals_big = lmp%extract_global('ndihedrals')
      ndihedrals = INT(ndihedrals_big)
  END IF

  CALL lmp%gather_dihedrals(dihedrals)
  dihedrals_array(1:5,1:SIZE(dihedrals)/5) => dihedrals
  count = 0
  DO i=1, ndihedrals
      count = count + check_dihedral(i, 2, 2, 1, 3, 6, dihedrals_array)
      count = count + check_dihedral(i, 2, 2, 1, 3, 4, dihedrals_array)
      count = count + check_dihedral(i, 3, 2, 1, 3, 5, dihedrals_array)
      count = count + check_dihedral(i, 1, 1, 3, 6, 8, dihedrals_array)
      count = count + check_dihedral(i, 1, 1, 3, 6, 7, dihedrals_array)
      count = count + check_dihedral(i, 5, 4, 3, 6, 8, dihedrals_array)
      count = count + check_dihedral(i, 5, 4, 3, 6, 7, dihedrals_array)
      count = count + check_dihedral(i, 5, 16, 10, 12, 13, dihedrals_array)
      count = count + check_dihedral(i, 5, 16, 10, 12, 14, dihedrals_array)
      count = count + check_dihedral(i, 5, 16, 10, 12, 15, dihedrals_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_dihedral(idx, datom1, datom2, datom3, datom4, dtype, darray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, datom1, datom2, datom3, datom4, dtype
    INTEGER(c_int), DIMENSION(:,:) :: darray
    check_dihedral = 0
    IF ((darray(1,idx) == datom1) .AND. (darray(2,idx) == datom2) &
        .AND. (darray(3,idx) == datom3) .AND. (darray(4,idx) == datom4)) THEN
        IF (darray(5,idx) == dtype) check_dihedral = 1
    END IF
  END FUNCTION check_dihedral

END FUNCTION f_lammps_test_gather_dihedrals_small

FUNCTION f_lammps_test_gather_dihedrals_big() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, ndihedrals
  INTEGER(c_int) :: count
  INTEGER(c_int64_t), DIMENSION(:), ALLOCATABLE, TARGET :: dihedrals
  INTEGER(c_int64_t), DIMENSION(:,:), POINTER :: dihedrals_array
  INTEGER(c_int64_t), POINTER :: ndihedrals_big

  ndihedrals_big = lmp%extract_global('ndihedrals')
  ndihedrals = INT(ndihedrals_big)
  CALL lmp%gather_dihedrals(dihedrals)
  dihedrals_array(1:5,1:SIZE(dihedrals)/5) => dihedrals
  count = 0
  DO i=1, ndihedrals
      count = count + check_dihedral(i, 2, 2, 1, 3, 6, dihedrals_array)
      count = count + check_dihedral(i, 2, 2, 1, 3, 4, dihedrals_array)
      count = count + check_dihedral(i, 3, 2, 1, 3, 5, dihedrals_array)
      count = count + check_dihedral(i, 1, 1, 3, 6, 8, dihedrals_array)
      count = count + check_dihedral(i, 1, 1, 3, 6, 7, dihedrals_array)
      count = count + check_dihedral(i, 5, 4, 3, 6, 8, dihedrals_array)
      count = count + check_dihedral(i, 5, 4, 3, 6, 7, dihedrals_array)
      count = count + check_dihedral(i, 5, 16, 10, 12, 13, dihedrals_array)
      count = count + check_dihedral(i, 5, 16, 10, 12, 14, dihedrals_array)
      count = count + check_dihedral(i, 5, 16, 10, 12, 15, dihedrals_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_dihedral(idx, datom1, datom2, datom3, datom4, dtype, darray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, datom1, datom2, datom3, datom4, dtype
    INTEGER(c_int64_t), DIMENSION(:,:) :: darray
    check_dihedral = 0
    IF ((darray(1,idx) == datom1) .AND. (darray(2,idx) == datom2) &
        .AND. (darray(3,idx) == datom3) .AND. (darray(4,idx) == datom4)) THEN
        IF (darray(5,idx) == dtype) check_dihedral = 1
    END IF
  END FUNCTION check_dihedral

END FUNCTION f_lammps_test_gather_dihedrals_big

FUNCTION f_lammps_test_gather_impropers_small() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, nimpropers, size_bigint
  INTEGER(c_int) :: count
  INTEGER(c_int), DIMENSION(:), ALLOCATABLE, TARGET :: impropers
  INTEGER(c_int), DIMENSION(:,:), POINTER :: impropers_array
  INTEGER(c_int), POINTER :: nimpropers_small
  INTEGER(c_int64_t), POINTER :: nimpropers_big

  size_bigint = lmp%extract_setting('bigint')
  IF (size_bigint == 4) THEN
      nimpropers_small = lmp%extract_global('nimpropers')
      nimpropers = nimpropers_small
  ELSE
      nimpropers_big = lmp%extract_global('nimpropers')
      nimpropers = INT(nimpropers_big)
  END IF

  CALL lmp%gather_impropers(impropers)
  impropers_array(1:5,1:SIZE(impropers)/5) => impropers
  count = 0
  DO i=1, nimpropers
      count = count + check_improper(i, 1, 6, 3, 8, 7, impropers_array)
      count = count + check_improper(i, 2, 8, 6, 10, 9, impropers_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_improper(idx, datom1, datom2, datom3, datom4, dtype, darray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, datom1, datom2, datom3, datom4, dtype
    INTEGER(c_int), DIMENSION(:,:) :: darray
    check_improper = 0
    IF ((darray(1,idx) == datom1) .AND. (darray(2,idx) == datom2) &
        .AND. (darray(3,idx) == datom3) .AND. (darray(4,idx) == datom4)) THEN
        IF (darray(5,idx) == dtype) check_improper = 1
    END IF
  END FUNCTION check_improper

END FUNCTION f_lammps_test_gather_impropers_small

FUNCTION f_lammps_test_gather_impropers_big() BIND(C) RESULT(count)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_int64_t
  USE LIBLAMMPS
  USE keepstuff, ONLY : lmp
  IMPLICIT NONE
  INTEGER :: i, nimpropers
  INTEGER(c_int) :: count
  INTEGER(c_int64_t), DIMENSION(:), ALLOCATABLE, TARGET :: impropers
  INTEGER(c_int64_t), DIMENSION(:,:), POINTER :: impropers_array
  INTEGER(c_int64_t), POINTER :: nimpropers_big

  nimpropers_big = lmp%extract_global('nimpropers')
  nimpropers = INT(nimpropers_big)
  CALL lmp%gather_impropers(impropers)
  impropers_array(1:5,1:SIZE(impropers)/5) => impropers
  count = 0
  DO i=1, nimpropers
      count = count + check_improper(i, 1, 6, 3, 8, 7, impropers_array)
      count = count + check_improper(i, 2, 8, 6, 10, 9, impropers_array)
  END DO

CONTAINS

  INTEGER FUNCTION check_improper(idx, datom1, datom2, datom3, datom4, dtype, darray)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx, datom1, datom2, datom3, datom4, dtype
    INTEGER(c_int64_t), DIMENSION(:,:) :: darray
    check_improper = 0
    IF ((darray(1,idx) == datom1) .AND. (darray(2,idx) == datom2) &
        .AND. (darray(3,idx) == datom3) .AND. (darray(4,idx) == datom4)) THEN
        IF (darray(5,idx) == dtype) check_improper = 1
    END IF
  END FUNCTION check_improper

END FUNCTION f_lammps_test_gather_impropers_big

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

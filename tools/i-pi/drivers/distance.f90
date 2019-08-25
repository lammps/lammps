! This contains the algorithms needed to calculate the distance between atoms.
! 
! Copyright (C) 2013, Joshua More and Michele Ceriotti
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!
! Functions:
!    vector_separation: Calculates the vector separating two atoms.
!    separation: Calculates the square distance between two vectors.
!    nearest_neighbours: Generates arrays to calculate the pairs of atoms within
!       a certain radius of each other.

      MODULE DISTANCE
      IMPLICIT NONE

      CONTAINS

         SUBROUTINE vector_separation(cell_h, cell_ih, ri, rj, rij, r2)
            ! Calculates the vector separating two atoms.
            !
            ! Note that minimum image convention is used, so only the image of
            ! atom j that is the shortest distance from atom i is considered.
            !
            ! Also note that while this may not work if the simulation
            ! box is highly skewed from orthorhombic, as
            ! in this case it is possible to return a distance less than the
            ! nearest neighbour distance. However, this will not be of 
            ! importance unless the cut-off radius is more than half the 
            ! width of the shortest face-face distance of the simulation box,
            ! which should never be the case.
            !
            ! Args:
            !    cell_h: The simulation box cell vector matrix.
            !    cell_ih: The inverse of the simulation box cell vector matrix.
            !    ri: The position vector of atom i.
            !    rj: The position vector of atom j
            !    rij: The vector separating atoms i and j.
            !    r2: The square of the distance between atoms i and j.

            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: rij 
            DOUBLE PRECISION, INTENT(OUT) :: r2

            INTEGER k
            DOUBLE PRECISION, DIMENSION(3) :: sij  
            ! The separation in a basis where the simulation box  
            ! is a unit cube.
                                                   
            sij = matmul(cell_ih, ri - rj)
            DO k = 1, 3
               ! Finds the smallest separation of all the images of atom i and j
               sij(k) = sij(k) - dnint(sij(k)) 
            ENDDO
            rij = matmul(cell_h, sij)
            r2 = dot_product(rij,rij)

         END SUBROUTINE

         SUBROUTINE separation(cell_h, cell_ih, ri, rj, r2)
            ! Calculates the squared distance between two position vectors.
            !
            ! Note that minimum image convention is used, so only the image of
            ! atom j that is the shortest distance from atom i is considered.
            !
            ! Also note that while this may not work if the simulation
            ! box is highly skewed from orthorhombic, as
            ! in this case it is possible to return a distance less than the
            ! nearest neighbour distance. However, this will not be of 
            ! importance unless the cut-off radius is more than half the 
            ! width of the shortest face-face distance of the simulation box,
            ! which should never be the case.
            !
            ! Args:
            !    cell_h: The simulation box cell vector matrix.
            !    cell_ih: The inverse of the simulation box cell vector matrix.
            !    ri: The position vector of atom i.
            !    rj: The position vector of atom j
            !    r2: The square of the distance between atoms i and j.

            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri 
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, INTENT(OUT) :: r2

            INTEGER k
            ! The separation in a basis where the simulation box  
            ! is a unit cube.
            DOUBLE PRECISION, DIMENSION(3) :: sij
            DOUBLE PRECISION, DIMENSION(3) :: rij
                                                   
            sij = matmul(cell_ih, ri - rj)
            DO k = 1, 3
               ! Finds the smallest separation of all the images of atom i and j
               sij(k) = sij(k) - dnint(sij(k)) 
            ENDDO
            rij = matmul(cell_h, sij)
            r2 = dot_product(rij, rij)

         END SUBROUTINE

         SUBROUTINE nearest_neighbours(rn, natoms, atoms, cell_h, cell_ih, index_list, n_list)
            ! Creates a list of all the pairs of atoms that are closer together
            ! than a certain distance.
            !
            ! This takes all the positions, and calculates which ones are
            ! shorter than the distance rn. This creates two vectors, index_list
            ! and n_list. index_list(i) gives the last index of n_list that
            ! corresponds to a neighbour of atom i.
            ! 
            !
            ! Args:
            !    rn: The nearest neighbour list cut-off parameter. This should
            !       be larger than the potential cut-off radius.
            !    natoms: The number of atoms in the system.
            !    atoms: A vector holding all the atom positions.
            !    cell_h: The simulation box cell vector matrix.
            !    cell_ih: The inverse of the simulation box cell vector matrix.
            !    index_list: A array giving the last index of n_list that 
            !       gives the neighbours of a given atom. Essentially keeps 
            !       track of how many atoms neighbour a given atom.
            !    n_list: An array giving the indices of the atoms that neighbour
            !       the atom determined by index_list. Essentially keeps track
            !       of which atoms neighbour a given atom.

            DOUBLE PRECISION, INTENT(IN) :: rn
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: atoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            INTEGER, DIMENSION(natoms), INTENT(OUT) :: index_list
            INTEGER, DIMENSION(natoms*(natoms-1)/2), INTENT(OUT) :: n_list

            INTEGER :: i, j
            DOUBLE PRECISION r2

         index_list(1) = 0

         DO i = 1, natoms - 1
            DO j = i + 1, natoms
               CALL separation(cell_h, cell_ih, atoms(i,:), atoms(j,:), r2)
               IF (r2 < rn*rn) THEN
                  ! We have found an atom that neighbours atom i, so the
                  ! i-th index of index_list is incremented by one, and a new
                  ! entry is added to n_list.
                  index_list(i) = index_list(i) + 1
                  n_list(index_list(i)) = j
               ENDIF
            ENDDO
            index_list(i+1) = index_list(i)
         ENDDO 

         END SUBROUTINE

      END MODULE

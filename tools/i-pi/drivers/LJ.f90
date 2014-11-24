! This performs the calculations necessary to run a Lennard-Jones (LJ) 
! simulation.
! 
! Copyright (C) 2013, Joshua More and Michele Ceriotti
! 
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! This contains the functions that calculate the potential, forces and
! virial tensor of a single-component LJ system.
! Includes functions which calculate the long-range correction terms for a
! simulation with a sharp nearest-neighbour cut-off.
!
! Functions:
!    LJ_functions: Calculates the LJ pair potential and the magnitude of the
!       forces acting on a pair of atoms.
!    LJ_fij: Calculates the LJ pair potential and force vector for the
!       interaction of a pair of atoms.
!    LJ_longrange: Calculates the long range correction to the potential
!       and virial.
!    LJ_getall: Calculates the potential and virial of the system and the
!       forces acting on all the atoms.

      MODULE LJ
         USE DISTANCE
      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: four_tau_by_3 = 8.3775804095727811d0 

      CONTAINS

         SUBROUTINE LJ_functions(sigma, eps, r, pot, force)
            ! Calculates the magnitude of the LJ force and potential between
            ! a pair of atoms at a given distance from each other.
            !
            ! Args:
            !    sigma: The LJ distance parameter.
            !    eps: The LJ energy parameter.
            !    r: The separation of the atoms.
            !    pot: The LJ interaction potential.
            !    force: The magnitude of the LJ force.

            DOUBLE PRECISION, INTENT(IN) :: sigma
            DOUBLE PRECISION, INTENT(IN) :: eps
            DOUBLE PRECISION, INTENT(IN) :: r
            DOUBLE PRECISION, INTENT(OUT) :: pot
            DOUBLE PRECISION, INTENT(OUT) :: force

            DOUBLE PRECISION sigma_by_r6

            sigma_by_r6 = sigma/r
            sigma_by_r6 = sigma_by_r6*sigma_by_r6*sigma_by_r6
            sigma_by_r6 = sigma_by_r6*sigma_by_r6

            pot = 4*eps*(sigma_by_r6*(sigma_by_r6 - 1))
            force = 24*eps*(sigma_by_r6*(2*sigma_by_r6 - 1)/r)

         END SUBROUTINE

         SUBROUTINE LJ_fij(sigma, eps, rij, r, pot, fij)
            ! This calculates the LJ potential energy and the magnitude and
            ! direction of the force acting on a pair of atoms.
            !
            ! Args:
            !    sigma: The LJ distance parameter.
            !    eps: The LJ energy parameter.
            !    rij: The vector joining the two atoms.
            !    r: The separation of the two atoms.
            !    pot: The LJ interaction potential.
            !    fij: The LJ force vector.

            DOUBLE PRECISION, INTENT(IN) :: sigma
            DOUBLE PRECISION, INTENT(IN) :: eps
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rij
            DOUBLE PRECISION, INTENT(IN) :: r
            DOUBLE PRECISION, INTENT(OUT) :: pot
            DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: fij
   
            DOUBLE PRECISION f_tot
   
            CALL LJ_functions(sigma, eps, r, pot, f_tot)
            fij = f_tot*rij/r
   
         END SUBROUTINE

         SUBROUTINE LJ_longrange(rc, sigma, eps, natoms, volume, pot_lr, vir_lr)
            ! Calculates the long range correction to the total potential and
            ! virial pressure.
            !
            ! Uses the tail correction for a sharp cut-off, with no smoothing
            ! function, as derived in Martyna and Hughes, Journal of Chemical
            ! Physics, 110, 3275, (1999).
            !
            ! Args:
            !    rc: The cut-off radius.
            !    sigma: The LJ distance parameter.
            !    eps: The LJ energy parameter.
            !    natoms: The number of atoms in the system.
            !    volume: The volume of the system box.
            !    pot_lr: The tail correction to the LJ interaction potential.
            !    vir_lr: The tail correction to the LJ virial pressure.

            DOUBLE PRECISION, INTENT(IN) :: rc
            DOUBLE PRECISION, INTENT(IN) :: sigma
            DOUBLE PRECISION, INTENT(IN) :: eps
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, INTENT(IN) :: volume
            DOUBLE PRECISION, INTENT(OUT) :: pot_lr
            DOUBLE PRECISION, INTENT(OUT) :: vir_lr

            DOUBLE PRECISION sbyr, s3byr3, s6byr3, s6byr6, prefactor

            sbyr = sigma/rc
            s3byr3 = sbyr*sbyr*sbyr
            s6byr6 = s3byr3*s3byr3
            prefactor = four_tau_by_3*natoms*natoms*eps/volume
            prefactor = prefactor*s3byr3*sigma*sigma*sigma

            pot_lr = prefactor*(s6byr6/3-1)
            vir_lr = prefactor*(s6byr6-1) + pot_lr

         END SUBROUTINE

         SUBROUTINE LJ_getall(rc, sigma, eps, natoms, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
            ! Calculates the LJ potential energy and virial and the forces 
            ! acting on all the atoms.
            !
            ! Args:
            !    rc: The cut-off radius.
            !    sigma: The LJ distance parameter.
            !    eps: The LJ energy parameter.
            !    natoms: The number of atoms in the system.
            !    atoms: A vector holding all the atom positions.
            !    cell_h: The simulation box cell vector matrix.
            !    cell_ih: The inverse of the simulation box cell vector matrix.
            !    index_list: A array giving the last index of n_list that 
            !       gives the neighbours of a given atom.
            !    n_list: An array giving the indices of the atoms that neighbour
            !       the atom determined by index_list.
            !    pot: The total potential energy of the system.
            !    forces: An array giving the forces acting on all the atoms.
            !    virial: The virial tensor, not divided by the volume.

            DOUBLE PRECISION, INTENT(IN) :: rc
            DOUBLE PRECISION, INTENT(IN) :: sigma
            DOUBLE PRECISION, INTENT(IN) :: eps
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(natoms,3), INTENT(IN) :: atoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            INTEGER, DIMENSION(natoms), INTENT(IN) :: index_list
            INTEGER, DIMENSION(natoms*(natoms-1)/2), INTENT(IN) :: n_list
            DOUBLE PRECISION, INTENT(OUT) :: pot
            DOUBLE PRECISION, DIMENSION(natoms,3), INTENT(OUT) :: forces
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: virial

            INTEGER i, j, k, l, start
            DOUBLE PRECISION, DIMENSION(3) :: fij, rij
            DOUBLE PRECISION r2, pot_ij, pot_lr, vir_lr, volume

            forces = 0.0d0
            pot = 0.0d0
            virial = 0.0d0

            start = 1

            DO i = 1, natoms - 1
               ! Only loops over the neigbour list, not all the atoms.
               DO j = start, index_list(i)
                  CALL vector_separation(cell_h, cell_ih, atoms(i,:), atoms(n_list(j),:), rij, r2)
                  IF (r2 < rc*rc) THEN ! Only calculates contributions between neighbouring particles.
                     CALL LJ_fij(sigma, eps, rij, sqrt(r2), pot_ij, fij)

                     forces(i,:) = forces(i,:) + fij
                     forces(n_list(j),:) = forces(n_list(j),:) - fij
                     pot = pot + pot_ij
                     DO k = 1, 3
                        DO l = k, 3
                           ! Only the upper triangular elements calculated.
                           virial(k,l) = virial(k,l) + fij(k)*rij(l)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               start = index_list(i) + 1
            ENDDO

            ! Assuming an upper-triangular vector matrix for the simulation box.
            volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)
            CALL LJ_longrange(rc, sigma, eps, natoms, volume, pot_lr, vir_lr)
            pot = pot + pot_lr
            DO k = 1, 3
               virial(k,k) = virial(k,k) + vir_lr
            ENDDO

         END SUBROUTINE

      END MODULE

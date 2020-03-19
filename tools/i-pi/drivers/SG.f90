! This performs the calculations necessary to run a simulation using a
! Silvera-Goldman (SG) potential for para-hydrogen. See I. Silvera and V.
! Goldman, J. Chem. Phys., 69, 4209 (1978).
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
! virial tensor of liquid para-hydrogen
! Includes functions which calculate the long-range correction terms for a
! simulation with a sharp nearest-neighbour cut-off.
!
! Functions:
!    f_c: Calculates the damping function for the dispersive interactions 
!       at short range.
!    exp_func: Calculates the short range repulsive part of the SG potential.
!    SG_functions: Calculates the SG pair potential and the magnitude of the
!       forces acting on a pair of atoms.
!    SG_fij: Calculates the SG pair potential and force vector for the
!       interaction of a pair of atoms.
!    SG_longrange: Calculates the long range correction to the potential
!       and virial.
!    SG_getall: Calculates the potential and virial of the system and the
!       forces acting on all the atoms.

      MODULE SG
         USE DISTANCE
      IMPLICIT NONE

      ! Parameters of the SG potential. This potential is of the form:
      ! V(r) = exp(alpha - beta*r - delta*r**2) - 
      !        (C_6/r**6 + C_8/r**8 - C_9/r**9 + C_10/r**10)*f_c(r)
      ! where f_c(r) = exp(-(rc_exp/r - 1)**2) if r <= rc_exp
      !            = 1 otherwise
      DOUBLE PRECISION, PARAMETER :: tau = 6.2831853071795862d0 !If you don't know why I used this name, you're not a real nerd
      DOUBLE PRECISION, PARAMETER :: alpha = 1.713d0
      DOUBLE PRECISION, PARAMETER :: beta = 1.5671d0
      DOUBLE PRECISION, PARAMETER :: delta = 0.00993d0
      DOUBLE PRECISION, PARAMETER :: delta_diff = delta*2.0d0
      DOUBLE PRECISION, PARAMETER :: rc_exp = 8.32d0
      DOUBLE PRECISION, PARAMETER :: C_6 = 12.14d0
      DOUBLE PRECISION, PARAMETER :: C_8 = 215.2d0
      DOUBLE PRECISION, PARAMETER :: C_9 = 143.1d0
      DOUBLE PRECISION, PARAMETER :: C_10 = 4813.9d0
      DOUBLE PRECISION, PARAMETER :: C_6_diff = C_6*6d0
      DOUBLE PRECISION, PARAMETER :: C_8_diff = C_8*8d0
      DOUBLE PRECISION, PARAMETER :: C_9_diff = C_9*9d0
      DOUBLE PRECISION, PARAMETER :: C_10_diff = C_10*10d0
      DOUBLE PRECISION, PARAMETER :: C_6_int = C_6/3d0
      DOUBLE PRECISION, PARAMETER :: C_8_int = C_8/5d0
      DOUBLE PRECISION, PARAMETER :: C_9_int = C_9/6d0
      DOUBLE PRECISION, PARAMETER :: C_10_int = C_10/7d0

      CONTAINS

         SUBROUTINE f_c(r, long_range, long_range_diff)
            ! Calculates the damping function for the dispersive interactions 
            ! at short range. 
            !
            ! Args:
            !    r: The separation of the atoms.
            !    long_range: The value of the damping function.
            !    long_range_diff: The differential of the damping function 
            !       with respect to r.

            DOUBLE PRECISION, INTENT(IN) :: r
            DOUBLE PRECISION, INTENT(OUT) :: long_range
            DOUBLE PRECISION, INTENT(OUT) :: long_range_diff

            DOUBLE PRECISION dist_frac

            IF (r > rc_exp) THEN
               long_range = 1.0d0
               long_range_diff = 0.0d0
            ELSE
               dist_frac = rc_exp/r - 1.0d0
               long_range = dexp(-(dist_frac)**2)
               long_range_diff = 2.0d0*dist_frac*rc_exp*long_range/(r*r)
            END IF

         END SUBROUTINE

         SUBROUTINE exp_func(r, pot, force)
            ! Calculates the repulsive part of the SG force and potential 
            ! between a pair of atoms at a given distance from each other.
            !
            ! Args:
            !    r: The separation of the atoms.
            !    pot: The repulsive part of the potential energy.
            !    force: The magnitude of the repulsive part of the force.

            DOUBLE PRECISION, INTENT(IN) :: r
            DOUBLE PRECISION, INTENT(OUT) :: pot
            DOUBLE PRECISION, INTENT(OUT) :: force

            pot = dexp(alpha - r*(beta + delta*r))
            force = (beta + delta_diff*r)*pot

         END SUBROUTINE

         SUBROUTINE SG_functions(r, pot, force)
            ! Calculates the magnitude of the SG force and potential between
            ! a pair of atoms at a given distance from each other.
            !
            ! Args:
            !    r: The separation of the atoms.
            !    pot: The SG interaction potential.
            !    force: The magnitude of the SG force.

            DOUBLE PRECISION, INTENT(IN) :: r
            DOUBLE PRECISION, INTENT(OUT) :: pot
            DOUBLE PRECISION, INTENT(OUT) :: force

            DOUBLE PRECISION long_range, long_range_diff, disp, disp_diff, exp_pot, exp_force
            DOUBLE PRECISION onr3, onr6, onr9, onr10

            onr3 = 1/(r*r*r)
            onr6 = onr3*onr3
            onr9 = onr6*onr3
            onr10 = onr9/r

            CALL exp_func(r, exp_pot, exp_force)
            CALL f_c(r, long_range, long_range_diff)
      
            disp = -(C_6*onr6 + C_8*onr9*r - C_9*onr9 + C_10*onr10)
            disp_diff = (C_6_diff*onr6/r + C_8_diff*onr9 - C_9_diff*onr10 + C_10_diff*onr10/r)

            pot = exp_pot + disp*long_range
            force = exp_force - disp_diff*long_range - disp*long_range_diff

         END SUBROUTINE

         SUBROUTINE SG_fij(rij, r, pot, fij)
            ! This calculates the SG potential energy and the magnitude and
            ! direction of the force acting on a pair of atoms.
            !
            ! Args:
            !    rij: The vector joining the two atoms.
            !    r: The separation of the two atoms.
            !    pot: The SG interaction potential.
            !    fij: The SG force vector.

            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rij
            DOUBLE PRECISION, INTENT(IN) :: r
            DOUBLE PRECISION, INTENT(OUT) :: pot
            DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: fij
   
            DOUBLE PRECISION f_tot
   
            CALL SG_functions(r, pot, f_tot)
            fij = f_tot*rij/r
   
         END SUBROUTINE

         SUBROUTINE SG_longrange(rc, natoms, volume, pot_lr, vir_lr)
            ! Calculates the long range correction to the total potential and
            ! virial pressure.
            !
            ! Uses the tail correction for a sharp cut-off, with no smoothing
            ! function, as derived in Martyna and Hughes, Journal of Chemical
            ! Physics, 110, 3275, (1999).
            !
            ! Note that we will assume that rc > rc_exp, and that 
            ! exp(alpha - beta*rc - delta*rc**2) << 0, so we can neglect the
            ! contribution of the repulsive potential and the dispersion 
            ! damping function in the long range correction terms.
            !
            ! Args:
            !    rc: The cut-off radius.
            !    natoms: The number of atoms in the system.
            !    volume: The volume of the system box.
            !    pot_lr: The tail correction to the SG interaction potential.
            !    vir_lr: The tail correction to the SG virial pressure.

            DOUBLE PRECISION, INTENT(IN) :: rc
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, INTENT(IN) :: volume
            DOUBLE PRECISION, INTENT(OUT) :: pot_lr
            DOUBLE PRECISION, INTENT(OUT) :: vir_lr

            DOUBLE PRECISION onr3, onr5, onr6, onr7, prefactor

            onr3 = 1/(rc*rc*rc)
            onr6 = onr3*onr3
            onr5 = onr6*rc
            onr7 = onr6/rc
            prefactor = tau*natoms*natoms/volume

            pot_lr = prefactor*(-C_6_int*onr3 - C_8_int*onr5 + C_9_int*onr6 - C_10_int*onr7)
            vir_lr = prefactor*(-C_6*onr3 - C_8*onr5 + C_9*onr6 - C_10*onr7)/3 + pot_lr

         END SUBROUTINE

         SUBROUTINE SG_getall(rc, natoms, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
            ! Calculates the SG potential energy and virial and the forces 
            ! acting on all the atoms.
            !
            ! Args:
            !    rc: The cut-off radius.
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
               ! Only loops over the neighbour list, not all the atoms.
               DO j = start, index_list(i)
                  CALL vector_separation(cell_h, cell_ih, atoms(i,:), atoms(n_list(j),:), rij, r2)
                  IF (r2 < rc*rc) THEN ! Only calculates contributions between neighbouring particles.
                     CALL SG_fij(rij, sqrt(r2), pot_ij, fij)

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
            CALL SG_longrange(rc, natoms, volume, pot_lr, vir_lr)
            pot = pot + pot_lr
            DO k = 1, 3
               virial(k,k) = virial(k,k) + vir_lr
            ENDDO

         END SUBROUTINE

      END MODULE

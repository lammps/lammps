! Create LAMMPS data file for collection of
!   polymer bead-spring chains of various lengths and bead sizes
! Syntax: chain < def.chain > data.file
!   def.chain is input file that specifies the chains
!   data.file is output file that will be input for LAMMPS
! includes image flags in data file so chains can be unraveled later

MODULE boxchain
  IMPLICIT NONE
  PUBLIC
  REAL(KIND=8) :: xprd,yprd,zprd,xboundlo,xboundhi,yboundlo,yboundhi,zboundlo,zboundhi

CONTAINS

  ! periodic boundary conditions - map atom back into periodic box

  SUBROUTINE pbc(x,y,z)
    REAL(KIND=8), INTENT(inout) :: x,y,z

    IF (x < xboundlo) x = x + xprd
    IF (x >= xboundhi) x = x - xprd
    IF (y < yboundlo) y = y + yprd
    IF (y >= yboundhi) y = y - yprd
    IF (z < zboundlo) z = z + zprd
    IF (z >= zboundhi) z = z - zprd

  END SUBROUTINE pbc
END MODULE boxchain

MODULE rngchain
  IMPLICIT NONE

CONTAINS

  ! *very* minimal random number generator

  REAL(KIND=8) FUNCTION random(iseed)
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: iseed
    REAL(KIND=8), PARAMETER :: aa=16807.0_8, mm=2147483647.0_8
    REAL(KIND=8) :: sseed

    sseed = REAL(iseed)
    sseed = MOD(aa*sseed,mm)
    random = sseed/mm
    iseed = INT(sseed)
  END FUNCTION random
END MODULE rngchain

PROGRAM chain
  USE boxchain
  USE rngchain
  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: nchain(:),nmonomer(:)
  INTEGER, ALLOCATABLE :: ntype(:),nbondtype(:)
  INTEGER, ALLOCATABLE :: atomtype(:),molecule(:)
  INTEGER, ALLOCATABLE :: imagex(:),imagey(:),imagez(:)
  REAL(KIND=8), ALLOCATABLE :: x(:),y(:),z(:)
  REAL(KIND=8), ALLOCATABLE :: bondlength(:),restrict(:)
  INTEGER :: i, n, m, nmolecule, natoms, ntypes, nbonds, nbondtypes
  INTEGER :: swaptype, iseed, nsets, iset, ichain, imonomer
  REAL(KIND=8) :: r, rhostar, volume, rsq, xinner, yinner, zinner, xsurf, ysurf, zsurf
  REAL(KIND=8) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, dx, dy, dz

  LOGICAL :: again

  ! read chain definitions

  READ (5,*)
  READ (5,*)
  READ (5,*) rhostar
  READ (5,*) iseed
  READ (5,*) nsets
  READ (5,*) swaptype

  ALLOCATE(nchain(nsets))
  ALLOCATE(nmonomer(nsets))
  ALLOCATE(ntype(nsets))
  ALLOCATE(nbondtype(nsets))
  ALLOCATE(bondlength(nsets))
  ALLOCATE(restrict(nsets))

  DO iset = 1,nsets
      READ (5,*)
      READ (5,*) nchain(iset)
      READ (5,*) nmonomer(iset)
      READ (5,*) ntype(iset)
      READ (5,*) nbondtype(iset)
      READ (5,*) bondlength(iset)
      READ (5,*) restrict(iset)
  ENDDO

  ! natoms = total # of monomers

  natoms = 0
  DO iset = 1,nsets
      natoms = natoms + nchain(iset)*nmonomer(iset)
  ENDDO

  ALLOCATE(x(natoms))
  ALLOCATE(y(natoms))
  ALLOCATE(z(natoms))
  ALLOCATE(atomtype(natoms))
  ALLOCATE(molecule(natoms))
  ALLOCATE(imagex(natoms))
  ALLOCATE(imagey(natoms))
  ALLOCATE(imagez(natoms))

  ! setup box size (sigma = 1.0)

  volume = natoms/rhostar
  xprd = volume**(1.0/3.0)
  yprd = xprd
  zprd = xprd

  xboundlo = -xprd/2.0
  xboundhi = -xboundlo
  yboundlo = xboundlo
  yboundhi = xboundhi
  zboundlo = xboundlo
  zboundhi = xboundhi

  ! generate random chains
  ! loop over sets and chains in each set

  n = 0
  nmolecule = 0

  DO iset = 1,nsets
      DO ichain = 1,nchain(iset)
          nmolecule = nmolecule + 1

          ! random starting point for the chain in the box

          x1 = 0.0
          y1 = 0.0
          z1 = 0.0
          x2 = xboundlo + random(iseed)*xprd
          y2 = yboundlo + random(iseed)*yprd
          z2 = zboundlo + random(iseed)*zprd

          ! store 1st monomer of chain
          ! 1st monomer is always in original box (image = 0)

          CALL pbc(x2,y2,z2)

          n = n + 1
          x(n) = x2
          y(n) = y2
          z(n) = z2
          atomtype(n) = ntype(iset)
          imagex(n) = 0
          imagey(n) = 0
          imagez(n) = 0
          IF (swaptype == 0) THEN
              molecule(n) = nmolecule
          ELSE
              molecule(n) = 1
          END IF

          ! generate rest of monomers in this chain

          DO imonomer = 2, nmonomer(iset)

              x0 = x1
              y0 = y1
              z0 = z1
              x1 = x2
              y1 = y2
              z1 = z2

              again = .TRUE.
              DO WHILE (again)
                  ! random point inside sphere of unit radius

                  xinner = 2.0*random(iseed) - 1.0
                  yinner = 2.0*random(iseed) - 1.0
                  zinner = 2.0*random(iseed) - 1.0
                  rsq = xinner*xinner + yinner*yinner + zinner*zinner
                  IF (rsq > 1.0) CYCLE

                  ! project point to surface of sphere of unit radius

                  r = SQRT(rsq)
                  xsurf = xinner/r
                  ysurf = yinner/r
                  zsurf = zinner/r

                  ! create new point by scaling unit offsets by bondlength (sigma = 1.0)

                  x2 = x1 + xsurf*bondlength(iset)
                  y2 = y1 + ysurf*bondlength(iset)
                  z2 = z1 + zsurf*bondlength(iset)

                  ! check that new point meets restriction requirement
                  ! only for 3rd monomer and beyond

                  dx = x2 - x0
                  dy = y2 - y0
                  dz = z2 - z0
                  r = SQRT(dx*dx + dy*dy + dz*dz)

                  IF (imonomer > 2 .AND. r <= restrict(iset)) CYCLE

                  ! store new point
                  again = .FALSE.

                  ! if delta to previous bead is large, then increment/decrement image flag

                  CALL pbc(x2,y2,z2)
                  n = n + 1
                  x(n) = x2
                  y(n) = y2
                  z(n) = z2
                  atomtype(n) = ntype(iset)

                  IF (ABS(x(n)-x(n-1)) < 2.0*bondlength(iset)) THEN
                      imagex(n) = imagex(n-1)
                  ELSE IF (x(n) - x(n-1) < 0.0) THEN
                      imagex(n) = imagex(n-1) + 1
                  ELSE IF (x(n) - x(n-1) > 0.0) THEN
                      imagex(n) = imagex(n-1) - 1
                  ENDIF

                  IF (ABS(y(n)-y(n-1)) < 2.0*bondlength(iset)) THEN
                      imagey(n) = imagey(n-1)
                  ELSE IF (y(n) - y(n-1) < 0.0) THEN
                      imagey(n) = imagey(n-1) + 1
                  ELSE IF (y(n) - y(n-1) > 0.0) THEN
                      imagey(n) = imagey(n-1) - 1
                  ENDIF

                  IF (ABS(z(n)-z(n-1)) < 2.0*bondlength(iset)) THEN
                      imagez(n) = imagez(n-1)
                  ELSE IF (z(n) - z(n-1) < 0.0) THEN
                      imagez(n) = imagez(n-1) + 1
                  ELSE IF (z(n) - z(n-1) > 0.0) THEN
                      imagez(n) = imagez(n-1) - 1
                  ENDIF

                  IF (swaptype == 0) THEN
                      molecule(n) = nmolecule
                  ELSE IF (swaptype == 1) THEN
                      molecule(n) = imonomer
                  ELSE IF (swaptype == 2) THEN
                      IF (imonomer <= nmonomer(iset)/2) THEN
                          molecule(n) = imonomer
                      ELSE
                          molecule(n) = nmonomer(iset)+1-imonomer
                      ENDIF
                  ENDIF
              ENDDO
          ENDDO

      ENDDO
  ENDDO

  ! compute quantities needed for LAMMPS file

  nbonds = 0
  ntypes = 0
  nbondtypes = 0
  DO iset = 1,nsets
      nbonds = nbonds + nchain(iset)*(nmonomer(iset)-1)
      IF (ntype(iset) > ntypes) ntypes = ntype(iset)
      IF (nbondtype(iset) > nbondtypes) nbondtypes = nbondtype(iset)
  ENDDO

  ! write out LAMMPS file

  WRITE (6,*) 'LAMMPS FENE chain data file'
  WRITE (6,*)

  WRITE (6,*) natoms,' atoms'
  WRITE (6,*) nbonds,' bonds'
  WRITE (6,*) 0,' angles'
  WRITE (6,*) 0,' dihedrals'
  WRITE (6,*) 0,' impropers'
  WRITE (6,*)

  WRITE (6,*) ntypes,' atom types'
  WRITE (6,*) nbondtypes,' bond types'
  WRITE (6,*) 0,' angle types'
  WRITE (6,*) 0,' dihedral types'
  WRITE (6,*) 0,' improper types'
  WRITE (6,*)

  WRITE (6,'(2F15.6,A)') xboundlo,xboundhi,' xlo xhi'
  WRITE (6,'(2F15.6,A)') yboundlo,yboundhi,' ylo yhi'
  WRITE (6,'(2F15.6,A)') zboundlo,zboundhi,' zlo zhi'

  WRITE (6,*)
  WRITE (6,*) 'Masses'
  WRITE (6,*)

  DO i = 1,ntypes
      WRITE (6,'(i3,f5.1)') i,1.0
  ENDDO

  WRITE (6,*)
  WRITE (6,*) 'Atoms # molecular'
  WRITE (6,*)

  DO i = 1,natoms
      WRITE (6,'(I10,I8,I8,3F10.4,3I4)') i,molecule(i),atomtype(i),x(i),y(i),z(i), &
          imagex(i),imagey(i),imagez(i)
  ENDDO

  IF (nbonds > 0) THEN
      WRITE (6,*)
      WRITE (6,*) 'Bonds'
      WRITE (6,*)

      n = 0
      m = 0
      DO iset = 1,nsets
          DO ichain = 1,nchain(iset)
              DO imonomer = 1,nmonomer(iset)
                  n = n + 1
                  IF (imonomer /= nmonomer(iset)) THEN
                      m = m + 1
                      WRITE (6,'(i9,i3,2i9)') m,nbondtype(iset),n,n+1
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
  ENDIF

  DEALLOCATE(nchain, nmonomer, ntype, nbondtype, bondlength, restrict)
  DEALLOCATE(x, y, z, atomtype, molecule, imagex, imagey, imagez)

END PROGRAM chain

! Create LAMMPS data file for 2d LJ simulation of micelles
!
! Syntax: micelle2d < def.micelle2d > data.file
!
! def file contains size of system and number to turn into surfactants
! attaches a random surfactant tail(s) to random head
! if nonflag is set, will attach 2nd-neighbor bonds in surfactant
! solvent atoms = type 1
! micelle heads = type 2
! micelle tails = type 3,4,5,etc.

MODULE boxmicelle
  IMPLICIT NONE
  PUBLIC
  REAL(KIND=8) :: xprd,yprd,zprd,xboundlo,xboundhi,yboundlo,yboundhi,zboundlo,zboundhi

CONTAINS

  ! periodic boundary conditions - map atom back into periodic box

  SUBROUTINE pbc(x,y)
    REAL(KIND=8), INTENT(inout) :: x,y

    IF (x < xboundlo) x = x + xprd
    IF (x >= xboundhi) x = x - xprd
    IF (y < yboundlo) y = y + yprd
    IF (y >= yboundhi) y = y - yprd

  END SUBROUTINE pbc
END MODULE boxmicelle

MODULE rngmicelle
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
END MODULE rngmicelle

PROGRAM micelle2d
  USE boxmicelle
  USE rngmicelle
  IMPLICIT NONE

  REAL(kind=8), ALLOCATABLE :: x(:,:)
  INTEGER, ALLOCATABLE :: atomtype(:), molecule(:)
  INTEGER, ALLOCATABLE :: bondatom(:,:),bondtype(:)
  INTEGER :: natoms, maxatom, ntypes, nbonds, nbondtypes, iseed
  INTEGER :: i, j, k, m, nx, ny, nsurf, ntails, nonflag
  REAL(kind=8) :: rhostar, rlattice, sigma, angle,r0
  REAL(kind=8), parameter :: pi = 3.14159265358979323846_8
  LOGICAL :: again

  READ(5,*)
  READ(5,*)
  READ(5,*) rhostar
  READ(5,*) iseed
  READ(5,*) nx,ny
  READ(5,*) nsurf
  READ(5,*) r0
  READ(5,*) ntails
  READ(5,*) nonflag

  natoms = nx*ny
  maxatom = natoms + nsurf*ntails
  ALLOCATE(x(2,maxatom), molecule(maxatom), atomtype(maxatom))

  nbonds = nsurf*ntails
  IF (nonflag.EQ.1) nbonds = nbonds + nsurf*(ntails-1)
  ALLOCATE(bondatom(2,nbonds), bondtype(nbonds))

! box size

  rlattice = (1.0_8/rhostar) ** 0.5_8

  xboundlo = 0.0_8
  xboundhi = nx*rlattice
  yboundlo = 0.0_8
  yboundhi = ny*rlattice
  zboundlo = -0.1_8
  zboundhi = 0.1_8

  sigma = 1.0_8

  xprd = xboundhi - xboundlo
  yprd = yboundhi - yboundlo
  zprd = zboundhi - zboundlo

! initial square lattice of solvents

  m = 0
  DO j = 1,ny
      DO i = 1,nx
          m = m + 1
          x(1,m) = xboundlo + (i-1)*rlattice
          x(2,m) = yboundlo + (j-1)*rlattice
          molecule(m) = 0
          atomtype(m) = 1
      ENDDO
  ENDDO

! turn some into surfactants with molecule ID
!  head changes to type 2
!  create ntails for each head of types 3,4,...
!  each tail is at distance r0 away in straight line with random orientation

  DO i = 1,nsurf

      again = .TRUE.
      DO WHILE(again)
          m = INT(random(iseed)*natoms + 1)
          IF (m > natoms) m = natoms
          IF (molecule(m) /= 0) CYCLE
          again = .FALSE.
      END DO
      molecule(m) = i
      atomtype(m) = 2

      angle = random(iseed)*2.0_8*pi
      DO j = 1,ntails
          k = (i-1)*ntails + j
          x(1,natoms+k) = x(1,m) + COS(angle)*j*r0*sigma
          x(2,natoms+k) = x(2,m) + SIN(angle)*j*r0*sigma
          molecule(natoms+k) = i
          atomtype(natoms+k) = 2+j
          CALL pbc(x(1,natoms+k),x(2,natoms+k))
          IF (j == 1) bondatom(1,k) = m
          IF (j /= 1) bondatom(1,k) = natoms+k-1
          bondatom(2,k) = natoms+k
          bondtype(k) = 1
      ENDDO

  ENDDO

! if nonflag is set, add (ntails-1) 2nd nearest neighbor bonds to end
!   of bond list
! k = location in bondatom list where nearest neighbor bonds for
!     this surfactant are stored

  IF (nonflag == 1) THEN

      nbonds = nsurf*ntails
      DO i = 1,nsurf
          DO j = 1,ntails-1
              k = (i-1)*ntails + j
              nbonds = nbonds + 1
              bondatom(1,nbonds) = bondatom(1,k)
              bondatom(2,nbonds) = bondatom(2,k+1)
              bondtype(nbonds) = 2
          ENDDO
      ENDDO

  ENDIF

! write LAMMPS data file

  natoms = natoms + nsurf*ntails
  nbonds = nsurf*ntails
  IF (nonflag == 1) nbonds = nbonds + nsurf*(ntails-1)
  ntypes = 2 + ntails
  nbondtypes = 1
  IF (nonflag == 1) nbondtypes = 2

  IF (nsurf == 0) THEN
      ntypes = 1
      nbondtypes = 0
  ENDIF

  WRITE (6,*) 'LAMMPS 2d micelle data file'
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

  WRITE (6,*) xboundlo,xboundhi,' xlo xhi'
  WRITE (6,*) yboundlo,yboundhi,' ylo yhi'
  WRITE (6,*) zboundlo,zboundhi,' zlo zhi'

  WRITE (6,*)
  WRITE (6,*) 'Masses'
  WRITE (6,*)

  DO i = 1,ntypes
      WRITE (6,*) i,1.0
  ENDDO

  WRITE (6,*)
  WRITE (6,*) 'Atoms # molecular'
  WRITE (6,*)

  DO i = 1,natoms
      WRITE (6,'(3I7,3F8.3)') i,molecule(i),atomtype(i),x(1,i),x(2,i),0.0
  ENDDO

  IF (nsurf > 0) THEN

      WRITE (6,*)
      WRITE (6,*) 'Bonds'
      WRITE (6,*)

      DO i = 1,nbonds
          WRITE (6,'(4I7)') i,bondtype(i),bondatom(1,i),bondatom(2,i)
      ENDDO

  ENDIF

  DEALLOCATE(x,molecule,atomtype,bondtype,bondatom)
END PROGRAM micelle2d

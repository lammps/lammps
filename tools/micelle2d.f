c Create LAMMPS 2003 input file for 2d LJ simulation of micelles
c
c Syntax: micelle2d < def.micelle2d > data.file
c
c def file contains size of system and number to turn into surfactants
c attaches a random surfactant tail(s) to random head
c if nonflag is set, will attach 2nd-neighbor bonds in surfactant
c solvent atoms = type 1
c micelle heads = type 2
c micelle tails = type 3,4,5,etc.

      program micelle2d
      parameter (maxatom=100000,maxbond=10000)
      real*4 x(2,maxatom)
      integer type(maxatom),molecule(maxatom)
      integer bondatom(2,maxbond),bondtype(maxbond)
      common xprd,yprd,zprd,xboundlo,xboundhi,yboundlo,yboundhi,
     $     zboundlo,zboundhi
 999  format (3i7,3f8.3)
 998  format (4i7)

      read (5,*)
      read (5,*)
      read (5,*) rhostar
      read (5,*) iseed
      read (5,*) nx,ny
      read (5,*) nsurf
      read (5,*) r0
      read (5,*) ntails
      read (5,*) nonflag
      
      natoms = nx*ny
      if (natoms+nsurf*ntails.gt.maxatom) then
        write (6,*) 'Too many atoms - boost maxatom'
        call exit(0)
      endif

      nbonds = nsurf*ntails
      if (nonflag.eq.1) nbonds = nbonds + nsurf*(ntails-1)
      if (nbonds.gt.maxbond) then
        write (6,*) 'Too many surfactants - boost maxbond'
        call exit(0)
      endif

c box size

      rlattice = (1.0/rhostar) ** 0.5

      xboundlo = 0.0
      xboundhi = nx*rlattice
      yboundlo = 0.0
      yboundhi = ny*rlattice
      zboundlo = -0.1
      zboundhi = 0.1
      
      sigma = 1.0

      xprd = xboundhi - xboundlo
      yprd = yboundhi - yboundlo
      zprd = zboundhi - zboundlo

c initial square lattice of solvents

      m = 0
      do j = 1,ny
        do i = 1,nx
          m = m + 1
          x(1,m) = xboundlo + (i-1)*rlattice
          x(2,m) = yboundlo + (j-1)*rlattice
          molecule(m) = 0
          type(m) = 1
        enddo
      enddo

c turn some into surfactants with molecule ID
c  head changes to type 2
c  create ntails for each head of types 3,4,...
c  each tail is at distance r0 away in straight line with random orientation

      do i = 1,nsurf

 10     m = random(iseed)*natoms + 1
        if (m.gt.natoms) m = natoms
        if (molecule(m) .ne. 0) goto 10

        molecule(m) = i
        type(m) = 2

        angle = random(iseed)*2.0*3.1415926
        do j = 1,ntails
          k = (i-1)*ntails + j 
          x(1,natoms+k) = x(1,m) + cos(angle)*j*r0*sigma
          x(2,natoms+k) = x(2,m) + sin(angle)*j*r0*sigma
          molecule(natoms+k) = i
          type(natoms+k) = 2+j
          call pbc(x(1,natoms+k),x(2,natoms+k))
          if (j.eq.1) bondatom(1,k) = m
          if (j.ne.1) bondatom(1,k) = natoms+k-1
          bondatom(2,k) = natoms+k
          bondtype(k) = 1
        enddo

      enddo

c if nonflag is set, add (ntails-1) 2nd nearest neighbor bonds to end
c   of bond list
c k = location in bondatom list where nearest neighbor bonds for
c     this surfactant are stored

      if (nonflag.eq.1) then

        nbonds = nsurf*ntails
        do i = 1,nsurf
          do j = 1,ntails-1
            k = (i-1)*ntails + j
            nbonds = nbonds + 1
            bondatom(1,nbonds) = bondatom(1,k)
            bondatom(2,nbonds) = bondatom(2,k+1)
            bondtype(nbonds) = 2
          enddo
        enddo

      endif

c write LAMMPS data file

      natoms = natoms + nsurf*ntails
      nbonds = nsurf*ntails
      if (nonflag.eq.1) nbonds = nbonds + nsurf*(ntails-1)
      ntypes = 2 + ntails
      nbondtypes = 1
      if (nonflag.eq.1) nbondtypes = 2

      if (nsurf.eq.0) then
        ntypes = 1
        nbondtypes = 0
      endif

      write (6,*) 'LAMMPS 2d micelle data file'
      write (6,*)

      write (6,*) natoms,' atoms'
      write (6,*) nbonds,' bonds'
      write (6,*) 0,' angles'
      write (6,*) 0,' dihedrals'
      write (6,*) 0,' impropers'
      write (6,*)

      write (6,*) ntypes,' atom types'
      write (6,*) nbondtypes,' bond types'
      write (6,*) 0,' angle types'
      write (6,*) 0,' dihedral types'
      write (6,*) 0,' improper types'
      write (6,*)

      write (6,*) xboundlo,xboundhi,' xlo xhi'
      write (6,*) yboundlo,yboundhi,' ylo yhi'
      write (6,*) zboundlo,zboundhi,' zlo zhi'

      write (6,*)
      write (6,*) 'Masses'
      write (6,*)

      do i = 1,ntypes
        write (6,*) i,1.0
      enddo

      write (6,*)
      write (6,*) 'Atoms'
      write (6,*)

      do i = 1,natoms
        write (6,999) i,molecule(i),type(i),x(1,i),x(2,i),0.0
      enddo

      if (nsurf.gt.0) then

        write (6,*)
        write (6,*) 'Bonds'
        write (6,*)

        do i = 1,nbonds
          write (6,998) i,bondtype(i),bondatom(1,i),bondatom(2,i)
        enddo
      
      endif

c write Xmovie bond geometry file

      open(1,file='bond.micelle')

      write (1,*) 'ITEM: BONDS'
      do i = 1,nbonds
        write (1,*) bondtype(i),bondatom(1,i),bondatom(2,i)
      enddo

      close(1)

      end


c ************
c Subroutines
c ************

c periodic boundary conditions - map atom back into periodic box

      subroutine pbc(x,y)
      common xprd,yprd,zprd,xboundlo,xboundhi,yboundlo,yboundhi,
     $     zboundlo,zboundhi

      if (x.lt.xboundlo) x = x + xprd
      if (x.ge.xboundhi) x = x - xprd
      if (y.lt.yboundlo) y = y + yprd
      if (y.ge.yboundhi) y = y - yprd

      return
      end


c RNG - compute in double precision, return single
      
      real*4 function random(iseed)
      real*8 aa,mm,sseed
      parameter (aa=16807.0D0,mm=2147483647.0D0)
      
      sseed = iseed
      sseed = mod(aa*sseed,mm)
      random = sseed/mm
      iseed = sseed

      return
      end

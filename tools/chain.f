c Create LAMMPS data file for collection of
c   polymer bead-spring chains of various lengths and bead sizes
c Syntax: chain < def.chain > data.file
c   def.chain is input file that specifies the chains
c   data.file is output file that will be input for LAMMPS
c includes image flags in data file so chains can be unraveled later

      program chain

      integer swaptype
      integer, allocatable :: nchain(:),nmonomer(:)
      integer, allocatable :: ntype(:),nbondtype(:)
      integer, allocatable :: type(:),molecule(:)
      integer, allocatable :: imagex(:),imagey(:),imagez(:)
      real*8, allocatable :: x(:),y(:),z(:)
      real*8, allocatable :: bondlength(:),restrict(:)
      common xprd,yprd,zprd,xboundlo,xboundhi,
     $     yboundlo,yboundhi,zboundlo,zboundhi
      real*8 random
 900  format(a)
 901  format(2f15.6,a)
 902  format(i3,f5.1)
 903  format(i10,i8,i8,3f10.4,3i4)
 904  format(i9,i3,2i9)

c read chain definitions

      read (5,*)
      read (5,*)
      read (5,*) rhostar
      read (5,*) iseed
      read (5,*) nsets
      read (5,*) swaptype

      allocate(nchain(nsets))
      allocate(nmonomer(nsets))
      allocate(ntype(nsets))
      allocate(nbondtype(nsets))
      allocate(bondlength(nsets))
      allocate(restrict(nsets))

      do iset = 1,nsets
        read (5,*)
        read (5,*) nchain(iset)
        read (5,*) nmonomer(iset)
        read (5,*) ntype(iset)
        read (5,*) nbondtype(iset)
        read (5,*) bondlength(iset)
        read (5,*) restrict(iset)
      enddo

c natoms = total # of monomers

      natoms = 0
      do iset = 1,nsets
        natoms = natoms + nchain(iset)*nmonomer(iset)
      enddo

      allocate(x(natoms))
      allocate(y(natoms))
      allocate(z(natoms))
      allocate(type(natoms))
      allocate(molecule(natoms))
      allocate(imagex(natoms))
      allocate(imagey(natoms))
      allocate(imagez(natoms))

c setup box size (sigma = 1.0)

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

c generate random chains
c loop over sets and chains in each set
      
      n = 0
      nmolecule = 0

      do iset = 1,nsets
        do ichain = 1,nchain(iset)
          nmolecule = nmolecule + 1

c random starting point for the chain in the box

          x1 = 0.0
          y1 = 0.0
          z1 = 0.0
          x2 = xboundlo + random(iseed)*xprd
          y2 = yboundlo + random(iseed)*yprd
          z2 = zboundlo + random(iseed)*zprd

c store 1st monomer of chain
c 1st monomer is always in original box (image = 0)

          call pbc(x2,y2,z2)
          n = n + 1
          x(n) = x2
          y(n) = y2
          z(n) = z2
          type(n) = ntype(iset)
          imagex(n) = 0
          imagey(n) = 0
          imagez(n) = 0
          if (swaptype == 0) then
            molecule(n) = nmolecule
          else
            molecule(n) = 1
          endif

c generate rest of monomers in this chain

          do imonomer = 2,nmonomer(iset)

            x0 = x1
            y0 = y1
            z0 = z1
            x1 = x2
            y1 = y2
            z1 = z2

c random point inside sphere of unit radius

 10         xinner = 2.0*random(iseed) - 1.0
            yinner = 2.0*random(iseed) - 1.0
            zinner = 2.0*random(iseed) - 1.0
            rsq = xinner*xinner + yinner*yinner + zinner*zinner
            if (rsq > 1.0) goto 10

c project point to surface of sphere of unit radius

            r = sqrt(rsq)
            xsurf = xinner/r
            ysurf = yinner/r
            zsurf = zinner/r

c create new point by scaling unit offsets by bondlength (sigma = 1.0)

            x2 = x1 + xsurf*bondlength(iset)
            y2 = y1 + ysurf*bondlength(iset)
            z2 = z1 + zsurf*bondlength(iset)

c check that new point meets restriction requirement
c only for 3rd monomer and beyond

            dx = x2 - x0
            dy = y2 - y0
            dz = z2 - z0
            r = sqrt(dx*dx + dy*dy + dz*dz)

            if (imonomer > 2 .and. r <= restrict(iset)) goto 10

c store new point
c if delta to previous bead is large, then increment/decrement image flag

            call pbc(x2,y2,z2)
            n = n + 1
            x(n) = x2
            y(n) = y2
            z(n) = z2
            type(n) = ntype(iset)

            if (abs(x(n)-x(n-1)) < 2.0*bondlength(iset)) then
              imagex(n) = imagex(n-1)
            else if (x(n) - x(n-1) < 0.0) then
              imagex(n) = imagex(n-1) + 1
            else if (x(n) - x(n-1) > 0.0) then
              imagex(n) = imagex(n-1) - 1
            endif

            if (abs(y(n)-y(n-1)) < 2.0*bondlength(iset)) then
              imagey(n) = imagey(n-1)
            else if (y(n) - y(n-1) < 0.0) then
              imagey(n) = imagey(n-1) + 1
            else if (y(n) - y(n-1) > 0.0) then
              imagey(n) = imagey(n-1) - 1
            endif

            if (abs(z(n)-z(n-1)) < 2.0*bondlength(iset)) then
              imagez(n) = imagez(n-1)
            else if (z(n) - z(n-1) < 0.0) then
              imagez(n) = imagez(n-1) + 1
            else if (z(n) - z(n-1) > 0.0) then
              imagez(n) = imagez(n-1) - 1
            endif

            if (swaptype == 0) then
              molecule(n) = nmolecule
            else if (swaptype == 1) then
              molecule(n) = imonomer
            else if (swaptype == 2) then
              if (imonomer <= nmonomer(iset)/2) then
                molecule(n) = imonomer
              else
                molecule(n) = nmonomer(iset)+1-imonomer
              endif
            endif

          enddo

        enddo
      enddo

c compute quantities needed for LAMMPS file

      nbonds = 0
      ntypes = 0
      nbondtypes = 0
      do iset = 1,nsets
        nbonds = nbonds + nchain(iset)*(nmonomer(iset)-1)
        if (ntype(iset) > ntypes) ntypes = ntype(iset)
        if (nbondtype(iset) > nbondtypes)
     $       nbondtypes = nbondtype(iset)
      enddo

c write out LAMMPS file

      write (6,900) 'LAMMPS FENE chain data file'
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

      write (6,901) xboundlo,xboundhi,' xlo xhi'
      write (6,901) yboundlo,yboundhi,' ylo yhi'
      write (6,901) zboundlo,zboundhi,' zlo zhi'

      write (6,*)
      write (6,900) 'Masses'
      write (6,*)

      do i = 1,ntypes
        write (6,902) i,1.0
      enddo

      write (6,*)
      write (6,900) 'Atoms'
      write (6,*)

      do i = 1,natoms
        write (6,903) i,molecule(i),type(i),x(i),y(i),z(i),
     $       imagex(i),imagey(i),imagez(i)
      enddo

      if (nbonds > 0) then

        write (6,*)
        write (6,900) 'Bonds'
        write (6,*)

        n = 0
        m = 0
        do iset = 1,nsets
          do ichain = 1,nchain(iset)
            do imonomer = 1,nmonomer(iset)
              n = n + 1
              if (imonomer /= nmonomer(iset)) then
                m = m + 1
                write (6,904) m,nbondtype(iset),n,n+1
              endif
            enddo
          enddo
        enddo

      endif

      end

c ************
c Subroutines
c ************

c periodic boundary conditions - map atom back into periodic box

      subroutine pbc(x,y,z)
      common xprd,yprd,zprd,xboundlo,xboundhi,
     $     yboundlo,yboundhi,zboundlo,zboundhi

      if (x < xboundlo) x = x + xprd
      if (x >= xboundhi) x = x - xprd
      if (y < yboundlo) y = y + yprd
      if (y >= yboundhi) y = y - yprd
      if (z < zboundlo) z = z + zprd
      if (z >= zboundhi) z = z - zprd

      return
      end


c RNG from Numerical Recipes
      
      real*8 function random(iseed)
      real*8 aa,mm,sseed
      parameter (aa=16807.0D0,mm=2147483647.0D0)
      
      sseed = iseed
      sseed = mod(aa*sseed,mm)
      random = sseed/mm
      iseed = sseed

      return
      end

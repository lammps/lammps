c -------------------------------------------------------------------------
c Code converts LAMMPS output to .pdb files
c Overlays coordinates from LAMMPS output onto template pdb
c Assumes atom order is the same between the two
c Converts from atom-based pbc to residue-based pbc
c Also assumes starting config fed to LAMMPS had residue-based pbc
c Paul Crozier, SNL, 2002
c -------------------------------------------------------------------------

      module global

      real*8 xprd,yprd,zprd,box(2,3)

      real*8, allocatable :: x(:,:),q(:),mass(:)
      real*8, allocatable :: comx(:),comy(:),comz(:),totmass(:)

      integer ntimestep,natoms,nframes,iframe,nper
      integer nbonds,nangles,ndihedrals,nimpropers,ntypes
      integer nbondtypes,nangletypes,ndihedtypes,nimprotypes
      integer nconfig,nconfig_skip,nskip,nframes_between_pdbs
      integer nmolecules,nprotein_residues

      integer, allocatable :: mytrue(:,:),type(:),molecule(:)
      integer, allocatable :: nboxx(:),nboxy(:),nboxz(:)

      character*200 data_file_path,config_file_path,pdb_file_path
      character*200 compare_file_path

      character*76, allocatable :: outbeg(:),outend(:)

      end module

c -------------------------------------------------------------------------
c -------------------------------------------------------------------------

      program mkpdb
      use global
      implicit none

      call read_in_mkpdb
      call read_data
      call mkpdb_start

      do iframe = nskip+1, nframes
         call find_config
         call read_config
         write(6,*) 'Frame # ', iframe
         if (mod(iframe,nframes_between_pdbs) == 0) call mk_pdb 
      enddo

      write(6,*) 'Done.'
      stop
      end

c -------------------------------------------------------------------------

      subroutine find_config
      use global
      implicit none

      integer l,m,n,i,j,ntotal

      real*8  buf(8)

      if (mod((iframe-1),nper) == 0) then
        n = (iframe-1)/nper + 1
        write(6,*) 'On config file # ', n
        close(21)
c        l = n/100
c        m = n/10 - l*10
c        n = mod(n,10)
        open(21,file=trim(config_file_path)
     $      //char(48+n),status='old')
        rewind 21

c skip the first frame of each config file

        read(21,*) 
        read(21,*) ntimestep
        read(21,*)
        read(21,*) ntotal
        read(21,*)
        read(21,*) box(1,1),box(2,1)
        read(21,*) box(1,2),box(2,2)
        read(21,*) box(1,3),box(2,3)
        read(21,*)
    
        if (ntotal /= natoms) write(6,*) 'Mismatch # of atoms'

        do i = 1, natoms
           read (21,*) (buf(j),j=1,5)
        enddo

      endif

      return
      end

c -------------------------------------------------------------------------

      logical function match(str1,str2,m)
      implicit none

      character*(*) str1,str2
      integer m

      match = .FALSE.
      m = len(str1) + 1
      if (len(str1).gt.len(str2)) return
      if (str1.eq.str2(1:len(str1))) match = .TRUE.

      return
      end

c -------------------------------------------------------------------------

      subroutine mk_pdb
      use global
      implicit none

      integer i,j,k,l,m,n,o,imolecule,ith_pdb
      real*8 xx,yy,zz,shiftx,shifty,shiftz,proteinmass

      ith_pdb = iframe/nframes_between_pdbs
      j = ith_pdb/1E4
      k = (ith_pdb - j*1E4)/1E3
      l = (ith_pdb - j*1E4 - k*1E3)/1E2
      m = (ith_pdb - j*1E4 - k*1E3 - l*1E2)/1E1
      n = (ith_pdb - j*1E4 - k*1E3 - l*1E2 - m*1E1)

      open(26,file=trim(pdb_file_path)//char(48+j)//char(48+k)//
     1              char(48+l)//char(48+m)//char(48+n)//'.pdb')

c Have to convert from pbc applied on an atomic basis to pbc applied
c on a residue basis.

c Step 1: Recenter system based on protein c.o.m.
 
      shiftx = 0.0
      shifty = 0.0
      shiftz = 0.0
      proteinmass = 0.0

      do i = 1, natoms
         imolecule = molecule(i)
         if (imolecule <= nprotein_residues) then
           shiftx = shiftx + (x(1,i) + mytrue(1,i)*xprd)*mass(type(i))
           shifty = shifty + (x(2,i) + mytrue(2,i)*yprd)*mass(type(i))
           shiftz = shiftz + (x(3,i) + mytrue(3,i)*zprd)*mass(type(i))
           proteinmass = proteinmass + mass(type(i))
         endif
      enddo

      shiftx = shiftx/proteinmass 
      shifty = shifty/proteinmass 
      shiftz = shiftz/proteinmass 

      do i = 1, natoms
         x(1,i) = x(1,i) - shiftx
         x(2,i) = x(2,i) - shifty 
         x(3,i) = x(3,i) - shiftz
      enddo

c Step 2: Find the c.o.m. of each residue --- "molecule"

      do i = 1, nmolecules
         comx(i) = 0.0
         comy(i) = 0.0
         comz(i) = 0.0
         totmass(i) = 0.0
      enddo

      do i = 1, natoms
         imolecule = molecule(i)
         comx(imolecule) = comx(imolecule) + 
     1                     (x(1,i) + mytrue(1,i)*xprd)*mass(type(i))
         comy(imolecule) = comy(imolecule) + 
     1                     (x(2,i) + mytrue(2,i)*yprd)*mass(type(i))
         comz(imolecule) = comz(imolecule) + 
     1                     (x(3,i) + mytrue(3,i)*zprd)*mass(type(i))
         totmass(imolecule) = totmass(imolecule) + mass(type(i))
      enddo

      do i = 1, nmolecules
         comx(i) = comx(i)/totmass(i)
         comy(i) = comy(i)/totmass(i)
         comz(i) = comz(i)/totmass(i)
      enddo

c Step 3: Decide how many boxes must be moved in each direction

      do i = 1, nmolecules
         nboxx(i) = nint(comx(i)/xprd)
         nboxy(i) = nint(comy(i)/yprd)
         nboxz(i) = nint(comz(i)/zprd)
      enddo

c Step 4: Apply moves to atoms. Write pdb file.

      do i = 1, natoms
         imolecule = molecule(i)
         xx = x(1,i) + (mytrue(1,i) - nboxx(imolecule))*xprd
         yy = x(2,i) + (mytrue(2,i) - nboxy(imolecule))*yprd
         zz = x(3,i) + (mytrue(3,i) - nboxz(imolecule))*zprd 
         write(26,100) outbeg(i),xx,yy,zz,outend(i)
      enddo   
  100 format(a30,3f8.3,a22)
      write(26,200) 'END'
  200 format(a3)

      close(26)

      return
      end

c -------------------------------------------------------------------------

      subroutine mkpdb_start
      use global
      implicit none

      integer i
      character*76 pdbline(natoms),str
 
      open(25,file=trim(compare_file_path),status='old')
      rewind 25
      do i = 1, natoms
        read(25,100) pdbline(i)
      enddo
     
  100 format (a)

      do i = 1, natoms
        str = pdbline(i)
        read (str(1:30),100) outbeg(i)
        read (str(55:76),100) outend(i)
      enddo

      return
      end

c -------------------------------------------------------------------------
c input data from config file

      subroutine read_config
      use global
      implicit none

c local variables

      integer i,j,itag,itrue,ntotal

      real*8  buf(8)

      read(21,*) 
      read(21,*) ntimestep
      read(21,*)
      read(21,*) ntotal
      read(21,*)
      read(21,*) box(1,1),box(2,1)
      read(21,*) box(1,2),box(2,2)
      read(21,*) box(1,3),box(2,3)
      read(21,*)
    
      if (ntotal /= natoms) write(6,*) 'Mismatch # of atoms'

      xprd = box(2,1) - box(1,1)
      yprd = box(2,2) - box(1,2)
      zprd = box(2,3) - box(1,3)

      do i = 1, natoms
         read (21,*) (buf(j),j=1,5)
         itag = nint(buf(1))
         type(itag)= nint(buf(2))
         x(1,itag) = buf(3)*xprd + box(1,1)
         x(2,itag) = buf(4)*yprd + box(1,2)
         x(3,itag) = buf(5)*zprd + box(1,3)
         mytrue(1,itag) = 0 
         mytrue(2,itag) = 0
         mytrue(3,itag) = 0
      enddo

      return
      end

c -------------------------------------------------------------------------
c read data from input file

      subroutine read_data
      use global
      implicit none

c local variables

      logical match
      integer i,j,jtmp,m,itag
      real*8  buf(7)
      character*80 str

 900  format (a)

      open(27,file=trim(data_file_path),status='old')
      rewind 27

      read (27,*)
      read (27,*)
      read (27,*) natoms
      read (27,*) nbonds
      read (27,*) nangles
      read (27,*) ndihedrals
      read (27,*) nimpropers
      read (27,*)
      read (27,*) ntypes
      if (nbonds.gt.0) read (27,*) nbondtypes
      if (nangles.gt.0) read (27,*) nangletypes
      if (ndihedrals.gt.0) read (27,*) ndihedtypes
      if (nimpropers.gt.0) read (27,*) nimprotypes
      read (27,*)
      read (27,*)
      read (27,*)
      read (27,*)

      allocate(q(natoms))
      allocate(type(natoms))
      allocate(molecule(natoms))
      allocate(mass(natoms))
      allocate(x(3,natoms))
      allocate(mytrue(3,natoms))
      allocate(outbeg(natoms))
      allocate(outend(natoms))

      do

        read (27,*,end=999,err=999)
        read (27,900,end=999,err=999) str
        read (27,*,end=999,err=999)

        if (match('All Done',str,m)) then

          goto 999

        else if (match('Masses',str,m)) then

          write (6,*) 'Masses ...'

          do i = 1,ntypes
             read (27,*) jtmp,mass(i)
          enddo

        else if (match('Atoms',str,m)) then

          write (6,*) 'Atoms ...'

          do i = 1,natoms
             read (27,*) (buf(j),j=1,7)
             itag = nint(buf(1))
             molecule(itag) = nint(buf(2))
             type(itag) = nint(buf(3))
             q(itag) = buf(4)
          enddo

        else if (match('Bonds',str,m)) then

          do i = 1,nbonds
             read (27,*)
          enddo

        else if (match('Angles',str,m)) then

          do i = 1,nangles
             read (27,*)
          enddo

        else if (match('Impropers',str,m)) then

          do i = 1,nimpropers
             read (27,*)
          enddo

        else if (match('Pair Coeffs',str,m)) then

          write (6,*) 'Pair Coeffs ...'

          do i = 1,ntypes
             read (27,*) 
          enddo

        else if (match('Bond Coeffs',str,m)) then

          do i = 1,nbondtypes
             read (27,*)
          enddo

        else if (match('Angle Coeffs',str,m)) then

          do i = 1,nangletypes
             read (27,*)
          enddo

        else if (match('Dihedral Coeffs',str,m)) then

          do i = 1,ndihedtypes
             read (27,*) 
          enddo

        else if (match('Dihedrals',str,m)) then

          do i = 1,ndihedrals
             read (27,*)  
          enddo

          goto 999

        else
          write (6,*) 'UNKNOWN: ',trim(str)
          write (6,*) 'Unknown identifier in data file'
        endif

      enddo

 999  continue

      close (27)

      nmolecules = molecule(natoms)
      allocate(nboxx(nmolecules))
      allocate(nboxy(nmolecules))
      allocate(nboxz(nmolecules))
      allocate(comx(nmolecules))
      allocate(comy(nmolecules))
      allocate(comz(nmolecules))
      allocate(totmass(nmolecules))

      return
      end

c -------------------------------------------------------------------------

c read data from in_mkpdb file

      subroutine read_in_mkpdb
      use global
      implicit none

 100  format (a)

      open(22,file='in_mkpdb')
      rewind 22

      read (22,*) nconfig
      read (22,*) nper
      read (22,*) nconfig_skip
      read (22,*) nframes_between_pdbs
      read (22,*) nprotein_residues
      read (22,100) data_file_path
      read (22,100) config_file_path
      read (22,100) pdb_file_path
      read (22,100) compare_file_path

      nframes = nconfig*nper
      nskip = nconfig_skip*nper
      iframe = nskip

      close (22)

      return
      end

c -------------------------------------------------------------------------

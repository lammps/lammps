*23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
c----------------------------------------------------------------------
      program lammps to cfg
c----------------------------------------------------------------------
*Programming by Jeff Greathouse and Ara Kooser
*Version 1.0 9/1/04
*Version 1.1 4/10/08 by Axel Kohlmeyer
*Sandia National Labs
*Converts LAMMPS atom dump to .cfg files for AtomEye
*This program is provided as is. Please see the README file.
c----------------------------------------------------------------------
      implicit none
      
      character*12 inhist,snapshot
      character*4 fftype(100),name,ciframe

      integer maxatom,maxtype
      parameter(maxatom=50000,maxtype=99)

      integer i,j,iframe,iframe1,iframe2,jframe
      integer natom,iatom,allatom,itype(maxatom)
      integer atype(maxtype),ntype,mass(maxtype)

      real*8 x(maxatom),y(maxatom),z(maxatom),amass(maxtype)
      real*8 xcell,ycell,zcell,upper,lower


c-------Reads in the user input file-------------------------------
      read(*,*) allatom
      if(allatom.gt.maxatom) then
        write(*,*) 'number of total atoms larger than:',maxatom
        write(*,*) 'change maxatom and recompile'
        STOP
      endif

      read(*,*) ntype
      if(ntype.gt.maxtype) then
        write(*,*) 'number of total types larger than:',maxtype
        write(*,*) 'change maxtype and recompile'
        STOP
      endif

      read(*,*) inhist
      read(*,*) iframe1
      read(*,*) iframe2

      do 1, i=1, ntype
        read(*,*) atype(i)
        read(*,*) amass(i)
        mass(i)=anint(amass(i))  
        read(*,*) fftype(i)
 1    continue

c-------Lammps output file is 9, reads in lmps header--------------         
      name=inhist(1:4)
      open(9,file=inhist,status='old',form='formatted')

      iatom=0
      iframe=0
      jframe=0

c---------------------------------------------------------------------
c----------This begins the frame by frame reading section-------------

 9999 continue
      read(9,*,end=999)
      read(9,*,end=999)
      read(9,*,end=999)
      read(9,*,end=999) natom
      read(9,*,end=999)
      read(9,*,end=999) lower,upper
      xcell=upper-lower
      read(9,*,end=999) lower,upper
      ycell=upper-lower
      read(9,*,end=999) lower,upper
      zcell=upper-lower
      read(9,*,end=999)
         
C clear data array.
      do 400 j=1,allatom
        itype(j)=0
        x(j)=0.0d0
        y(j)=0.0d0
        z(j)=0.0d0
 400  continue

      do 440 j=1,natom
        read(9,*,end=999)iatom,itype(iatom),x(iatom),y(iatom),z(iatom)
 440  continue

      jframe=jframe+1
      if(jframe.lt.iframe1)goto 9999
      iframe=iframe+1
      
      if(iframe.ge.iframe2) goto 999
          

c--------------------------------------------------------------------
c-------This section writes each ts to a seperate .cfg file----------
      ciframe='.cfg'
      write(snapshot,'(i5.5,a4)')iframe,ciframe
      open(unit=iframe+20,file=snapshot,status='new',
     *    form='formatted')

      write((iframe+20),'(a22,i7)')'Number of particles = ',natom
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a16)')'A = 1.0 Angstrom'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),435)'H0(1,1) = ',xcell,' A'
      write((iframe+20),'(a14)')'H0(1,2) = 0 A'
      write((iframe+20),'(a14)')'H0(1,3) = 0 A'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a14)')'H0(2,1) = 0 A'
      write((iframe+20),435)'H0(2,2) = ',ycell,' A'
      write((iframe+20),'(a14)')'H0(2,3) = 0 A'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a14)')'H0(3,1) = 0 A'
      write((iframe+20),'(a14)')'H0(3,2) = 0 A'
      write((iframe+20),435)'H0(3,3) = ',zcell,' A'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a12)')'eta(1,1) = 0'
      write((iframe+20),'(a12)')'eta(1,2) = 0'
      write((iframe+20),'(a12)')'eta(1,3) = 0'
      write((iframe+20),'(a12)')'eta(2,2) = 0'
      write((iframe+20),'(a12)')'eta(2,3) = 0'
      write((iframe+20),'(a12)')'eta(3,3) = 0'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),*)
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
      write((iframe+20),'(a1)')'#'
          
 435  format(a11,f10.5,a2)

      do 460, j=1,allatom
        do 450,i=1,ntype
          if(itype(j).eq.atype(i))
     *        write((iframe+20),445)mass(i),fftype(i),x(j),
     *        y(j),z(j),' 0',' 0',' 0'
c---445 is the format for writing atom data to .cfg file------------
 445      format(i3.3,1x,a2,1x,3f9.6,3a2)

 450    continue
 460  continue

      go to 9999

 999  continue
      close(9)      
      end

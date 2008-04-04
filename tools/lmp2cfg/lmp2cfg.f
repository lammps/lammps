*23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
c----------------------------------------------------------------------
      program lammps to cfg
c----------------------------------------------------------------------
*Programming by Jeff Greathouse and Ara Kooser
*Version 1.0 9/1/04
*Sandia National Labs
*Converts LAMMPS atom dump to .cfg files for AtomEye
*This program is provided as is. Please see the README file.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      
      character*12 inhist,snapshot
      character*4 fftype(100),name,ciframe

      integer natom,iatom,itype(50000),q
      integer atype(99),itycon,ntype,mass(99)

      dimension x(50000),y(50000),z(50000)
      dimension amass(99)


c-------Reads in the user input file-------------------------------

          read(*,*) ntype
          read(*,*) inhist
          read(*,*) iframe1
          read(*,*) iframe2
c          write(*,*) ntype
c          write(*,*) inhist
       do 1, i=1, ntype

          read(*,*) atype(i)
          read(*,*) amass(i)
          mass(i)=anint(amass(i))  
          read(*,*) fftype(i)

       
c          write(*,*) atype(i)
c          write(*,*) amass(i)
c          write(*,*) fftype(i)
 1     continue

c-------Lammps output file is 9, reads in lmps header--------------         

         name=inhist(1:4)
         open(9,file=inhist,status='old',form='formatted')
c         open(2,status='new',form='formatted')


          iatom=0
          iframe=0
          jframe=0

c---------------------------------------------------------------------
c----------This begins the frame by frame reading section-------------

 9999    continue

         read(9,*,end=999)
         read(9,*,end=999)
         read(9,*,end=999)
         read(9,*,end=999) natom
         read(9,*,end=999)
         read(9,*,end=999) xlower,xupper
         read(9,*,end=999) ylower,yupper
         read(9,*,end=999) zlower,zupper
         read(9,*,end=999)
 50      format(2f12.5)
         xcell=xupper-xlower
         ycell=yupper-ylower
         zcell=zupper-zlower
          



 1000    format(1x,i5,1x,i2,3f9.5)
         




 

         do 440 j=1,natom
 420      read(9,*,end=999)iatom,itype(iatom),x(iatom),y(iatom),z(iatom)


*23456789|123456789|123456789|123456789|123456789|123456789|123456789|12 


 

 440      continue

          jframe=jframe+1
          if(jframe.lt.iframe1)goto 9999
          iframe=iframe+1
      
          if(iframe.ge.iframe2) goto 999
          

c--------------------------------------------------------------------
c-------This section writes each ts to a seperate .cfg file----------

c         ciframe=char(iframe)
c         snapshot(iframe)=name//ciframe//'.cfg'
c         write(*,*)ciframe
         
c         write(snapshot,'("Cfgs/",i7.7,".cfg")') iframe
         ciframe='.cfg'
c         write(*,*)ciframe 
         write(snapshot,'(i5.5,a4)')iframe,ciframe
c         write(*,*)snapshot
         open(unit=iframe+20,file=snapshot,status='new',
     *   form='formatted')

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
          
 435      format(a11,f10.5,a2)

          do 460, j=1,natom
          
           do 450,i=1,ntype
c            write(*,*)i,amass(i),fftype(i) 
c             write(*,*)ntype
            if(itype(j).eq.atype(i))
     *        write((iframe+20),445)mass(i),fftype(i),x(j),
     *      y(j),z(j),' 0',' 0',' 0'

c---445 is the format for writing atom data to .cfg file------------
 445     format(i3.3,1x,a2,1x,3f9.6,3a2)

 450      continue

 460      continue

          go to 9999

 999      continue
          close(9)      


          end

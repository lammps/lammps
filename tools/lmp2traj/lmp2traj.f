*23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
*********************************************************************
      program densprof 

*********************************************************************
*** Density profile calculations for orthogonal cells
***
***  modified for more atom types (up to 20) and with generalized 
***  surface map arrays - 1 array 40x50x50 for 20 atom types x2 surfaces
***  A.K. 8-13-02
*** Modified for reading in the x direction Ara Kooser
*** Modified for reading in LAMMPS files Jeff Greathouse and Ara Kooser
*** 7/19/04
*********************************************************************

      implicit real*8 (a-h,o-z)

      LOGICAL       exists

      character*12 inhist,outout,outeang,outhb,
     *    outsurf,outacf,outdat,outang
      character name*4,file*10,text*6
      character*5 mapat(20),rtype
      character*8 colname(40)
      character*40 mapcom(40)
      character dummy*18, dummy1*12, dummy2*100

      integer natom, ntype, nmaps, iconn(12), ndip(1000), mapid(20)
      integer otype, htype, atype(1000)
      integer itype(50000), icount(50), iwato(50000), iwath(50000,2)

      dimension  x(50000),y(50000),z(50000),q(50000)
      dimension  g(50,1000),dipx(1000),dipy(1000),dipz(1000),dip(1000)
      dimension  dipa(1000),dipc(1000),angdis(60,100)
      dimension  smap(40,50,50)
      dimension  amapmin(40),amapmax(40),amamap(40),psmap(40)



      dummy1='123456789|123'

c************************************
      read (*,*) inhist
      read (*,*) iframe1
      read (*,*) iframe2
      read (*,*) acell
      read (*,*) bcell
      read (*,*) ccell
      read (*,*) alpha
      read (*,*) beta
      read (*,*) gamma
      read (*,*) vol
      read (*,*) otype
      read (*,*) htype
      read (*,*) watmiz
      read (*,*) watmaz
      read (*,*) zshift
      read (*,*) nmaps
      do 1 i=1,nmaps
         read (*,*) mapcom(i)
         read (*,*) mapid(i)
         read (*,*) colname(i)
         read (*,*) amapmin(i)
         read (*,*) amapmax(i)
  1   continue
**************************************
c      write(*,*) (mapcom(i),mapid(i),colname(i),amapmin(i),
c     *            amapmax(i),i=1,nmaps)
      name=inhist(1:4)
      outout=name//'-out.txt'
      outsurf=name//'-srf.txt'
      outdat=name//'-dat.bin'
      outang=name//'-ang.txt'
      open(11,file=outout,status='new',form='formatted',err=999)
      open(12,file=outsurf,status='new',form='formatted',err=999)
      open(13,file=outang,status='new',form='formatted',err=999)
      write(11,'(a35,a15)') 'Analysis output file - ', outout
c---- removed     
      write(11,'(a35,a15)') 'Input trajectory (*.asc) - ', inhist
      write(11,'(a35,i7)') 'First frame to analyze - ', iframe1
      write(11,'(a35,i7)') 'Total frames to analyze - ', iframe2 
      write(11,'(a35,i7)') 'Total number of atoms - ', natom
      write(11,'(a35,i7)') 'Total number of atom types - ', ntype
      if(nwat.ne.(iwat-1)) write(*,*) 'nwat.ne.(iwat-1)',nwat,(iwat-1)
      write(11,'(a35,i7)/') 'Total number of H2O molecules - ',nwat
      write(11,'(/a13)') 'Atom types: '
      write(11,1112) (i,'=',atype(i),i=1,ntype)
      write(11,'(/a18)') 'Cell parameters: '
      write(11,1111) 'a=',acell,'b=',bcell,'c=',ccell
      write(11,1111) 'alpha=',alpha,'beta=',beta,'gamma=',gamma
      write(11,1113) 'Vol=',vol
 1111 format(3(3x,a6,f10.6))
 1113 format(3x,a6,f13.6/)
 1112 format(6(3x,i3,a1,a5))

c***********************************************
c   physical constants and conversion factors
      pi=4*atan(1.)
      raddeg=180./pi
      rg=1.98719e-3
      r3=3.e3*rg
      calkj=4.184
      an=0.60223
      barcal=an/calkj*1.e-4
      ted=ccell
c      yfn=1.0/dfloat(natom)
c      yfan=an*yfn
      yted=1./ted
      dipwat=0.41*1.60217733*10./3.33564
c**********************************************
      yrc=1000./ccell
      yrc1=100./ccell
      yra=50./acell
      yrb=50./bcell
      ywz=100./(watmaz-watmiz)
      yang=60./180.
c------ i=1,50 corresponds to 50 possible atom types
      do 202 i=1,50
         icount(i)=0
c            write(*,*)(mapcom(k),colname(k),amapmin(k),
c     *        amapmax(k),mapid(k),k=1,nmaps)
         do 201 j=1,1000
            g(i,j)=0.
            if(i.eq.1) then 
               ndip(j)=0
               dipx(j)=0.
               dipy(j)=0.
               dipz(j)=0.
               dip (j)=0.
               dipa(j)=0.
               dipc(j)=0.
            endif
  201    continue
  202 continue

      do 205 k=1,40
       do 204 i=1,50
         do 203 j=1,50
            smap(k,i,j)=0.
  203    continue
  204  continue
  205 continue

      do 302 i=1,50000
         x(i)=0.
         y(i)=0.
         z(i)=0.
  302 continue

      do 305 k=1,60
       do 304 i=1,100
          angdis(k,i)=0.
  304  continue
  305 continue

c***************************************************************************
      iframe=0
      jframe=0

      open(9,file=inhist,status='old',form='formatted',err=999)
c      open(9,file=inhist,status='old',form='unformatted',err=999)

****************************************************************************
 9999 continue
****************************************************************************        
         iatom=0
         itycon=0
c----Sets counter to zero


c---------------------------------TESTING.F 
         read(9,*,end=999)
         read(9,*,end=999) CurrentTime
         read(9,*,end=999)
         read(9,*,end=999) natom
         read(9,*,end=999)
         read(9,*,end=999)
         read(9,*,end=999)
         read(9,*,end=999)
         read(9,*,end=999)

 400     format(i6,i3,3f9.5)
 410     format(1x,i5,1x,i2,3f9.5)
         

      
         do 440 j=1,natom
c----Comment this for a z mineral
c 420        read(9,*,end=999) iatom,itype(iatom),z(iatom),x(iatom),
c     *      y(iatom)

c---- Uncomment this for a z mineral
 420        read(9,*,end=999) iatom,itype(iatom),x(iatom),y(iatom),
     *      z(iatom)          

c              write(11,400)iatom, itype(iatom), x(iatom),y(iatom),
c     *      z(iatom)

c---- itype(iatom) = current atom being read
c---- atype(i) = list of unique atom types

           do 430 i=1,itycon
              if (itype(iatom).eq.atype(i)) goto 440
 430        continue 

c----    increases itycon by one
           itycon=itycon+1
c----   icount(itycon)=0 prevents summation over all frames 
           icount(itycon)=0
           atype(itycon)=itype(iatom)
           ntype=itycon
          
 440    continue
          
c---------------------------------TESTING.F

     
         jframe=jframe+1
         if(jframe.lt.iframe1) goto 9999
         iframe=iframe+1
         if(mod(iframe,100).eq.0) write(*,*) iframe, jframe

         do 9997 i=1,natom

            
c----       Added to change LAMMPS fractional coords to Cart.      
            xx=x(i) * acell
            yy=y(i) * bcell
            zz=(z(i)+zshift) * ccell
c----       This section rescales atoms into the box

            sx=acell*dnint(xx/acell-0.5)
            xx=xx-sx
            sy=bcell*dnint(yy/bcell-0.5) 
            yy=yy-sy
            sz=ccell*dnint(zz/ccell-0.5) 
            zz=zz-sz

            j=int(zz*yrc)+1
            g(itype(i),j)=g(itype(i),j)+1.
            icount(itype(i))=icount(itype(i))+1

***************************************************************
***     This part is now made universal
***vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            ix=int(xx*yra)+1 
            iy=int(yy*yrb)+1 
            do 777 imap=1,nmaps
              if(itype(i).eq.mapid(imap)) then
                if((zz.ge.amapmin(imap)).and.(zz.le.amapmax(imap))) then
                   smap(imap,ix,iy)=smap(imap,ix,iy)+1.

            

c---------binning occurs here with smap count of atoms at given height
                endif
              endif 
  777        continue

             if(itype(i).eq.otype) then

                xo=x(i)
                yo=y(i)
                zo=z(i)+zshift
                xh1=x(i+1)
                xh2=x(i+2)
                yh1=y(i+1)
                yh2=y(i+2)
                zh1=z(i+1)+zshift
                zh2=z(i+2)+zshift

c                sx=acell*dnint(xo/acell-0.5)
c                xo=xo-sx
c                xh1=xh1-sx
c                xh2=xh2-sx
c                sy=bcell*dnint(yo/bcell-0.5) 
c                yo=xo-sy
c                yh1=yh1-sy
c                yh2=yh2-sy
c                sz=ccell*dnint(zo/ccell-0.5) 
c                zo=zo-sz
c                zh1=zh1-sz
c                zh2=zh2-sz

                j=int(zz*yrc)+1
                dipxj=dipwat*(xh1+xh2-2*xo)
                dipyj=dipwat*(yh1+yh2-2*yo)
                dipzj=dipwat*(zh1+zh2-2*zo)
                dipxyj=sqrt(dipxj**2+dipyj**2)
                dipj=sqrt(dipxyj**2+dipzj**2)
                if(abs(dipj).gt.3.5) write(*,*) 'dip',j,dipj
                if(dipj.ne.0.) then
                   dipcj=dipxyj/dipj

c                   if(dipcj.ge.1.) then 
c                      write(*,*) ' dipcj ',dipcj
c                      dipcj=0.999999
c                   elseif(dipcj.le.-1.) then 
c                      write(*,*) ' dipcj ',dipcj
c                      dipcj=-0.999999
c                   endif

                   adipzj=abs(dipzj)
                   if(adipzj.le.0.000001) then
                      write(*,*) j,' dipzj ',dipzj
                      dipaj=0.
                   else
                      dipaj=acos(dipcj)*dipzj/adipzj
c                      dipaj=acosd(dipcj)*dipzj/adipzj

c                      if(abs(dipaj).ge.89.) then
c                         write(*,*) ' vertical dip mom:', dipaj
c                         write(*,*) xo,yo,zo
c                         write(*,*) xh1,yh1,zh1
c                         write(*,*) xh2,yh2,zh2
c                      endif
                   endif
                endif

                dipx(j)=dipx(j)+dipxj
                dipy(j)=dipy(j)+dipyj
                dipz(j)=dipz(j)+dipzj
                dip (j)=dip (j)+dipj
                dipc(j)=dipc(j)+dipcj
                dipa(j)=dipa(j)+dipaj
                ndip(j)=ndip(j)+1

                jj=int(zz*yrc1)+1
                kk=int((dipaj+90.)*yang)+1
                angdis(kk,jj)=angdis(kk,jj)+1.

             endif 

***^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
***     This part is now made universal
***************************************************************

 9997    continue

c      endif
      if(iframe.ge.iframe2) goto 999
      go to 9999

 999  continue
      close(9)

  667 format(3x,' After analysis of ',i6,' configurations: ')

      write(11,667) iframe
      write(12,667) iframe
      write(13,667) iframe

       fact=vol*yrc/(acell*bcell)
      write(11,*) 'Normalization factor:',fact

      write(11,800) (i,i=1,ntype)

      write(11,'(70x,20i10)') (icount(i),i=1,ntype)

 800  format('       r  ',
     *   '   dipX   ','   dipY   ','   dipZ   ',
     *   '   |dip|  ','   cos_dip','   ang_dip',
     *    27i10)
 
      

      do 206 i=1,1000
         fi=dfloat(i)
         rs1=(fi-0.5)/yrc
         if(ndip(i).ne.0) then
            dipx(i)=dipx(i)/ndip(i)
            dipy(i)=dipy(i)/ndip(i)
            dipz(i)=dipz(i)/ndip(i)
            dip (i)=dip (i)/ndip(i)
            dipc(i)=dipc(i)/ndip(i)
            dipa(i)=dipa(i)/ndip(i)
         endif
c         write(11,444) rs1,(g(j,i)*fact/icount(j),j=1,ntype),
         write(11,444) rs1,dipx(i),dipy(i),dipz(i),
     *                dip(i),dipc(i),dipa(i),
     *                (g(j,i)*fact/iframe/1000.,j=1,ntype)     
  206 continue
  444 format(3x, f7.3,26f10.5)

      close(11)

      do 1206 ii=1,nmaps
         amamap(ii)=0.
 1206 continue 

      do 1205 k=1,nmaps
       do 1204 i=1,50
         do 1203 j=1,50
            if (amamap(k).lt.smap(k,i,j)) amamap(k)=smap(k,i,j)
 1203    continue
 1204  continue
 1205 continue

c        write(12,446) acell,bcell,(amamap(i),i=1,nmaps)

      write(12,336) (mapcom(i),mapid(i),colname(i),i=1,nmaps)
      write(12,335) (colname(ii),ii=1,nmaps)

      do 209 j=1,50
       fj=float(j)
       ra=(fj-0.5)/yra
       do 210 k=1,50
          fk=float(k)
          rb=(fk-0.5)/yrb
          do 211 i=1,nmaps
            if(amamap(i).ne.0.) then
               psmap(i)=smap(i,j,k)/amamap(i)
            else
               psmap(i)=0.
            endif
  211     continue

          write(12,446) ra,rb,(psmap(i),i=1,nmaps)

 210     continue
 209  continue

  446 format(42f8.3)
  335 format(//'   ra      rb      ',40a8)
  336 format(1x,a40,':   at.type - ',i5,';   col.title - ',a8)
         
         write(12,*) (mapcom(i),i=1,nmaps)
         write(12,'(a40)')(mapcom(i),i=1,nmaps)

      close(12)

      write(13,*) '  angle ','   z    ',' angdis '
      do 309 j=1,60
       fj=float(j)
       ang=(fj-0.5)/yang-90.
       do 310 k=1,100
          fk=float(k)
          z1=(fk-0.5)/yrc1

          write(13,446) ang,z1,angdis(j,k)/iframe

 310     continue
 309  continue


      close(13)

      stop
      end

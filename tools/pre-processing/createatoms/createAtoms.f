c X. W. Zhou, xzhou@sandia.gov
      include 'createAtoms.h'
      call latgen()
      call delete()
      call hitvalue()
      call disturb()
      call systemshift()
      call colect()
      if (ilatseed .ne. 0) call randomize()
      call velgen()
      call outputfiles()
      if (float(nntype(1)) .ne. 0.0)
     *   write(6,*)(float(nntype(i))/float(nntype(1)),i=1,ntypes)
      write(6,*)'atoms of different species ',(nntype(i),i=1,ntypes)
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine latgen()
      include 'createAtoms.h'
      character*8 lattype
      common /latinfo/ lattype,alat(3),xrot(3),yrot(3),zrot(3),
     *   periodicity(3),xbound(2),ybound(2),zbound(2),strain(3),
     *   delx(3),rcell(3),ccell(nelmax)
      common /iran/ iseed
      namelist /maincard/ ntypes,amass,ielement,ilatseed,
     *   perlb,perub,iseed,xy,xz,yz
      namelist /latcard/ lattype,alat,xrot,yrot,zrot,
     *   periodicity,xbound,ybound,zbound,
     *   strain,delx
      iseed=3245566
      natoms=0
      ilatseed=0
      xy=1.0e12
      xz=1.0e12
      yz=1.0e12
      read(5,maincard) 
      do i = 1,3
         perlen(i)=perub(i)-perlb(i)
      enddo
      do i=1,ntypes
         nntype(i)=0
      enddo
1     continue
      strain(1)=0.0
      strain(2)=0.0
      strain(3)=0.0
      xbound(1)=0.0
      xbound(2)=0.0
      ybound(1)=0.0
      ybound(2)=0.0
      zbound(1)=0.0
      zbound(2)=0.0
      delx(1)=0.0
      delx(2)=0.0
      delx(3)=0.0
      lattype='none'
      read(5,latcard)
      if (lattype .eq. 'none') goto 2
      call latgen1()
      call defvalue()
      goto 1
2     continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine latgen1()
      include 'createAtoms.h'
      character*8 lattype
      common /latinfo/ lattype,alat(3),xrot(3),yrot(3),zrot(3),
     *   periodicity(3),xbound(2),ybound(2),zbound(2),strain(3),
     *   delx(3),rcell(3),ccell(nelmax)
      common /tmpvar/ avec(3),bvec(3),cvec(3),xold(3),yold(3),zold(3),
     *   avecp(3),bvecp(3),cvecp(3),rcellp(3),roter(3,3),rvs(3,100)
      data xold/1.0,0.0,0.0/,yold/0.0,1.0,0.0/,zold/0.0,0.0,1.0/
      namelist /subcard/ rcell,ccell

      call storelat(lattype,avec,bvec,cvec)
      if (xbound(1) .eq. xbound(2)) then
         xbound(1)=perlb(1)
         xbound(2)=perub(1)
      endif
      if (ybound(1) .eq. ybound(2)) then
         ybound(1)=perlb(2)
         ybound(2)=perub(2)
      endif
      if (zbound(1) .eq. zbound(2)) then
         zbound(1)=perlb(3)
         zbound(2)=perub(3)
      endif
      vol=alat(1)*alat(2)*alat(3)
      xnorm=sqrt((xrot(1)*alat(2)*alat(3))**2+
     *           (xrot(2)*alat(1)*alat(3))**2+
     *           (xrot(3)*alat(1)*alat(2))**2)
      ynorm=sqrt((yrot(1)*alat(2)*alat(3))**2+
     *           (yrot(2)*alat(1)*alat(3))**2+
     *           (yrot(3)*alat(1)*alat(2))**2)
      znorm=sqrt((zrot(1)*alat(2)*alat(3))**2+
     *           (zrot(2)*alat(1)*alat(3))**2+
     *           (zrot(3)*alat(1)*alat(2))**2)
      periodicity(1)=vol*periodicity(1)/xnorm
      periodicity(2)=vol*periodicity(2)/ynorm
      periodicity(3)=vol*periodicity(3)/znorm
      do 5 i=1,3
         xrot(i)=vol*xrot(i)/(alat(i)*xnorm)
         yrot(i)=vol*yrot(i)/(alat(i)*ynorm)
         zrot(i)=vol*zrot(i)/(alat(i)*znorm)
5     continue
      do 10 i=1,3
      do 10 j=1,3
         roter(i,j)=xold(i)*xrot(j)+yold(i)*yrot(j)+zold(i)*zrot(j)
10    continue
      do 20 i=1,3
         avecp(i)=0.0
         bvecp(i)=0.0
         cvecp(i)=0.0
      do 20 j=1,3
         avecp(i)=avecp(i)+roter(i,j)*avec(j)*0.5*alat(j)
         bvecp(i)=bvecp(i)+roter(i,j)*bvec(j)*0.5*alat(j)
         cvecp(i)=cvecp(i)+roter(i,j)*cvec(j)*0.5*alat(j)
20    continue
      do 30 i=1,3
         avec(i)=avecp(i)
         bvec(i)=bvecp(i)
         cvec(i)=cvecp(i)
30    continue
      call applystrain()
      nrange=20
      nsmall=0
      do 300 ii=-nrange,nrange
      do 300 jj=-nrange,nrange
      do 300 kk=-nrange,nrange
         xt=ii*avec(1)+jj*bvec(1)+kk*cvec(1)
         if (xt .ge. periodicity(1)-0.1 .or. xt .lt. -0.1) goto 300
         yt=ii*avec(2)+jj*bvec(2)+kk*cvec(2)
         if (yt .ge. periodicity(2)-0.1 .or. yt .lt. -0.1) goto 300
         zt=ii*avec(3)+jj*bvec(3)+kk*cvec(3)
         if (zt .ge. periodicity(3)-0.1 .or. zt .lt. -0.1) goto 300
         nsmall=nsmall+1
         rvs(1,nsmall)=xt
         rvs(2,nsmall)=yt
         rvs(3,nsmall)=zt
300   continue
      nxa=xbound(1)/periodicity(1)-2
      nxb=xbound(2)/periodicity(1)+2
      nya=ybound(1)/periodicity(2)-2
      nyb=ybound(2)/periodicity(2)+2
      nza=zbound(1)/periodicity(3)-2
      nzb=zbound(2)/periodicity(3)+2
1     rcell(1)=-100.0
      read(5,subcard)
      if (rcell(1) .eq. -100.0) then
         return
      endif
      do 190 m = 2,ntypes
         ccell(m)=min(1.0,ccell(m-1)+ccell(m))
190   continue
      do 800 i=1,3
         rcellp(i)=0.0
      do 800 j=1,3
         rcellp(i)=rcellp(i)+roter(i,j)*rcell(j)*alat(j)
800   continue
      do 801 i=1,3
         rcell(i)=rcellp(i)*(1.0+strain(i))
801   continue
      do 250 nn=1,nsmall
      do 251 ii=nxa,nxb
         xt=ii*periodicity(1)+rvs(1,nn)+rcell(1)
         if (xt .ge. xbound(2) .or. xt .lt. xbound(1)) goto 251
      do 252 jj=nya,nyb
         yt=jj*periodicity(2)+rvs(2,nn)+rcell(2)
         if (yt .ge. ybound(2) .or. yt .lt. ybound(1)) goto 252
      do 253 kk=nza,nzb
         zt=kk*periodicity(3)+rvs(3,nn)+rcell(3)
         if (zt .ge. zbound(2) .or. zt .lt. zbound(1)) goto 253
         tst=ranl()
         it=0
      do 201 m=ntypes,1,-1
         if (tst .le. ccell(m)) it=m
201   continue
         natoms=natoms+1
         if (natoms .gt. natmax) then
            write(6,9601)natmax
9601        format('number of atoms exceeds maximum dimension: ',i10)
         endif
         itype(natoms)=it
         nhittag(natoms)=0
         nntype(it)=nntype(it)+1
         rv(1,natoms)=xt+delx(1)
         rv(2,natoms)=yt+delx(2)
         rv(3,natoms)=zt+delx(3)
         if (rv(1,natoms) .lt. perlb(1))
     *      rv(1,natoms)=rv(1,natoms)+perlen(1)
         if (rv(1,natoms) .ge. perub(1))
     *      rv(1,natoms)=rv(1,natoms)-perlen(1)
         if (rv(2,natoms) .lt. perlb(2))
     *      rv(2,natoms)=rv(2,natoms)+perlen(2)
         if (rv(2,natoms) .ge. perub(2))
     *      rv(2,natoms)=rv(2,natoms)-perlen(2)
         if (rv(3,natoms) .lt. perlb(3))
     *      rv(3,natoms)=rv(3,natoms)+perlen(3)
         if (rv(3,natoms) .ge. perub(3))
     *      rv(3,natoms)=rv(3,natoms)-perlen(3)
253   continue
252   continue
251   continue
250   continue
      goto 1
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defvalue()
      include 'createAtoms.h'
      integer oldtype
      character*8 lattype
      common /latinfo/ lattype,alat(3),xrot(3),yrot(3),zrot(3),
     *   periodicity(3),xbound(2),ybound(2),zbound(2),strain(3),
     *   delx(3),rcell(3),ccell(nelmax)
      common /tmpvar/ avec(3),bvec(3),cvec(3),xold(3),yold(3),zold(3),
     *   avecp(3),bvecp(3),cvecp(3),rcellp(3),roter(3,3),rvs(3,100)
      dimension cent(3),plane(3),planec(3),surface(3)
      namelist /defcard/ xmin,xmax,ymin,ymax,zmin,zmax,
     *   oldtype,newtype,prob,cent,plane,dist,surface,
     *   nramp,ndirec,dismax,theta,phi,nscrew,anglemin,anglemax,
     *   xbottom,xheight,zbottom,zheight,nrepeatsx,nrepeatsz
      data pi/3.1415927/

1     continue
      xmin=-1.0e-20
      xmax=1.0e20
      ymin=-1.0e20
      ymax=1.0e20
      zmin=-1.0e20
      zmax=1.0e20
      oldtype=0
      dist=-1.0
      newtype=-11
      cent(1)=0.5*(perlb(1)+perub(1))
      cent(2)=0.5*(perlb(2)+perub(2))
      cent(3)=0.5*(perlb(3)+perub(3))
      surface(1)=0.0
      surface(2)=0.0
      surface(3)=0.0
      nramp=1
      nscrew=0
      anglemin=0.0
      anglemax=360.0
      ndirec=3
      dismax=0.0
      theta=5000.0
      phi=5000.0
      nrepeatsx=0
      nrepeatsz=0
      read(5,defcard)
      if (newtype .lt. -10 .and. dismax .eq. 0.0) return
      if (newtype .gt. ntypes) then
         ntypes=ntypes+1
         if (oldtype .ne. 0) then
            amass(ntypes)=amass(oldtype)
         else
            amass(ntypes)=amass(ntypes-1)
         endif 
      endif
      if (surface(1) .eq. 1.0) then
         perlb(1)=perlb(1)-10.0
         perub(1)=perub(1)+10.0
         perlen(1)=perub(1)-perlb(1)
      endif
      if (surface(2) .eq. 1.0) then
         perlb(2)=perlb(2)-10.0
         perub(2)=perub(2)+10.0
         perlen(2)=perub(2)-perlb(2)
      endif
      if (surface(3) .eq. 1.0) then
         perlb(3)=perlb(3)-10.0
         perub(3)=perub(3)+10.0
         perlen(3)=perub(3)-perlb(3)
      endif
      vol=alat(1)*alat(2)*alat(3)
      if (dismax .ne. 0.0) then
         astart=anglemin
         afinish=anglemax
         if (anglemax .lt. anglemin) then
            anglemin=afinish
            anglemax=astart
         endif
         if (nramp .eq. 1) then
            pstart=xmin
            pfinish=xmax
            if (xmax .lt. xmin) then
               xmin=pfinish
               xmax=pstart
            endif
         endif
         if (nramp .eq. 2) then
            pstart=ymin
            pfinish=ymax
            if (ymax .lt. ymin) then
               ymin=pfinish
               ymax=pstart
            endif
         endif
         if (nramp .eq. 3) then
            pstart=zmin
            pfinish=zmax
            if (zmax .lt. zmin) then
               zmin=pfinish
               zmax=pstart
            endif
         endif
         do 33 i=1,natoms
            if (rv(1,i) .lt. xmin) goto 33
            if (rv(1,i) .gt. xmax) goto 33
            if (rv(2,i) .lt. ymin) goto 33
            if (rv(2,i) .gt. ymax) goto 33
            if (rv(3,i) .lt. zmin) goto 33
            if (rv(3,i) .gt. zmax) goto 33
            if (nscrew .eq. 0) then
               if (nramp .eq. 0) then
                  rv(ndirec,i)=rv(ndirec,i)+dismax
               else
                  rv(ndirec,i)=rv(ndirec,i)+
     *               (rv(nramp,i)-pstart)/(pfinish-pstart)*dismax
               endif
            else
               if (ndirec .eq. 1) then
                  dx=rv(2,i)-cent(2)
                  dy=rv(3,i)-cent(3)
                  dxmax=abs(ymax-cent(2))
                  dxmin=abs(ymin-cent(2))
                  rmax=min(dxmin,dxmax)
               endif
               if (ndirec .eq. 2) then
                  dx=rv(3,i)-cent(3)
                  dy=rv(1,i)-cent(1)
                  dxmax=abs(zmax-cent(3))
                  dxmin=abs(zmin-cent(3))
                  rmax=min(dxmin,dxmax)
               endif
               if (ndirec .eq. 3) then
                  dx=rv(1,i)-cent(1)
                  dy=rv(2,i)-cent(2)
                  dxmax=abs(xmax-cent(1))
                  dxmin=abs(xmin-cent(1))
                  rmax=min(dxmin,dxmax)
               endif
               dr=sqrt(dx**2+dy**2)
               if (dr .eq. 0.0) then
                  angle=0.0
               else
                  angle=acos(dx/dr)*180.0/pi
                  if (dy .lt. 0.0) angle=360.0-angle
               endif
               if (angle .lt. anglemin) goto 33
               if (angle .gt. anglemax) goto 33
               if (nscrew .eq. 1) then
                  astart1=astart
                  afinish1=afinish
               else if (nscrew .eq. 2) then
                  if (astart .eq. 0.0 .and.
     *               afinish .eq. 180.0) then
                     if (dr .le. dxmax) then
                        astart1=astart
                     else
                        astart1=acos(dxmax/dr)*180.0/pi
                     endif
                     if (dr .le. dxmin) then
                        afinish1=afinish
                     else
                        afinish1=180.0-acos(dxmin/dr)*180.0/pi
                     endif
                  else if (astart .eq. 180.0 .and.
     *               afinish .eq. 0.0) then
                     if (dr .le. dxmin) then
                        astart1=astart
                     else
                        astart1=180.0-acos(dxmin/dr)*180.0/pi
                     endif
                     if (dr .le. dxmax) then
                        afinish1=afinish
                     else
                        afinish1=acos(dxmax/dr)*180.0/pi
                     endif
                  else if (astart .eq. 180.0 .and.
     *               afinish .eq. 360.0) then
                     if (dr .le. dxmin) then
                        astart1=astart
                     else
                        astart1=180.0+acos(dxmin/dr)*180.0/pi
                     endif
                     if (dr .le. dxmax) then
                        afinish1=afinish
                     else
                        afinish1=360.0-acos(dxmax/dr)*180.0/pi
                     endif
                  else if (astart .eq. 360.0 .and.
     *               afinish .eq. 180.0) then
                     if (dr .le. dxmax) then
                        astart1=astart
                     else
                        astart1=360.0-acos(dxmax/dr)*180.0/pi
                     endif
                     if (dr .le. dxmin) then
                        afinish1=afinish
                     else
                        afinish1=180.0+acos(dxmin/dr)*180.0/pi
                     endif
                  else
                     write(6,*)'problem in nscrew=2'
                     stop
                  endif
               endif
               if (dr .gt. rmax) dr=rmax
               rv(ndirec,i)=rv(ndirec,i)+dr/rmax*
     *            (angle-astart1)/(afinish1-astart1)*dismax
            endif
33       continue
      else if (theta .eq. 5000.0 .and. phi .eq. 5000.0 .and.
     *   nrepeatsx*nrepeatsz .eq. 0.0) then
         if (dist .lt. 0.0) then
            do 2 i=1,natoms
               if (rv(1,i) .lt. xmin) goto 2
               if (rv(1,i) .gt. xmax) goto 2
               if (rv(2,i) .lt. ymin) goto 2
               if (rv(2,i) .gt. ymax) goto 2
               if (rv(3,i) .lt. zmin) goto 2
               if (rv(3,i) .gt. zmax) goto 2
               if (oldtype .eq. itype(i) .or. oldtype .eq. 0) then
                  if (ranl() .lt. prob) itype(i)=newtype
               endif
2           continue
         else
            if (theta .eq. 5000.0) then
               pnorm=sqrt((plane(1)*alat(2)*alat(3))**2+
     *                    (plane(2)*alat(1)*alat(3))**2+
     *                    (plane(3)*alat(1)*alat(2))**2)
               do 30 n=1,3
                  plane(n)=vol*plane(n)/(alat(n)*pnorm)
30             continue
               do 20 n=1,3
                  planec(n)=0.0
               do 20 m=1,3
                  planec(n)=planec(n)+roter(n,m)*plane(m)
20             continue
            else
               planec(1)=sin(pi/180.0*theta)*cos(pi/180.0*phi)
               planec(3)=sin(pi/180.0*theta)*sin(pi/180.0*phi)
               planec(2)=cos(pi/180.0*theta)
            endif
            do 22 i=1,natoms
               dis1=rv(1,i)-cent(1)
               dis2=rv(2,i)-cent(2)
               dis3=rv(3,i)-cent(3)
               if (dis1*planec(1)+dis2*planec(2)+dis3*planec(3)
     *         .gt. dist) then
                  if (oldtype .eq. itype(i) .or. oldtype .eq. 0)
     *               itype(i)=newtype
               endif
22          continue
         endif
      else if (nrepeatsx*nrepeatsz .ne. 0) then
         do 4 i=1,natoms
            if (rv(1,i) .le. xmin .or. rv(1,i) .ge. xmax) then
               ylocalx = xbottom
            else
               ylocalx = xbottom + 0.5*xheight*(1.0+
     *            cos(2.0*pi*nrepeatsx/(xmax-xmin)*
     *            (rv(1,i)-(xmin+xmax)/(2.0*nrepeatsx))))
            endif
            if (rv(3,i) .le. zmin .or. rv(3,i) .ge. zmax) then
               ylocalz = zbottom
            else
               ylocalz = zbottom + 0.5*zheight*(1.0+
     *            cos(2.0*pi*nrepeatsz/(zmax-zmin)*
     *            (rv(3,i)-(zmin+zmax)/(2.0*nrepeatsz))))
            endif
            ylocal = ylocalx
            if (ylocal .gt. ylocalz) ylocal = ylocalz
            if (rv(2,i) .gt. ylocal) then
               itype(i)=newtype
            endif
4        continue
      endif

      goto 1
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine storelat(latty,avec,bvec,cvec)
      character*8 lstored(10),latty
      dimension avecs(3,10),bvecs(3,10),cvecs(3,10)
      dimension avec(3),bvec(3),cvec(3)
      data lstored/'fcc     ','bcc     ','sc      ',7*'        '/
      data avecs/1.0,1.0,0.0, 1.0,1.0,-1.0, 2.0,0.0,0.0, 21*0.0/
      data bvecs/0.0,1.0,1.0, -1.0,1.0,1.0, 0.0,2.0,0.0, 21*0.0/
      data cvecs/1.0,0.0,1.0, 1.0,-1.0,1.0, 0.0,0.0,2.0, 21*0.0/
      imatch=0
      do 10 i=1,10
         if (latty .ne. lstored(i)) goto 10
         imatch = 1
         do 5 j=1,3
            avec(j)=avecs(j,i)
            bvec(j)=bvecs(j,i)
            cvec(j)=cvecs(j,i)
 5       continue
10    continue
      if (imatch .eq. 1) return
      write(6,15) latty
15    format('could not find lattice type ')
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine applystrain()
      include 'createAtoms.h'
      character*8 lattype
      common /latinfo/ lattype,alat(3),xrot(3),yrot(3),zrot(3),
     *   periodicity(3),xbound(2),ybound(2),zbound(2),strain(3),
     *   delx(3),rcell(3),ccell(nelmax)
      common /tmpvar/ avec(3),bvec(3),cvec(3),xold(3),yold(3),zold(3),
     *   avecp(3),bvecp(3),cvecp(3),rcellp(3),roter(3,3),rvs(3,100)
      do i=1,3
         periodicity(i)=(1.0+strain(i))*periodicity(i)
         avec(i)=(1.0+strain(i))*avec(i)
         bvec(i)=(1.0+strain(i))*bvec(i)
         cvec(i)=(1.0+strain(i))*cvec(i)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hitvalue()
      include 'createAtoms.h'
      namelist /hitcard/ xmin,xmax,ymin,ymax,zmin,zmax,nhit
      nhitcards=0
1001  nhit=0
      xmin=-10000.0
      xmax=10000.0
      ymin=-10000.0
      ymax=10000.0
      zmin=-10000.0
      zmax=10000.0
      read(5,hitcard)
      if (nhit .eq. 0) then
         return
      endif
      nchange=0
      do 1002 i=1,natoms
         if (rv(1,i) .lt. xmin) goto 1002
         if (rv(1,i) .gt. xmax) goto 1002
         if (rv(2,i) .lt. ymin) goto 1002
         if (rv(2,i) .gt. ymax) goto 1002
         if (rv(3,i) .lt. zmin) goto 1002
         if (rv(3,i) .gt. zmax) goto 1002
         nchange=nchange+1
         nhittag(i)=nhit
1002  continue
      if (nchange .ne. 0) nhitcards=nhitcards+1
      goto 1001
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine systemshift()
      include 'createAtoms.h'
      namelist /shiftcard/ mode
      mode=0
      read(5,shiftcard)
      if (mode .eq. 0) return
      if (mode .eq. 1) then
         top=perub(2)
         bottom=perlb(2)
         sdely=0.0
         do i=1,natoms
            if (nhittag(i) .eq. 3 .and. top .gt. rv(2,i)) top=rv(2,i)
            if (nhittag(i) .eq. 4 .and. bottom .lt. rv(2,i))
     *         bottom=rv(2,i)
         enddo
         sdely=-0.5*(top+bottom)
         do i=1,natoms
            rv(2,i)=rv(2,i)+sdely
         enddo     
      endif
      if (mode .eq. 2) then
c     This shift only assumes a maximum of two plane spacings.
         xmin=10.0
         do i=1,natoms
            if (rv(1,i) .lt. xmin) xmin=rv(1,i)
         enddo
         shiftneg=-10.0
         shiftpos=10.0
         do i=1,natoms
            shift=rv(1,i)-xmin
            shift=shift-perlen(1)*nint(shift/perlen(1))
            if (shift .gt. 0.01 .and. shift .lt. shiftpos)
     *         shiftpos=shift
            if (shift .lt. -0.01 .and. shift .gt. shiftneg)
     *         shiftneg=shift
         enddo
         if (-shiftneg .gt. shiftpos) then
            shift=0.5*shiftneg
         else
            shift=0.5*shiftpos
         endif
         perlb(1)=xmin+shift
         perub(1)=perlb(1)+perlen(1)
      endif 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine colect()
      include 'createAtoms.h'
      do i=1,natoms
         if (rv(1,i) .lt. perlb(1)) rv(1,i)=rv(1,i)+perlen(1)
         if (rv(1,i) .ge. perub(1)) rv(1,i)=rv(1,i)-perlen(1)
         if (rv(2,i) .lt. perlb(2)) rv(2,i)=rv(2,i)+perlen(2)
         if (rv(2,i) .ge. perub(2)) rv(2,i)=rv(2,i)-perlen(2)
         if (rv(3,i) .lt. perlb(3)) rv(3,i)=rv(3,i)+perlen(3)
         if (rv(3,i) .ge. perub(3)) rv(3,i)=rv(3,i)-perlen(3)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function ranl()
      common /iran/ iseed
      iseed=iseed*125
      iseed=iseed-2796203*(iseed/2796203)
      ranl=abs(iseed/2796203.0)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine delete()
      include 'createAtoms.h'
      nw_del=0
      do i=1,ntypes
         nntype(i)=0
      enddo
      natoms0=natoms
      natoms=0
      do 1 i=1,natoms0
         if (itype(i) .eq. -1) nw_del=nw_del+1
         if (itype(i) .le. 0) goto 1
         natoms=natoms+1
         ntag(natoms)=natoms
         itype(natoms)=itype(i)
         nhittag(natoms)=nhittag(i)
         nntype(itype(natoms))=nntype(itype(natoms))+1
         rv(1,natoms)=rv(1,i)
         rv(2,natoms)=rv(2,i)
         rv(3,natoms)=rv(3,i)
1     continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine randomize()
      include 'createAtoms.h'
      dimension rv0(6,natmax),itype0(natmax),nhittag0(natmax)
      do 1 i=1,natoms
         itype0(i)=itype(i)
         if (nhitcards .gt. 0) nhittag0(i)=nhittag(i)
         rv0(1,i)=rv(1,i)
         rv0(2,i)=rv(2,i)
         rv0(3,i)=rv(3,i)
         rv0(4,i)=rv(4,i)
         rv0(5,i)=rv(5,i)
         rv0(6,i)=rv(6,i)
1     continue
      im=1
10    im=im*2
      if (im .lt. natoms) goto 10
      ia=365
      ic=8161
      ivalue=mod(ilatseed,natoms)
      do 2 i=1,natoms
20       ivalue=mod(ia*ivalue+ic,im)
         if (ivalue .ge. natoms) goto 20
         ntmp=ivalue+1
         if (ntmp .ge. 1 .and. ntmp .le. natoms) then
            ntag(ntmp)=i
            itype(ntmp)=itype0(i)
            if (nhitcards .gt. 0) nhittag(ntmp)=nhittag0(i)
            rv(1,ntmp)=rv0(1,i)
            rv(2,ntmp)=rv0(2,i)
            rv(3,ntmp)=rv0(3,i)
            rv(4,ntmp)=rv0(4,i)
            rv(5,ntmp)=rv0(5,i)
            rv(6,ntmp)=rv0(6,i)
         endif
2     continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write0(dumpfile)
      include 'createAtoms.h'
      character*32 dumpfile, typ
      open (unit=2,file=dumpfile)
      write(2,*) natoms
      write(2,100) 'GaN lattice'
      do i=1,natoms
         if (itype(i).eq.1) then
            typ = 'Ga'
         else
            typ = 'N'
         endif
         write(2,150) typ(1:length(typ)),rv(1,i),rv(2,i),rv(3,i)
      enddo
 50   format(i8)
 100  format(a)
 150  format(a3,3f13.5)
      close(2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write1(dumpfile)
      include 'createAtoms.h'
      character*32 dumpfile
      zero=0.0
      open (unit=2,file=dumpfile)
      write(2,*) 'ITEM: NUMBER OF ATOMS'
      write(2,*) natoms
      write(2,*) 'ITEM: BOX BOUNDS'
      write(2,*) perlb(1),perub(1)
      write(2,*) perlb(2),perub(2)
      write(2,*) perlb(3),perub(3)
      write(2,*) 'ITEM: TIME'
      write(2,*) 0.0
      write(2,*) 'ITEM: ATOMS'
      if (nhitcards .eq. 0) then
         do i=1,natoms
            write(2,*) ntag(i),itype(i),rv(1,i),rv(2,i),rv(3,i)
         enddo
      else
         do i=1,natoms
            write(2,*) ntag(i),itype(i),nhittag(i),
     *         rv(1,i),rv(2,i),rv(3,i)
         enddo
      endif
      close(2)
      open(unit=2,file=dumpfile(1:index(dumpfile,' ')-1)//'.vel')
      write(2,*) 'ITEM: NUMBER OF ATOMS'
      write(2,*) natoms
      write(2,*) 'ITEM: BOUNDARY VELOCITY'
      write(2,*) zero,zero,zero
      write(2,*) 'ITEM: TIME'
      write(2,*) zero
      write(2,*) 'ITEM: VELOCITIES'
      do i=1,natoms
         write(2,*) ntag(i),rv(4,i),rv(5,i),rv(6,i)
      enddo
      close(2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write2(dumpfile)
      include 'createAtoms.h'
      data conmas/1.0365e-4/
      character*32 dumpfile
      open (unit=2,file=dumpfile)
      zero=0.0
      write(2,*) 'using dynamo'
      write(2,9502) natoms,ntypes,zero
9502  format(2i10,e15.8)
      write(2,9503) (perub(i),i=1,3),(perlb(i),i=1,3)
9503  format(3e25.16)
      write(2,9504) (amass(i)*conmas,ielement(i),i=1,ntypes)
9504  format(e25.16,i10)
      write(2,9505) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
9505  format(3e25.16/3e25.16/i10)
      write(2,9503) zero,zero,zero
      close(2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write3(dumpfile)
      include 'createAtoms.h'
      character*32 dumpfile
      open (unit=2,file=dumpfile)
      area=(1.0-float(nw_del)/float(natoms0))*perlen(2)*perlen(3)
      write(2,*) 'Start File for LAMMPS ',area
      write(2,*)
      write(2,*) natoms,' atoms'
      write(2,*)
      write(2,*) ntypes,' atom types'
      write(2,*)
      write(2,*) perlb(1),perub(1),' xlo xhi'
      write(2,*) perlb(2),perub(2),' ylo yhi'
      write(2,*) perlb(3),perub(3),' zlo zhi'
      if (xy .ge. 1.0e8 .or. xz .ge. 1.0e8 .or. yz .ge. 1.0e8) then
         write(2,*)
      else
         write(2,*)
         write(2,*)xy,xz,yz,' xy xz yz'
         write(2,*)
      endif
      write(2,*) 'Masses'
      write(2,*)
      do i=1,ntypes
         write(2,*)i,amass(i)
      enddo
      write(2,*)
      write(2,*) 'Atoms'
      write(2,*)
      do i=1,natoms
         write(2,*) i,itype(i),rv(1,i),rv(2,i),rv(3,i)
      enddo
      write(2,*)
      write(2,*) 'Velocities'
      write(2,*)
      zero=0.0
      do i=1,natoms
         write(2,*) i,rv(4,i),rv(5,i),rv(6,i)
      enddo
      close(2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write4(dumpfile)
      include 'createAtoms.h'
      character*32 dumpfile,full
      dumpfile='ens'
      full=dumpfile(1:index(dumpfile,' ')-1)//'.case'
      open(unit=2,file=full)
      write(2,101)
101   format('FORMAT')
      write(2,102)
102   format('type: ensight')
      write(2,103)
103   format('GEOMETRY')
      write(2,104)'ens.case.****.geo'
104   format('model: 1                     ',a30)
      write(2,105)
105   format('VARIABLE')
      write(2,111)
111   format('TIME')
      write(2,112)
112   format('time set: 1')
      write(2,113)
113   format('number of steps: 1')
      write(2,114)
114   format('filename start number: 1')
      write(2,115)
115   format('filename increment: 1')
      write(2,116)
116   format('time values:')
      write(2,117)
117   format('1')
      close(2)
      full=dumpfile(1:index(dumpfile,' ')-1)//'.case.0001.geo'
      open(unit=2,file=full)
      write(2,201)
      write(2,202)
      write(2,203)
      write(2,204)
      write(2,205)
      write(2,206)natoms
      do j=1,natoms
         write(2,207)j,rv(1,j),rv(2,j),rv(3,j)
      enddo
      do j=1,ntypes
         write(2,208)j
         write(2,209)j
         write(2,210)
         write(2,206)nntype(j)
         do k=1,natoms
            if (itype(k) .eq. j) write(2,211)k,k
         enddo
      enddo
201   format('spparks input')
202   format('geometry for element group 1')
203   format('node id given')
204   format('element id given')
205   format('coordinates')
206   format(i8)
207   format(i8,3e12.5)
208   format('part',i2)
209   format('atombox',i2)
210   format('point')
211   format(2i8)
      close(2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write5(dumpfile)
      include 'createAtoms.h'
      character*32 dumpfile,full
      character*200 sentence
      full=dumpfile(1:index(dumpfile,' ')-1)//'.lat'
      open (unit=2,file=full)
         call neigen()
         nmax=numneigh(1)
         do i=1,natoms
            if (nmax .lt. numneigh(i)) nmax=numneigh(i)
         enddo
         write(2,*) 'spparks input lattice'
         write(2,*) 
         write(2,*) '3 dimension'
         write(2,*) natoms,' vertices'
         write(2,*) nmax,' max connectivity'
         write(2,*) perlb(1),perub(1),' xlo xhi'
         write(2,*) perlb(2),perub(2),' ylo yhi'
         write(2,*) perlb(3),perub(3),' zlo zhi'
         write(2,*)
         write(2,*) 'Vertices'
         write(2,*)
         do i=1,natoms
            write(2,*) i,rv(1,i),rv(2,i),rv(3,i)
         enddo
         write(2,*)
         write(2,*) 'Edges'
         write(2,*)
         do i=1,natoms
            write(sentence,*) i,(neigh(k,i),k=1,numneigh(i))
            write(2,*)sentence(1:len_trim(sentence))
         enddo
      close(2)
      full=dumpfile(1:index(dumpfile,' ')-1)//'.conf'
      open (unit=2,file=full)
      nvalues=1
      write(2,*) 'spparks input configuration'
      write(2,*) natoms,nvalues
      write(2,*) 
      do i=1,natoms
         write(2,*) i,itype(i)
      enddo
      close(2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine disturb()
      include 'createAtoms.h'
      dimension strain(3)
      namelist /disturbcard/ dismax,strain
      dismax=0.0
      strain(1)=0.0
      strain(2)=0.0
      strain(3)=0.0
      read(5,disturbcard)
      if (dismax .ne. 0.0) then
         do 1 i=1,natoms
            rv(1,i)=rv(1,i)+dismax*(ranl()-0.5)
            rv(2,i)=rv(2,i)+dismax*(ranl()-0.5)
            rv(3,i)=rv(3,i)+dismax*(ranl()-0.5)
1        continue
      endif
      if (strain(1)**2+strain(2)**2+strain(3)**2 .ne. 0.0) then
         do 2 i=1,3
            perlen(i)=perlen(i)*(1.0+strain(i))
            perub(i)=perlb(i)+perlen(i)
2        continue
         do 3 i=1,natoms
         do 3 j=1,3
            rv(j,i)=perlb(j)+(rv(j,i)-perlb(j))*(1.0+strain(j))
3        continue
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outputfiles()
      character*32 dynamo,paradyn,lammps,xyz,ensight,spparks
      namelist /filecard/ dynamo,paradyn,lammps,xyz,ensight,spparks
      xyz="none"
      paradyn="none"
      dynamo="none"
      lammps="none"
      ensight="none"
      spparks="none"
      read(5,filecard)
      if (xyz .ne. "none") call write0(xyz)
      if (paradyn .ne. "none") call write1(paradyn)
      if (dynamo .ne. "none") call write2(dynamo)
      if (lammps .ne. "none") call write3(lammps)
      if (ensight .ne. "none") call write4(ensight)
      if (spparks .ne. "none") call write5(spparks)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function length(str)
      character*(*) str
      n = len(str)
 10   if (n.gt.0) then
        if (str(n:n).eq.' ') then
          n = n - 1
          goto 10
        endif
      endif
      length = n
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine velgen()
      include 'createAtoms.h'
      data boltz/8.617e-5/,conmas/1.0365e-4/
      namelist /velcard/ Tbnd,Tmid,naxis,nmesh,equal,ensureTave,
     *   steeper
      dimension Tmesh(200),pmesh(200),vcm(3,200),tmass(200),ekin(200),
     *   scale(200),ncount(200)
      Tbnd=0.0
      Tmid=0.0
      naxis=1
      nmesh=100
      equal=1.0
      ensureTave=0.0
      steeper=0.0
      read(5,velcard)
      Tbnd=equal*Tbnd
      Tmid=equal*Tmid
      Tave=0.5*(Tbnd+Tmid)
      Tbnd=Tave+(Tbnd-Tave)*(1.0+steeper)
      Tmid=Tave+(Tmid-Tave)*(1.0+steeper)

      if (Tbnd+Tmid .le. 0.0) then 
         do i=1,natoms
            rv(4,i)=0.0
            rv(5,i)=0.0
            rv(6,i)=0.0
         enddo
      else
         wmesh=perlen(naxis)/nmesh
         do i=1,nmesh
            pmesh(i)=(i-0.5)*wmesh
            Tmesh(i)=Tmid+(Tbnd-Tmid)*abs(pmesh(i)-0.5*perlen(naxis))/
     *         (0.5*perlen(naxis))
            vcm(1,i)=0.0
            vcm(2,i)=0.0
            vcm(3,i)=0.0
            tmass(i)=0.0
            ekin(i)=0.0
            ncount(i)=0
         enddo
         do i=1,natoms
            index=(rv(naxis,i)-perlb(naxis))/wmesh+1
            it=itype(i)
            vnorm=sqrt(Tmesh(index)*boltz/(amass(it)*conmas))
            do j=1,3
               rv(j+3,i)=vnorm*rgauss()
               vcm(j,index)=vcm(j,index)+amass(it)*conmas*rv(j+3,i)
            enddo
            tmass(index)=tmass(index)+amass(it)*conmas
            ncount(index)=ncount(index)+1
         enddo
         do i=1,nmesh
            do j=1,3
               vcm(j,i)=vcm(j,i)/tmass(i)
            enddo
         enddo
         do i=1,natoms
            index=(rv(naxis,i)-perlb(naxis))/wmesh+1
            do j=1,3
               rv(j+3,i)=rv(j+3,i)-vcm(j,index)
            enddo
         enddo
         do i=1,natoms
            index=(rv(naxis,i)-perlb(naxis))/wmesh+1
            it=itype(i)
            ekin(index)=ekin(index)+0.5*amass(it)*conmas*
     *         (rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)
         enddo
         do i=1,nmesh
            treal=2.0/(3.0*boltz)*ekin(i)/ncount(i)
            if (treal .eq. 0.0) then
               scale(i)=0.0
            else
               scale(i)=sqrt(Tmesh(i)/treal)
            endif
         enddo
         do i=1,natoms
            index=(rv(naxis,i)-perlb(naxis))/wmesh+1
            do j=4,6
               rv(j,i)=rv(j,i)*scale(index)
            enddo
         enddo
         if (ensureTave .eq. 1.0) then
            ekint=0.0
            do i=1,natoms
               it=itype(i)
               ekint=ekint+0.5*amass(it)*conmas*
     *            (rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)
            enddo
            trealt=2.0/(3.0*boltz)*ekint/natoms
            if (trealt .eq. 0.0) then
               scalet=0.0
            else
               scalet=sqrt(Tave/trealt)
            endif
            do i=1,natoms
               do j=4,6
                  rv(j,i)=rv(j,i)*scalet
               enddo
            enddo
         endif
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function rgauss()
      data twopi/6.283185308/
1     continue
      a1 = ranl()
      if (a1 .eq. 0.0) goto 1
      a2 = ranl()
      rgauss = sqrt(-2.0*log(a1))*cos(twopi*a2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine neigen()
      include 'createAtoms.h'
      dimension idbox(100,100,100,100),nbox(100,100,100)
      cutoff=100.0
      xmin=rv(1,1)
      xmax=rv(1,1)
      ymin=rv(2,1)
      ymax=rv(2,1)
      zmin=rv(3,1)
      zmax=rv(3,1)
      do i=2,natoms
         if (xmin .gt. rv(1,i)) xmin=rv(1,i)
         if (xmax .lt. rv(1,i)) xmax=rv(1,i)
         if (ymin .gt. rv(2,i)) ymin=rv(2,i)
         if (ymax .lt. rv(2,i)) ymax=rv(2,i)
         if (zmin .gt. rv(3,i)) zmin=rv(3,i)
         if (zmax .lt. rv(3,i)) zmax=rv(3,i)
         dx=rv(1,i)-rv(1,1)
         dy=rv(2,i)-rv(2,1)
         dz=rv(3,i)-rv(3,1)
         dx=dx-perlen(1)*nint(dx/perlen(1))
         dy=dy-perlen(2)*nint(dy/perlen(2))
         dz=dz-perlen(3)*nint(dz/perlen(3))
         tmp=dx**2+dy**2+dz**2
         if (tmp .lt. cutoff) cutoff=tmp
      enddo
      cutoff=sqrt(cutoff)+0.1
      nx=(xmax-xmin)/cutoff+1
      ny=(ymax-ymin)/cutoff+1
      nz=(zmax-zmin)/cutoff+1
      if (nx .gt. 100) nx=100
      if (ny .gt. 100) ny=100
      if (nz .gt. 100) nz=100
      delx=(xmax-xmin)/nx
      dely=(ymax-ymin)/ny
      delz=(zmax-zmin)/nz
      do n1=1,nx
      do n2=1,ny
      do n3=1,nz
         nbox(n1,n2,n3)=0
      enddo
      enddo
      enddo
      do j=1,natoms
         numneigh(j)=0
         n1=(rv(1,j)-xmin)/delx+1
         if (n1 .lt. 1) n1=1
         if (n1 .gt. nx) n1=nx
         n2=(rv(2,j)-ymin)/dely+1
         if (n2 .lt. 1) n2=1
         if (n2 .gt. ny) n2=ny
         n3=(rv(3,j)-zmin)/delz+1
         if (n3 .lt. 1) n3=1
         if (n3 .gt. nz) n3=nz
         nbox(n1,n2,n3)=nbox(n1,n2,n3)+1
         idbox(n1,n2,n3,nbox(n1,n2,n3))=j
      enddo
      do j=1,natoms
         n1=(rv(1,j)-xmin)/delx+1
         if (n1 .lt. 1) n1=1
         if (n1 .gt. nx) n1=nx
         n2=(rv(2,j)-ymin)/dely+1
         if (n2 .lt. 1) n2=1
         if (n2 .gt. ny) n2=ny
         n3=(rv(3,j)-zmin)/delz+1
         if (n3 .lt. 1) n3=1
         if (n3 .gt. nz) n3=nz
      do na=n1-1,n1+1
         if (na .lt. 1) then
            naa=nx
         else if (na .gt. nx) then
            naa=1
         else
            naa=na
         endif
      do nb=n2-1,n2+1
         if (nb .lt. 1) then
            nbb=ny
         else if (nb .gt. ny) then
            nbb=1
         else
            nbb=nb
         endif
      do nc=n3-1,n3+1
         if (nc .lt. 1) then
            ncc=nz
         else if (nc .gt. nz) then
            ncc=1
         else
            ncc=nc
         endif
      do 100 ndd=1,nbox(naa,nbb,ncc)
         k=idbox(naa,nbb,ncc,ndd)
         if (k .ge. j) goto 100
         dx=rv(1,k)-rv(1,j)
         dy=rv(2,k)-rv(2,j)
         dz=rv(3,k)-rv(3,j)
         dx=dx-perlen(1)*nint(dx/perlen(1))
         dy=dy-perlen(2)*nint(dy/perlen(2))
         dz=dz-perlen(3)*nint(dz/perlen(3))
         tmp=sqrt(dx**2+dy**2+dz**2)
         if (tmp .le. cutoff) then
            numneigh(j)=numneigh(j)+1
            neigh(numneigh(j),j)=k
            numneigh(k)=numneigh(k)+1
            neigh(numneigh(k),k)=j
         endif
100   continue
      enddo
      enddo
      enddo
      enddo
      return
      end

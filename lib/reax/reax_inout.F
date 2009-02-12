**********************************************************************
*                                                                    *
*     REAXFF Reactive force field program                            *
*                                                                    *
*     Developed and written by Adri van Duin, duin@wag.caltech.edu   *
*                                                                    *
*     Copyright (c) 2001-2010 California Institute of Technology     *
*                                                                    *
*     This is an open-source program. Feel free to modify its        *
*     contents. Please keep me informed of any useful modification   *
*     or addition that you made. Please do not distribute this       *
*     program to others; if people are interested in obtaining       *
*     a copy of this program let them contact me first.              *
*                                                                    *
**********************************************************************
**********************************************************************

      subroutine ffinpt

**********************************************************************
#include "cbka.blk"
#include "cbkboncor.blk"
#include "cbkconst.blk"
#include "cbkcovbon.blk"
#include "cbkff.blk"
#include "cbkfftorang.blk"
#include "cbknonbon.blk"
#include "cbksrthb.blk"
#include "cbktorsion.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "valang.blk"
#include "cbksrtbon1.blk"
#include "cbkchb.blk"
      dimension rcore2(nsort),ecore2(nsort),acore2(nsort)
**********************************************************************
*                                                                    *
*     Read in force field                                            *
*                                                                    *
**********************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In ffinpt'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      open (4,file='ffield.reax',status='old')
      rewind (4)
      iline=0
      read (4,'(a40)',end=990,err=990)qffield
      iline=iline+1
**********************************************************************
*                                                                    *
*     Read in general force field parameters                         *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990)npar
      iline=iline+1
      do i1=1,npar
      read (4,1300,end=990,err=990)vpar(i1)
      iline=iline+1
      end do
      cutoff=0.01*vpar(30)
      swa=vpar(12)
      if (abs(swa).gt.0.01) write (*,*)
     $'Warning: non-zero value for lower Taper-radius cutoff'
      swb=vpar(13)
      if (swb.lt.zero) stop 
     $'Negative value for upper Taper-radius cutoff'
      if (swb.lt.5.0) write (*,*)
     $'Warning: very low value for upper Taper-radius cutoff:',swb
**********************************************************************
*                                                                    *
*     Read in atom type data                                         *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990) nso
      iline=iline+1
      read (4,*,end=990,err=990)
      iline=iline+1
      read (4,*,end=990,err=990)
      iline=iline+1
      read (4,*,end=990,err=990)
      iline=iline+1
      if (nso.gt.nsort) stop 'Maximum number of atom types exceeded'
      do i1=1,nso
      read (4,1200,end=990,err=990)qas(i1),rat(i1),aval(i1),amas(i1),
     $rvdw(i1),eps(i1),gam(i1),rapt(i1),stlp(i1)
      iline=iline+1
      read (4,1250,end=990,err=990)alf(i1),vop(i1),valf(i1),
     $valp1(i1),valp2(i1),chi(i1),eta(i1),vnphb
      iline=iline+1
      read (4,1250,end=990,err=990)vnq(i1),vlp1(i1),vincr(i1),
     $bo131(i1),bo132(i1),bo133(i1),sigqeq(i1),default
      iline=iline+1
      read (4,1250,end=990,err=990)vovun(i1),vval1(i1),vrom,
     $vval3(i1),vval4(i1),rcore2(i1),ecore2(i1),acore2(i1)
      iline=iline+1
      idef(i1)=int(default)
      nphb(i1)=int(vnphb)
      end do
**********************************************************************
*                                                                    *
*     Calculate van der Waals and Coulomb pair-parameters            *
*                                                                    *
**********************************************************************
      do i1=1,nso
      do i2=1,nso
      rcore(i1,i2)=sqrt(rcore2(i1)*rcore2(i2))
      ecore(i1,i2)=sqrt(ecore2(i1)*ecore2(i2))
      acore(i1,i2)=sqrt(acore2(i1)*acore2(i2))
      p1co(i1,i2)=sqrt(4.0*rvdw(i1)*rvdw(i2))
      p2co(i1,i2)=sqrt(eps(i1)*eps(i2))
      p3co(i1,i2)=sqrt(alf(i1)*alf(i2))
      gamwh=sqrt(vop(i1)*vop(i2))
      gamwco(i1,i2)=1.0/gamwh**vpar(29)
      gamch=sqrt(gam(i1)*gam(i2))
      gamcco(i1,i2)=1.0/gamch**3
      rob1(i1,i2)=0.50*(rat(i1)+rat(i2))
      rob2(i1,i2)=0.50*(rapt(i1)+rapt(i2))
      rob3(i1,i2)=0.50*(vnq(i1)+vnq(i2))
      end do
      end do
**********************************************************************
*                                                                    *
*     Read in bond type data                                         *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990)nboty
      iline=iline+1
      read (4,*,end=990,err=990)
      iline=iline+1
      if (2*nboty.gt.nbotym) stop 'Maximum nr. of bond types exceeded'
      ih=0
      do i1=1,nboty
      ih=ih+1
      read (4,1400,end=990,err=990)nbs(ih,1),nbs(ih,2),de1(ih),
     $de2(ih),de3(ih),psi(ih),pdo(ih),v13cor(ih),popi(ih),vover(ih)
      iline=iline+1
      read (4,1450,end=990,err=990)psp(ih),pdp(ih),ptp(ih),
     $bom(ih),bop1(ih),bop2(ih),ovc(ih),vuncor(ih)
      iline=iline+1
      if (nbs(ih,1).ne.nbs(ih,2)) then
      ih=ih+1
      nbs(ih,1)=nbs(ih-1,2)
      nbs(ih,2)=nbs(ih-1,1)
      de1(ih)=de1(ih-1)
      de2(ih)=de2(ih-1)
      de3(ih)=de3(ih-1)
      psi(ih)=psi(ih-1)
      pdo(ih)=pdo(ih-1)
      v13cor(ih)=v13cor(ih-1)
      vover(ih)=vover(ih-1)
      psp(ih)=psp(ih-1)
      pdp(ih)=pdp(ih-1)
      ptp(ih)=ptp(ih-1)
      bop1(ih)=bop1(ih-1)
      bop2(ih)=bop2(ih-1)
*     bop3(ih)=bop3(ih-1)
*     bop4(ih)=bop4(ih-1)
      bom(ih)=bom(ih-1)
      popi(ih)=popi(ih-1)
      ovc(ih)=ovc(ih-1)
      end if
      end do
      nboty2=ih
**********************************************************************
*                                                                    *
*     Read in off-diagonal parameters                                *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990)nodmty
      iline=iline+1
      if (nodmty.gt.nodmtym)
     $stop 'Maximum nr. of off-diagonal Morse types exceeded'
      ih=0
      do i1=1,nodmty
      ih=ih+1
      read (4,1400,end=990,err=990)nodm1,nodm2,deodmh,rodmh,godmh,
     $rsig,rpi,rpi2
      iline=iline+1
      if (rsig.gt.zero) rob1(nodm1,nodm2)=rsig
      if (rsig.gt.zero) rob1(nodm2,nodm1)=rsig
      if (rpi.gt.zero) rob2(nodm1,nodm2)=rpi
      if (rpi.gt.zero) rob2(nodm2,nodm1)=rpi
      if (rpi2.gt.zero) rob3(nodm1,nodm2)=rpi2
      if (rpi2.gt.zero) rob3(nodm2,nodm1)=rpi2
      if (rodmh.gt.zero) p1co(nodm1,nodm2)=2.0*rodmh
      if (rodmh.gt.zero) p1co(nodm2,nodm1)=2.0*rodmh
      if (deodmh.gt.zero) p2co(nodm1,nodm2)=deodmh
      if (deodmh.gt.zero) p2co(nodm2,nodm1)=deodmh
      if (godmh.gt.zero) p3co(nodm1,nodm2)=godmh
      if (godmh.gt.zero) p3co(nodm2,nodm1)=godmh
      end do
**********************************************************************
*                                                                    *
*     Read in valency angle and conjugation type data                *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990)nvaty
      iline=iline+1
      if (nvaty.gt.nvatym)
     $stop 'Maximum nr. of valency angle types exceeded'
      do i1=1,nvaty
      read (4,1500,end=990,err=990)nvs(i1,1),nvs(i1,2),
     $nvs(i1,3),th0(i1),vka(i1),vka3(i1),vka8(i1),vkac(i1),vkap(i1),
     $vval2(i1)
      iline=iline+1
      end do
**********************************************************************
*                                                                    *
*     Read in torsion angle type data                                *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990)ntoty
      iline=iline+1
      if (ntoty.gt.ntotym)
     $stop 'Maximum nr. of torsion angle types exceeded'
      do i1=1,ntoty
      read (4,1600,end=990,err=990)nts(i1,1),nts(i1,2),nts(i1,3),
     $nts(i1,4),v1(i1),
     $v2(i1),v3(i1),v4(i1),vconj(i1),v2bo(i1),v3bo(i1)
      iline=iline+1
      end do
**********************************************************************
*                                                                    *
*     Read in hydrogen bond type data                                *
*                                                                    *
**********************************************************************
      read (4,1100,end=990,err=990)nhbty
      iline=iline+1
      if (nhbty.gt.nhbtym)
     $stop 'Maximum nr. of hydrogen bond types exceeded'
      do i1=1,nhbty
      read (4,1500,end=990,err=990)nhbs(i1,1),nhbs(i1,2),
     $nhbs(i1,3),rhb(i1),dehb(i1),vhb1(i1),vhb2(i1)
      iline=iline+1
      end do
**********************************************************************
*                                                                    *
*     Calculate vdWaals interaction parameters                       *
*                                                                    *
**********************************************************************
      do i1=1,nso
      do i2=1,nso
      rr=(rvdw(i1)+rvdw(i2))
      rr2=rr*rr
      eps2=sqrt(eps(i1)*eps(i2))
      rr6=rr2*rr2*rr2
      pvdw1(i1,i2)=eps2*rr6*rr6
      pvdw1(i2,i1)=eps2*rr6*rr6
      pvdw2(i1,i2)=2.0*eps2*rr6
      pvdw2(i2,i1)=2.0*eps2*rr6
      end do
      end do
**********************************************************************
*                                                                    *
*     Error part                                                     *
*                                                                    *
**********************************************************************
      goto 999
  990 write (*,*)'Error or end-of-file reading unit 4 on line:',iline
      stop
  999 continue
      close(4)
**********************************************************************
*                                                                    *
*     Format part                                                    *
*                                                                    *
**********************************************************************
 1100 format (i3,2x,a2,3x,3d22.15)
 1200 format (1x,a2,10f9.4)
 1250 format (3x,10f9.4)
 1300 format (f10.4)
 1400 format (2i3,8f9.4)
 1450 format (6x,8f9.4)
 1500 format (3i3,7f9.4)
 1600 format (4i3,7f9.4)
      return
      end
**********************************************************************
***********************************************************************

      subroutine mdsav(node)

***********************************************************************
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkatomcoord.blk"
#include "cbkbo.blk"
#include "cbkc.blk"
#include "cbkch.blk"
#include "cbkconst.blk"
#include "cbkdistan.blk"
#include "cbkenergies.blk"
#include "cbkia.blk"
#include "cbkinit.blk"
#include "cbklonpar.blk"
#include "cbknubon2.blk"
#include "cbkqa.blk"
#include "cbktregime.blk"
#include "cbksrtbon1.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"

      dimension idum(mbond+3),bodum(mbond+3),qat2(2)
      character*25 qfileh
      character*33 qfile2
      character*4 qext
      character*6 qmdfi
      character *7 var
      character *3 qat2,pepname
      character *1 qrtemp
************************************************************************
*                                                                      *
*     Save coordinates, velocities and accelerations of MD-system      *
*                                                                      *
************************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In mdsav'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

************************************************************************
c
c     This is just for test purposes
c
************************************************************************
c$$$      write(6,*) '***************************'
c$$$      write(6,*) 'mdsav node number is ',node
c$$$      write(6,*) '***************************'
      return


      qfileh='Unknown'
      qmdfi='moldyn'
      pepname='   '
      ipeptide=0
      if (ni.eq.2) qmdfi='molsav'

      if (iopt.eq.0) then

      do i1=1,mbond+3
      idum(i1)=nzero
      bodum(i1)=zero
      end do
C      if (napp.eq.1) 
C     $open (7,file='fort.7',status='unknown',access='append')
      if (napp.ne.1) 
     $open (7,file='fort.7',status='unknown')
      nsbmaxh=5*((nsbmax/5)+1)
      write (7,100)na,qmol,mdstep,nsbmaxh
      if (nbiolab.eq.1) write (67,101)na,qmol
      do i1=1,na
      bosum=0.0
      do i3=1,nsbmax
      if (iag(i1,2+i3).gt.0) bosum=bosum+bo(nubon1(i1,i3))
      end do
      if (nsbmax.lt.5) then
      write (7,200)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)),
     $(idum(i2),i2=1,5-iag(i1,2)),
     $iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)),
     $(bodum(i2),i2=1,5-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
      if (nbiolab.eq.1) then      !Delphi-connection table output
      write (67,201)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2))
      end if
      else if (nsbmax.lt.10) then
      write (7,210)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)),
     $(idum(i2),i2=1,10-iag(i1,2)),
     $iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)),
     $(bodum(i2),i2=1,10-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbmax.lt.15) then
      write (7,220)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)),
     $(idum(i2),i2=1,15-iag(i1,2)),
     $iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)),
     $(bodum(i2),i2=1,15-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbmax.lt.20) then
      write (7,230)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)),
     $(idum(i2),i2=1,20-iag(i1,2)),
     $iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)),
     $(bodum(i2),i2=1,20-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbmax.lt.25) then
      write (7,240)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)),
     $(idum(i2),i2=1,25-iag(i1,2)),
     $iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)),
     $(bodum(i2),i2=1,25-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbmax.gt.25) then
      write (7,250)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)),
     $(idum(i2),i2=1,35-iag(i1,2)),
     $iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)),
     $(bodum(i2),i2=1,35-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
      end if
      end do
      boss=zero
      vlps=0.0
C      if (napp.eq.1)
C     $open (8,file='fort.8',status='unknown',access='append')
      if (napp.ne.1)
     $open (8,file='fort.8',status='unknown')
      nsbmaxh=5*((nsbma2/5)+1)
      write (8,100)na,qmol,mdstep,nsbmaxh
      chsum=0.0
      do i1=1,na
      bosum=0.0
      do i3=1,nsbma2
      if (ia(i1,2+i3).gt.0) bosum=bosum+bo(nubon2(i1,i3))
      end do
      if (nsbma2.lt.5) then
      write (8,200)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)),
     $(idum(i2),i2=1,5-ia(i1,2)),
     $ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)),
     $(bodum(i2),i2=1,5-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbma2.lt.10) then
      write (8,210)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)),
     $(idum(i2),i2=1,10-ia(i1,2)),
     $ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)),
     $(bodum(i2),i2=1,10-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbma2.lt.15) then
      write (8,220)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)),
     $(idum(i2),i2=1,15-ia(i1,2)),
     $ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)),
     $(bodum(i2),i2=1,15-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbma2.lt.20) then
      write (8,230)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)),
     $(idum(i2),i2=1,20-ia(i1,2)),
     $ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)),
     $(bodum(i2),i2=1,20-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
      else if (nsbma2.lt.25) then
      write (8,240)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)),
     $(idum(i2),i2=1,25-ia(i1,2)),
     $ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)),
     $(bodum(i2),i2=1,25-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
      end if
      boss=boss+bosum/2.0
      vlps=vlps+vlp(i1)
      chsum=chsum+ch(i1)
      end do
      write (7,*)2.0*boss,vlps,2.0*boss+2.0*vlps,chsum
      close(8)
      close(7)

      end if

      if (noutpt.eq.0) then
      write (var,'(f7.4)')float(mdstep/nsav)/1d4
      if (ni.eq.0) open (unit=67,file=qmdfi//var(3:7),
     $status='unknown')
      write (67,300)qmol
      do i1=1,na
      write (67,400)i1,qa(i1),(c(i1,i2),i2=1,3)
      end do
      write (67,*)
      close(67)
      end if

      if (noutpt.eq.2) then
C      open (88,file='moldyn.bgf',status='unknown',access='append')
      call writebgf(88)
      close (88)
      end if

      if ((ni.eq.1.and.iopt.eq.0).or.(ni.eq.1.and.iopt.eq.1.and.
     $iflga.eq.1)) then
      qrtemp=qr
      if (qr.eq.'I') qr='C'
      if (qfileh.eq.' ') then
      write (*,*)'Warning: no file name given; use Unknown' 
      qfileh='Unknown'
      end if
      qfile2=qfileh
      if (imodfile.eq.0) then
      istart=1
      qstrana1(1:25)=qfileh
      call stranal(istart,iend,vout,iout,1)
      qfile2=qfileh(istart:iend-1)//".geo"
      end if
      call writegeo(98)

      if (imodfile.eq.1.or.iopt.eq.0) then
      open (88,file=qfile2,status='unknown')
      call writegeo(88)
      close (88)
      end if 

      qr=qrtemp

      if (iopt.eq.0) then

      do i1=1,na
      write (56,410) i1,ch(i1)
      write (55,410) i1,chgbgf(i1)
      end do
**********************************************************************
*                                                                    *
*     Write .pdb output file                                         *
*                                                                    *
**********************************************************************
      open (unit=47,file='output.pdb',status='unknown')
      do i1=1,na
      write (47,412)'ATOM  ',i1,qa(i1),pepname,ipeptide,c(i1,1),
     $c(i1,2),c(i1,3),1.0,2.2,qa(i1)
      end do
      write (47,*) 'TER'
      write (47,*) 'END'
      close (47)

      if (nsurp.eq.0) then
      if (kx.gt.0.or.ky.gt.0.or.kz.gt.0) then
      qrtemp=qr
**********************************************************************
*                                                                    *
*     Write crystal structure including periodic images              *
*                                                                    *
**********************************************************************
*     mux=(1+kx+kx)
*     muy=(1+ky+ky)
*     muz=(1+kz+kz)
*     qr='F'
*     write (86,'(2x,a1,1x,a60)')qr,qmol
*     qr=qrtemp
*     write (86,'(3f10.4)')mux*axiss(1),muy*axiss(2),muz*axiss(3)
*     write (86,'(3f10.4)')angle(1),angle(2),angle(3)
*     do i1=1,na
*     write (86,'(i4,1x,a2,3x,3d22.15)')i1,qa(i1),(c(i1,i2),i2=1,3)
*     end do
*     nhulp=na+1
*     do k1=-kx,kx
*     do k2=-ky,ky
*     do k3=-kz,kz
*     if (k1.ne.0.or.k2.ne.0.or.k3.ne.0) then
*     do i1=1,na
*     cx=c(i1,1)+k1*tm11
*     cy=c(i1,2)+k1*tm21+k2*tm22
*     cz=c(i1,3)+k1*tm31+k2*tm32+k3*tm33
*     write (86,'(i4,1x,a2,3x,3d22.15)')nhulp,qa(i1),cx,cy,cz
*     nhulp=nhulp+1
*     end do
*     end if
*     end do
*     end do
*     end do
*     write (86,*)
**********************************************************************
*                                                                    *
*     Write crystal structure with extra unit cells                  *
*                                                                    *
**********************************************************************
      mux=1+iexx
      muy=1+iexy
      muz=1+iexz
      qr='F'
      write (85,'(2x,a1,1x,a60)')qr,qmol
      qr=qrtemp
      write (85,'(3f10.4)')mux*axiss(1),muy*axiss(2),muz*axiss(3)
      write (85,'(3f10.4)')angle(1),angle(2),angle(3)
      do i1=1,na
      write (85,'(i4,1x,a2,3x,3d22.15)')i1,qa(i1),(c(i1,i2),i2=1,3)
      end do
      nhulp=na+1
      do k1=0,iexx
      do k2=0,iexy
      do k3=0,iexz
      if (k1.ne.0.or.k2.ne.0.or.k3.ne.0) then
      do i1=1,na
      cx=c(i1,1)+k1*tm11
      cy=c(i1,2)+k1*tm21+k2*tm22
      cz=c(i1,3)+k1*tm31+k2*tm32+k3*tm33
      write (85,'(i4,1x,a2,3x,3d22.15)')nhulp,qa(i1),cx,cy,cz
      nhulp=nhulp+1
      end do
      end if
      end do
      end do
      end do
      write (85,*)

      end if
      end if
      end if

      end if

      if (ni.eq.0.or.ni.eq.2) then
********************************************************************** 
*                                                                    *
*     Write ASCII trajectory file                                    *
*                                                                    *
********************************************************************** 
      if (ni.eq.0) open(unit=66,file=qmdfi//'.vel',status='unknown')
      if (ni.eq.2) then
      write (var,'(f7.4)')float(mdstep/nsav3)/1d4
      open (unit=66,file=qmdfi//var(3:7),status='unknown')
      end if
      write (66,500)axis(1),axis(2),axis(3)
      write (66,550)angle(1),angle(2),angle(3)
      write (66,600)na,((c(i,j),j=1,3),qlabel(i),i=1,na)
      write (66,700)((vel(j,i),j=1,3),i=1,na)
      write (66,800)((accel(j,i),j=1,3),i=1,na)
      write (66,900)((aold(j,i),j=1,3),i=1,na)
      write (66,1000)tempmd
      write (66,1050)
      close (66)
      end if
      if (ni.ne.2.and.iopt.eq.0) then

C      open (unit=68,file='xmolout',status='unknown',access='append')
      write (68,1200)na
      write (68,1300)qmol,mdstep+nit+nprevrun,estrc,
     $axis(1),axis(2),axis(3),angle(1),angle(2),angle(3)
      do i1=1,na
      if (ixmolo.eq.0) write (68,1400)qa(i1),(c(i1,i2),i2=1,3)
      if (ixmolo.eq.1) write (68,1400)qa(i1),(c(i1,i2),i2=1,3),
     $(vel(i2,i1)/1e+10,i2=1,3),iag(i1,3+mbond)
      if (ixmolo.eq.2) write (68,1401)qa(i1),(c(i1,i2),i2=1,3),
     $iag(i1,3+mbond)
      end do
      close (68)

      if (itrout.ne.0) then
C      open (unit=69,file='xmolout2',status='unknown',access='append')
      write (69,1200)na
      write (69,1300)qmol,mdstep+nit+nprevrun,estrc,
     $axis(1),axis(2),axis(3),angle(1),angle(2),angle(3)
      do i1=1,na
      if (ixmolo.eq.0) write (69,1400)qa(i1),(cp(i1,i2),i2=1,3)
      if (ixmolo.eq.1) write (69,1400)qa(i1),(cp(i1,i2),i2=1,3),
     $(vel(i2,i1)/1e+10,i2=1,3),iag(i1,3+mbond)
      if (ixmolo.eq.2) write (68,1401)qa(i1),(c(i1,i2),i2=1,3),
     $iag(i1,3+mbond)
      end do
      close (69)
      end if

      call molanal
      end if
********************************************************************** 
*                                                                    *
*     Generate BIOGRAF output-file                                   *
*                                                                    *
********************************************************************** 
      if ((ni.eq.1.and.iopt.eq.0).or.(ni.eq.1.and.iopt.eq.1.and.
     $iflga.eq.1)) then

      if (qfileh.eq.' ') then
      write (*,*)'Warning: no file name given; use Unknown' 
      qfileh='Unknown'
      end if
      qfile2=qfileh
      if (imodfile.eq.0) then
      istart=1
      qstrana1(1:25)=qfileh
      call stranal(istart,iend,vout,iout,1)
      qfile2=qfileh(istart:iend-1)//".bgf"
      end if
      call writebgf(90)

      if (imodfile.eq.1.or.iopt.eq.0) then
      open (88,file=qfile2,status='unknown')
      call writebgf(88)
      close (88)
      end if

      end if
     
      return
********************************************************************** 
*                                                                    *
*     Format part                                                    *
*                                                                    *
********************************************************************** 
  100 format (i4,1x,a40,'Iteration:',i8,' #Bonds:',i4)       
  101 format (i3,2x,a40)
  200 format (8i4,8f7.3)
  201 format (8i3)
  210 format (13i4,13f7.3)
  220 format (18i4,18f7.3)
  230 format (23i4,23f7.3)
  240 format (28i4,28f7.3)
  250 format (38i4,38f7.3)
  300 format (2x,a1,1x,a60)
  301 format (2x,a1,1x,f6.2,a60)
  302 format (2x,a1,1x,2f6.2,a60)
  310 format (2x,a1,1x,a60)
  320 format (3f10.4)
  400 format (i4,1x,a2,3x,3(d21.14,1x),1x,a5,1x,i5)
  410 format (i4,f12.6)
  412 format(A6,I5,1x,A2,3x,A3,2x,i4,4x,3f8.3,f6.2,f6.2,4x,2x,A6)
  500 format (1x,'Lattice parameters:',/(3f15.8))
  550 format (3f15.8)
  600 format (i4,1x,'Atom coordinates (Angstrom):',/
     $(3d24.15,1x,a5))
  700 format (1x,'Atom velocities (Angstrom/s):',/(3d24.15))
  800 format (1x,'Atom accelerations (Angstrom/s**2):',/(3d24.15))
  900 format (1x,'Previous atom accelerations:',/(3d24.15))
 1000 format (1x,'MD-temperature (K):',/(1d24.15))
 1050 format (1x,'Connections, bond orders and lone pairs:')
 1100 format (8i3,8f8.4)
 1200 format (i4)
 1300 format (a40,i6,f12.4,6f7.2)
 1400 format (a2,3f10.5,3f15.5,i6)
 1401 format (a2,3f10.5,i6)
 1500 format ('BIOGRF',i4)
 1600 format ('XTLGRF',i4)
 1700 format ('DESCRP ',a60)
 1800 format ('REMARK ',a60)
 1900 format ('FFIELD ',a40)
 2000 format ('RUTYPE ',a40)
 2100 format ('CRYSTX ',6f11.5)
 2200 format ('CELLS ',6i5)
 2300 format ('#              At1 At2   R12    Force1  Force2  ',
     $'dR12/dIteration(MD only)')
 2400 format ('BOND RESTRAINT ',2i4,f8.4,f8.2,f8.5,f10.7)
 2500 format ('#               At1 At2 At3 Angle   Force1  Force2',
     $'  dAngle/dIteration (MD only)')
 2600 format ('ANGLE RESTRAINT ',3i4,2f8.2,f8.4,f9.6)
 2700 format ('#                 At1 At2 At3 At3 Angle   Force1  ',
     $'Force2  dAngle/dIteration (MD only)')
 2800 format ('TORSION RESTRAINT ',4i4,2f8.2,f8.4,f9.6)
 2900 format ('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,',
     $'3f10.5,1x,a5,i3,i2,1x,f8.5)')
 3000 format ('HETATM',1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,
     $a5,i3,i2,1x,f8.5)
 3100 format ('FORMAT CONECT (a6,12i6)')
 3200 format ('CONECT',12i6)
 3300 format ('UNIT ENERGY   kcal')
 3400 format ('ENERGY',5x,f14.6)
 3500 format ('END')
      end

************************************************************************
************************************************************************

      subroutine readc

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkcha.blk"
#include "cbkconst.blk"
#include "cbkdistan.blk"
#include "cbkinit.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"

      character*6 qident
      character*20 qhulp
*     dimension qident(100)
************************************************************************
*                                                                      *
*     Read control file                                                *
*                                                                      *
************************************************************************
c$$$c      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$c      write (65,*) 'In readc'
c$$$c      call timer(65)
c$$$c      close (65)
c$$$c      end if
      if (mdstep.gt.0.or.nit.gt.0) nmmsav=nmm
************************************************************************
*                                                                      *
*     Set default values                                               *
*                                                                      *
************************************************************************
      nreac=0
      axis1=200.0d0
      axis2=200.0d0
      axis3=200.0d0
      cutof2=0.001d0
      cutof3=0.300d0
      tsetor=298.0d0
      tset2=298.0d0
      pset=0.0d0
      tincr=0.0d0
      tstep=0.5d0
*     swa=0.0   !Moved to force field
*     swb=12.5  !Moved to force field
      taut=2.5d0
      taut2=2.5d0
      ndtau=50000
      taup=500.0d0
      vqnd=100.0d0
      errnh=1.0d0
      range=2.5d0
      maxstp=1000
      nequi=0
      nmethod=3
      ncha=3
      ncha2=1
      nchaud=1
      nvlist=25
      nrep1=5
      nsav=50
      icheck=0
      ivels=0
      itfix=0
      ncontrol=25
      noutpt=0
      napp=0
      nsurp=0
      ncons=2
      nrand=0
      nmm=0
      endpo=1.0d0
      endpo2=1.0d0
      nfc=50
      nsav2=50
      nmmax=50
      i5758=0
      parc1=1.0d0
      parc2=0.001d0
      icell=0
      ingeo=1
      ccpar=1.0005d0
      icelo2=0
      nrdd=0
      nrddf=200000
      nbiolab=0
c      ngeofor=0
      nincrop=0
      accerr=2.50d0
      vrange=2.50d0
      vlbora=5.00d0
      nsav3=1000
      nhop2=25
      nprevrun=0
      ndebug=0
      volcha=10.00d0
      ixmolo=0
      inpt=0
      iconne=0
      imolde=0
      ianaly=0
      icentr=0
      itrans=0
      itrout=0
      tpnrad=300.0d0
      ityrad=3
      iexx=1
      iexy=1
      iexz=1
      syscha=0.00d0
      inmov1=0
      inmov2=0
      vfield=0.00d0
      itstep=0
      ifreq=0
      isymm=1
      icpres=0
      delvib=0.0001d0
c     shock variables
      shock_vel = 2.d0 ! impact velocity for shock simulations (nm/ps)
      shock_z_sep = 10.0d0 ! separation z value to apply initial velocities in shocks
      ishock_type = 0.0d0 ! shock type. 0: simple impact; 1: compressing c axis
c     Hugoniostat variables
      Hug_E0 = 0.d0 ! Reference energy
      Hug_P0 = 0.d0 ! Reference pressure
      Hug_V0 = 0.d0 ! Reference volume
c     Shear flow simulations for viscosity
      xImpVcm = 1.d0 ! velocity applied in shear simulations (in nm/ps), left half mover at -xImpVcm and right at +xImpVcm
c$$$************************************************************************
c$$$*                                                                      *
c$$$*     Read control-file                                                *
c$$$*                                                                      *
c$$$************************************************************************
c$$$      open (10,file='control',status='old')
c$$$   10 read (10,'(a20)',end=20,err=30)qhulp
c$$$      if (qhulp(1:1).eq.'#') goto 10
c$$$      read (qhulp,*,err=30)vhulp
c$$$      read (qhulp,'(8x,a6)',err=30)qident
c$$$      if (qident.eq.'Hug_V0') Hug_P0=vhulp
c$$$      if (qident.eq.'Hug_P0') Hug_V0=vhulp
c$$$      if (qident.eq.'Hug_E0') Hug_E0=vhulp
c$$$      if (qident.eq.'shea_v') xImpVcm=vhulp
c$$$      if (qident.eq.'shok_t') ishock_type=int(vhulp)
c$$$      if (qident.eq.'shok_z') shock_z_sep=vhulp
c$$$      if (qident.eq.'shok_v') shock_vel=vhulp
c$$$      if (qident.eq.'nreac') nreac=int(vhulp)
c$$$      if (qident.eq.'axis1') axis1=vhulp
c$$$      if (qident.eq.'axis2') axis2=vhulp
c$$$      if (qident.eq.'axis3') axis3=vhulp
c$$$      if (qident.eq.'cutof2') cutof2=vhulp
c$$$      if (qident.eq.'cutof3') cutof3=vhulp
c$$$      if (qident.eq.'mdtemp') tsetor=vhulp
c$$$      if (qident.eq.'mdtem2') tset2=vhulp
c$$$      if (qident.eq.'mdpres') pset=vhulp*0.001
c$$$      if (qident.eq.'tincr') tincr=vhulp
c$$$      if (qident.eq.'tstep') tstep=vhulp
c$$$*     if (qident.eq.'lowtap') swa=vhulp !Moved to force field
c$$$*     if (qident.eq.'uptap') swb=vhulp  !Moved to force field
c$$$      if (qident.eq.'tdamp1') taut=vhulp
c$$$      if (qident.eq.'tdamp2') taut2=vhulp
c$$$      if (qident.eq.'ntdamp') ndtau=int(vhulp)
c$$$      if (qident.eq.'pdamp1') taup=vhulp
c$$$      if (qident.eq.'tdhoov') vqnd=vhulp
c$$$      if (qident.eq.'achoov') errnh=vhulp/100.0
c$$$      if (qident.eq.'range') range=vhulp
c$$$      if (qident.eq.'nmdit') maxstp=int(vhulp)
c$$$      if (qident.eq.'nmdeqi') nequi=int(vhulp)
c$$$      if (qident.eq.'imdmet') nmethod=int(vhulp)
c$$$      if (qident.eq.'icharg') ncha=int(vhulp)
          nchaold=ncha
c$$$      if (qident.eq.'ichaen') ncha2=int(vhulp)
c$$$      if (qident.eq.'ichupd') nchaud=int(vhulp)
c$$$      if (qident.eq.'iout1') nrep1=int(vhulp)
c$$$      if (qident.eq.'iout2') nsav=int(vhulp)
c$$$      if (qident.eq.'icheck') ntest=int(vhulp)
c$$$      if (qident.eq.'ivels') nvel=int(vhulp)
c$$$      if (qident.eq.'itfix') ntscale=int(vhulp)
c$$$      if (qident.eq.'irecon') ncontrol=int(vhulp)
c$$$      if (qident.eq.'iout3') noutpt=int(vhulp)
c$$$      if (qident.eq.'iappen') napp=int(vhulp)
c$$$      if (qident.eq.'isurpr') nsurp=int(vhulp)
c$$$      if (qident.eq.'itdmet') ncons=int(vhulp)
c$$$      if (qident.eq.'iravel') nrand=int(vhulp)
c$$$      if (qident.eq.'imetho') nmm=int(vhulp)
c$$$      if (qident.eq.'endmm') endpo=vhulp
          endpoold=endpo
c$$$      if (qident.eq.'endmd') endpo2=vhulp
c$$$      if (qident.eq.'imaxmo') nfc=int(vhulp)
          nfcold=nfc
c$$$      if (qident.eq.'iout4') nsav2=int(vhulp)
c$$$      if (qident.eq.'imaxit') nmmax=int(vhulp)
          nmmaxold=nmmax
c$$$      if (qident.eq.'iout5') i5758=int(vhulp)
c$$$      if (qident.eq.'parsca') parc1=vhulp
c$$$      if (qident.eq.'parext') parc2=vhulp
c$$$      if (qident.eq.'icelop') icell=int(vhulp)
          icellold=icell
c$$$      if (qident.eq.'igeopt') ingeo=int(vhulp)
c$$$      if (qident.eq.'celopt') ccpar=vhulp
c$$$      if (qident.eq.'icelo2') icelo2=int(vhulp)
          icelo2old=icelo2
c$$$      if (qident.eq.'ideve1') nrdd=int(vhulp)
c$$$      if (qident.eq.'ideve2') nrddf=int(vhulp)
c$$$      if (qident.eq.'ibiola') nbiolab=int(vhulp)
c$$$c      if (qident.eq.'igeofo') ngeofor=int(vhulp)
c$$$      if (qident.eq.'iincop') nincrop=int(vhulp)
c$$$      if (qident.eq.'accerr') accincr=vhulp
c$$$      if (qident.eq.'iout6') nsav3=int(vhulp)
c$$$      if (qident.eq.'irten') nhop2=int(vhulp)
c$$$      if (qident.eq.'npreit') nprevrun=int(vhulp)
c$$$      if (qident.eq.'idebug') ndebug=int(vhulp)
c$$$      if (qident.eq.'volcha') volcha=vhulp
c$$$      if (qident.eq.'ixmolo') ixmolo=int(vhulp)
c$$$      if (qident.eq.'inpt') inpt=int(vhulp)
c$$$      if (qident.eq.'iconne') iconne=int(vhulp)
c$$$      if (qident.eq.'imolde') imolde=int(vhulp)
c$$$      if (qident.eq.'ianaly') ianaly=int(vhulp)
c$$$      if (qident.eq.'icentr') icentr=int(vhulp)
c$$$      if (qident.eq.'itrans') itrans=int(vhulp)
c$$$      if (qident.eq.'itrout') itrout=int(vhulp)
c$$$      if (qident.eq.'nvlist') nvlist=int(vhulp)
c$$$      if (qident.eq.'vrange') vrange=vhulp
c$$$      if (qident.eq.'vlbora') vlbora=vhulp
c$$$      if (qident.eq.'tpnrad') tpnrad=vhulp
c$$$      if (qident.eq.'ityrad') ityrad=int(vhulp)
c$$$      if (qident.eq.'iexx') iexx=int(vhulp)
c$$$      if (qident.eq.'iexy') iexy=int(vhulp)
c$$$      if (qident.eq.'iexz') iexz=int(vhulp)
c$$$      if (qident.eq.'syscha') syscha=vhulp
c$$$      if (qident.eq.'inmov1') inmov1=int(vhulp)
c$$$      if (qident.eq.'inmov2') inmov2=int(vhulp)
c$$$      if (qident.eq.'itstep') itstep=int(vhulp)
c$$$      if (qident.eq.'ifreq') ifreq=int(vhulp)
c$$$      if (qident.eq.'isymm') isymm=int(vhulp)
c$$$      if (qident.eq.'icpres') icpres=int(vhulp)
c$$$      if (qident.eq.'delvib') delvib=vhulp
c$$$      goto 10
c$$$   20 continue
      close (10)
      axis(1)=axis1
      axis(2)=axis2
      axis(3)=axis3
      if (axiss(1).gt.zero) then
      axis(1)=axiss(1)
      axis(2)=axiss(2)
      axis(3)=axiss(3)
      end if
      if (tincr.lt.0.0001.and.tincr.gt.-0.0001) tset=tsetor
      iequi=1
      if (nequi.gt.0) iequi=0
      if (iopt.eq.1.and.napp.eq.1) then
      stop 'No fort.7 and fort.8 append with iopt=1 !'
      end if
      if (mdstep.gt.0.or.nit.gt.0) nmm=nmmsav
      if (mdstep.gt.0.and.itstep.eq.1) then
      tstepmax=tstep
      tstep=tstep*(tsetor/tempmd)
      if (tstep.gt.tstepmax) tstep=tstepmax
      end if
      tstep=1.0d-15*tstep
      taus=taut
      taut=1.0d-15*taut
      taut2=1.0d-15*taut2
      taup=1.0d-15*taup
      ts2=tstep/2.0
      ts22=tstep*ts2
      return
   30 continue
      write (*,*)'Error reading control-file'
      stop 'Error reading control-file'
************************************************************************
*                                                                      *
*     Format part                                                      *
*                                                                      *
************************************************************************
 1050 format (f7.3)
 1055 format (f7.4)
 1056 format (f9.4)
 1060 format (i8)
 1070 format (f7.5)
      end
************************************************************************
************************************************************************

      subroutine staint

************************************************************************
#include "cbka.blk"
#include "cbkdcell.blk"
#include "cbkqa.blk"
#include "control.blk"
#include "small.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
      dimension bvt(nat,4)
************************************************************************
*                                                                      *
*     Generate cartesian coordinates from internal coordinate input    *
*                                                                      *
************************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In staint'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      k=0
   10 read (3,1200,end=20,err=20)(ijk(k+1,k1),k1=1,3),k2,qa(k+1),
     $bvt(k+1,3),bvt(k+1,2),bvt(k+1,1)
      qlabel(k+1)=qa(k+1)
      qresi1(k+1)='   '
      qresi2(k+1)=' '
      qresi3(k+1)='     '
      qffty(k+1)='     '
      if (k2.ne.k+1) then
      write (*,*)'Wrong order in internal coordinates at atom:',k2
      goto 20
*     stop 'Wrong order in internal coordinates'
      end if
      k=k+1
      if (k.gt.nat) then
      write (*,*)na,nat
      stop 'Maximum number of atoms exceeded'
      end if
      goto 10
   20 continue
      na=k
      
************************************************************************
*                                                                      *
*     CALCULATION OF CARTESIAN COORDINATES FROM INTERNAL COORDINAATES  *
*                                                                      *
************************************************************************

   12 C(1,1)=ZERO
      C(1,2)=ZERO
      C(1,3)=ZERO
      C(2,1)=BVT(2,1)
      C(2,2)=ZERO
      C(2,3)=ZERO
      HR=(BVT(3,2)-90.0D0)*DGRRDN
      C(3,1)=C(2,1)+BVT(3,1)*SIN(HR)
      C(3,2)=BVT(3,1)*COS(HR)
      C(3,3)=ZERO
      DO 32 K1=4,NA
      J=IJK(K1,2)
      KB=K1-1
      XH=C(J,1)
      YH=C(J,2)
      ZH=C(J,3)
      DO 13 K2=1,KB
      C(K2,1)=C(K2,1)-XH
      C(K2,2)=C(K2,2)-YH
      C(K2,3)=C(K2,3)-ZH
      DO 13 K3=1,3
   13 IF (ABS(C(K2,K3)).LT.1.0D-15) C(K2,K3)=1.0D-15
      K=IJK(K1,3)
      P2=C(K,2)*C(K,2)+C(K,3)*C(K,3)
      IF (P2.NE.ZERO) THEN
      P=SQRT(P2)
      Q=SQRT(C(K,1)*C(K,1)+P2)
      SA=C(K,2)/P
      CA=C(K,3)/P
      SB=-C(K,1)/Q
      CB=P/Q
      ELSE
      SA=ZERO
      CA=ONE
      SB=ONE
      CB=ZERO
      ENDIF
      DO 16 K2=1,KB
      AZ=C(K2,1)
      BZ=C(K2,2)
      C(K2,1)=AZ*CB+BZ*SB*SA+C(K2,3)*SB*CA
      C(K2,2)=BZ*CA-C(K2,3)*SA
   16 C(K2,3)=-AZ*SB+BZ*CB*SA+C(K2,3)*CB*CA
      IF (C(K,3).LE.ZERO) THEN
      DO 17 K2=1,KB
   17 C(K2,3)=-C(K2,3)
      ENDIF
      I=IJK(K1,1)
      IF (1.0D5*ABS(C(I,1)).LE.ABS(C(I,2))) THEN
      T1=HALF*PI
      ELSE
      YX=ABS(C(I,2)/C(I,1))
      T1=ATAN(YX)
      ENDIF
      IF (C(I,1).GE.ZERO.AND.C(I,2).LT.ZERO) T1=TWO*PI-T1
      IF (C(I,1).LT.ZERO.AND.C(I,2).GE.ZERO) T1=PI-T1
      IF (C(I,1).LT.ZERO.AND.C(I,2).LT.ZERO) T1=T1+PI
      DO 31 K2=1,KB
      IF (C(K2,1).EQ.ZERO.AND.C(K2,2).EQ.ZERO) GOTO 31
      IF (1.0D5*ABS(C(K2,1)).LT.ABS(C(K2,2))) THEN
      T2=HALF*PI
      ELSE
      YX=ABS(C(K2,2)/C(K2,1))
      T2=ATAN(YX)
      ENDIF
      IF (C(K2,1).GE.ZERO.AND.C(K2,2).LT.ZERO) T2=TWO*PI-T2
      IF (C(K2,1).LT.ZERO.AND.C(K2,2).GE.ZERO) T2=PI-T2
      IF (C(K2,1).LT.ZERO.AND.C(K2,2).LT.ZERO) T2=T2+PI
      T3=T2-T1
      IF (T3.LT.ZERO)T3=T3+TWO*PI
      RZ=SQRT(C(K2,1)*C(K2,1)+C(K2,2)*C(K2,2))
      C(K2,1)=RZ*COS(T3)
      C(K2,2)=RZ*SIN(T3)
   31 CONTINUE
      HR=(BVT(K1,2)-90.0D0)*DGRRDN
      HT=BVT(K1,3)*DGRRDN
      CHR=COS(HR)
      C(K1,1)=BVT(K1,1)*CHR*COS(HT)
      C(K1,2)=BVT(K1,1)*CHR*SIN(HT)
   32 C(K1,3)=C(IJK(K1,3),3)+BVT(K1,1)*SIN(HR)
      
      return
 1200 FORMAT(4I3,1X,A2,3F10.5,4X,I1,F10.5)
      end
************************************************************************
************************************************************************

      subroutine outint

************************************************************************
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkbo.blk"
#include "cbkconst.blk"
#include "cbkia.blk"
#include "cbkinit.blk"
#include "cbknubon2.blk"
#include "cbkqa.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "small.blk"
#include "cbksrtbon1.blk"
************************************************************************
*                                                                      *
*     Output internal coordinates                                      *
*                                                                      *
************************************************************************
      dimension dvdc(3,3),dargdc(3,3)
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In outint'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      write (91,50)qmol
      open (82,file='output.MOP',status='unknown')
      write (82,*)
      write (82,'(a40)')qmol
      write (82,*)
      close (82)

*     IF (NMOLO.GT.1) THEN
*     WRITE(6,*)' OUTPUT INTERNAL COORDINATES NOT POSSIBLE FOR ',
*    $'CALCULATION ON MORE THAN ONE MOLECULE'
*     RETURN
*     END IF

************************************************************************
*                                                                      *
*     Output of internal coordinates.                                  *
*     First 3 atoms of other input file.                               *
*                                                                      *
************************************************************************
      N1=1
      N2=2
      N3=3
C      open (82,file='output.MOP',status='unknown',access='append')
      write(91,100)N1,qa(n1)
      write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(n1),
     $zero,nzero,zero,nzero,zero,nzero,nzero,nzero,nzero
      call dista2(n1,n2,rr,dx,dy,dz)
      write(91,200)N1,N2,qa(n2),RR
      write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(n2),
     $rr,none,zero,nzero,zero,nzero,n1,nzero,nzero
      close (82)

      call dista2(n2,n3,rr,dx,dy,dz)
      hv=zero
      call calvalres(n1,n2,n3,arg,hv,dvdc,dargdc)
      WRITE(91,300)N1,N2,N3,qa(n3),rdndgr*HV,RR
C      open (82,file='output.MOP',status='unknown',access='append')
      write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(n3),
     $rr,none,rdndgr*hv,none,zero,nzero,n2,n1,nzero
      close (82)

      naih=3
      
      do i1=naih+1,na
      bomax=zero
      j1=0
      do i2=1,ia(i1,2)
      iob=ia(i1,2+i2)
      ncubo=nubon2(i1,i2)
      if (bo(ncubo).gt.bomax.and.iob.lt.i1) then
      bomax=bo(ncubo)
      j1=iob
      end if
      end do
      if (j1.eq.0) j1=i1-1
      call dista2(j1,i1,rr,dx,dy,dz)
      
      bomax=zero
      j2=0
      do i2=1,ia(j1,2)
      iob=ia(j2,2+i2)
      ncubo=nubon2(j1,i2)
      if (bo(ncubo).gt.bomax.and.iob.lt.i1.and.
     $abo(iob).gt.bo(ncubo)+0.2) then
      bomax=bo(ncubo)
      j2=iob
      end if
      end do
      if (j2.eq.0) j2=i1-2
      if (j2.eq.j1) j2=j2+1

      call calvalres(j2,j1,i1,arg,hh,dvdc,dargdc)

      bomax=zero
      j3=0
      do i2=1,ia(j2,2)
      iob=ia(j2,2+i2)
      ncubo=nubon2(j2,i2)
      if (bo(ncubo).gt.bomax.and.iob.lt.i1.and.iob.ne.j1) then
      bomax=bo(ncubo)
      j3=iob
      end if
      end do
      if (j3.eq.0) j3=i1-3
      if (j3.eq.j2.and.j3.ne.j1-1) j3=j3+1
      if (j3.eq.j2.and.j3.ne.j1-2) j3=j3+2
      if (j3.eq.j1.and.j3.ne.j2-1) j3=j3+1
      if (j3.eq.j1.and.j3.ne.j2-2) j3=j3+2

      call caltor(j3,j2,j1,i1,ht)

      write(91,400)j3,j2,j1,i1,qa(i1),ht,rdndgr*hh,rr
C      open (82,file='output.MOP',status='unknown',access='append')
      write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(i1),
     $rr,none,rdndgr*hh,none,ht,none,j1,j2,j3
      close (82)
      end do

      close(82)
      return
   50 format ('  I',2x,a60)
  100 FORMAT(9X,I3,1x,a2)
  200 FORMAT(6X,2I3,1x,a2,20X,F10.5)
  300 FORMAT(3X,3I3,1x,a2,10X,2F10.5)
  400 FORMAT(4I3,1x,a2,3F10.5)

      end
************************************************************************
************************************************************************

      subroutine outres

************************************************************************
#include "cbka.blk"
#include "cbkbo.blk"
#include "cbkch.blk"
#include "cbkd.blk"
#include "cbkenergies.blk"
#include "cbkh.blk"
#include "cbkimove.blk"
#include "cbkrbo.blk"
#include "cbktorang.blk"
#include "cbktorsion.blk"
#include "cbktregime.blk"
#include "cbkvalence.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
#include "cbkinit.blk"

************************************************************************
*                                                                      *
*     Output molecular data                                            *
*                                                                      *
************************************************************************
      dimension isort(100),iad1(100),iad2(100),iad3(100),iad4(100)
      character*60 qm2
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In outres'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      read (9,100,end=50)idata,qm2
*     if (qm2.ne.qmol) then
*     write (*,*)'Wrong molecule in outres-file'
*     write (*,*)qmol
*     write (*,*)qm2
*     return
*     end if
      do 25 i1=1,idata
      read (9,200)isort(i1),iad1(i1),iad2(i1),iad3(i1),iad4(i1)
      ndata2=ndata2+1

      if (isort(i1).eq.1) then
*     do i2=1,nbon
*     if (ib(i2,2).eq.iad1(i1).and.ib(i2,3).eq.iad2(i1)) then
*     if (iopt.ne.1) write (81,*)iad1(i1),iad2(i1),rbo(i2)
*     caldat(ndata2)=rbo(i2)
*     end if
*     end do
      call dista2(iad1(i1),iad2(i1),dish,dx,dy,dz)
      write (81,*)iad1(i1),iad2(i1),dish
      caldat(ndata2)=dish
      end if

      if (isort(i1).eq.2) then
      do i2=1,nval
      if (iv(i2,2).eq.iad1(i1).and.iv(i2,3).eq.iad2(i1).and.
     $iv(i2,4).eq.iad3(i1)) then
      if (iopt.ne.1) write (81,*)iad1(i1),iad2(i1),
     $iad3(i1),h(i2)*rdndgr
      caldat(ndata2)=h(i2)*rdndgr
      end if
      end do
      end if

      if (isort(i1).eq.3) then
      do i2=1,ntor
      if (it(i2,2).eq.iad1(i1).and.it(i2,3).eq.iad2(i1).and.
     $it(i2,4).eq.iad3(i1).and.it(i2,5).eq.iad4(i1)) then
      if (iopt.ne.1) write (81,*)iad1(i1),iad2(i1),iad3(i1),iad4(i1),
     $abs(thg(i2))
      caldat(ndata2)=abs(thg(i2))
      end if
      end do
      end if

      if (isort(i1).eq.4) then
      if (iopt.ne.1) write (81,*)estrmin
      caldat(ndata2)=estrmin
      end if

      if (isort(i1).eq.5) then
      if (iopt.ne.1) write (81,*)estrmin
      caldat(ndata2)=estrmin
      end if

      if (isort(i1).eq.6) then
      if (iopt.ne.1) write (81,*)iad1(i1),axiss(iad1(i1))
      caldat(ndata2)=axiss(iad1(i1))
      end if

      if (isort(i1).eq.7) then
      if (iopt.ne.1) write (81,*)eco
      caldat(ndata2)=eco
      end if

      if (isort(i1).eq.8) then
      do i2=1,nbon
      if (ib(i2,2).eq.iad1(i1).and.ib(i2,3).eq.iad2(i1)) then
      if (iopt.ne.1) write (81,*)iad1(i1),iad2(i1),bo(i2)
      caldat(ndata2)=bo(i2)
      end if
      end do
      end if

      if (isort(i1).eq.9) then
      if (iopt.ne.1) write (81,*)ch(iad1(i1))
      caldat(ndata2)=ch(iad1(i1))
      end if
      
      if (isort(i1).eq.10) then
      rmsg=0.0
      nmovh=0
      do i2=1,na
      do i3=1,3
      rmsg=rmsg+imove(i2)*d(i3,i2)*d(i3,i2)
      nmovh=nmovh+imove(i2)
      end do
      end do
      rmsg=sqrt(rmsg/float(nmovh*3))

      if (iopt.ne.1) write (81,*)rmsg
      caldat(ndata2)=rmsg
      end if

      if (isort(i1).eq.11) then
      if (iopt.ne.1) write (81,*)1000.0*pressu
      caldat(ndata2)=1000.0*pressu
      end if

   25 continue
      
   50 return
************************************************************************
*                                                                      *
*     Format part                                                      *
*                                                                      *
************************************************************************
  100 format (i3,a60)
  200 format (5i3)
      end
************************************************************************
************************************************************************

      subroutine readgeo

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkdistan.blk"
#include "cbkinit.blk"
#include "cbkqa.blk"
#include "cbksrtbon1.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"
      character*80 qromb
      character*25 qfileh
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readgeo'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      if (ngeofor.eq.-1) return

********************************************************************** 
*                                                                    *
*     Read in system geometry                                        *
*                                                                    *
********************************************************************** 
      if (ngeofor.eq.0) then
      call readdelphi (qfileh,iend,naold)
      namov=na
      end if
 
      if (ngeofor.eq.1) then
      call readbgf(iend,naold)
      end if
 


      if (ngeofor.eq.2) then
**********************************************************************
*                                                                    *
*     Read in free format (xmol) geometry                            *
*                                                                    *
**********************************************************************
      qr='1'
      read (3,'(i6)')na
      namov=na
      read (3,'(a60)')qmol
      do i1=1,na
      read (3,'(a80)')qromb
      ifirstchar=80
      do i2=1,80
      if (qromb(i2:i2).ne.' '.and.i2.lt.ifirstchar) ifirstchar=i2
      end do
      read (qromb(ifirstchar:80),'(a2)')qa(i1)
      read (qromb(ifirstchar+2:80),*)c(i1,1),c(i1,2),c(i1,3)
      qlabel(i1)=qa(i1)
      qresi1(i1)='   '
      qresi2(i1)=' '
      qresi3(i1)='     '
      qffty(i1)='     '
      end do
      ibity=1
      axiss(1)=-1.0
      end if
 
 
      if (ngeofor.eq.3) then
**********************************************************************
*                                                                    *
*     Read in ChemDraw CC1-file                                      *
*                                                                    *
**********************************************************************
      qr='1'
      read (3,*)na
      namov=na
      read (3,'(a60)')qmol
      do i1=1,na
      read (3,'(2x,a2,5x,3f12.6)')qa(i1),c(i1,1),c(i1,2),c(i1,3)
      end do
      end if

      if (ngeofor.eq.4) then
**********************************************************************
*                                                                    *
*     Read in .pdb-format                                            *
*                                                                    *
**********************************************************************
      qr='C'
      call readpdb(iendf)
      namov=na
      ibity=1
      axiss(1)=-1.0
      qfile(nprob)=qmol
      if (iendf.eq.1) stop 'End-of-file while reading in .pdb'
      end if

**********************************************************************
*                                                                    *
*     Set up periodic system                                         *
*                                                                    *
**********************************************************************
      axis(1)=axiss(1)
      axis(2)=axiss(2)
      axis(3)=axiss(3)
      angle(1)=angles(1)
      angle(2)=angles(2)
      angle(3)=angles(3)
      if (axiss(1).lt.zero) then
      axis(1)=axis1
      axis(2)=axis2
      axis(3)=axis3
      angle(1)=90.0
      angle(2)=90.0
      angle(3)=90.0
      end if
      halfa=angle(1)*dgrrdn
      hbeta=angle(2)*dgrrdn
      hgamma=angle(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axis(1)*sinbet*sinphi
      tm21=axis(1)*sinbet*cosphi
      tm31=axis(1)*cosbet
      tm22=axis(2)*sinalf
      tm32=axis(2)*cosalf
      tm33=axis(3)

      return
      end
************************************************************************
************************************************************************

      subroutine readdelphi (qfileh,iend,naold)

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkdcell.blk"
#include "cbkdistan.blk"
#include "cbkff.blk"
#include "cbkh.blk"
#include "cbkinit.blk"
#include "cbkqa.blk"
#include "cbkrestr.blk"
#include "cbksrtbon1.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
      character*25 qfileh
********************************************************************** 
*                                                                    *
*     Read in geometries in Delphi-format (xyz)                      *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readdelphi'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      if (imodfile.eq.1) then
      open (3,file=qfileh,status='old')
      end if
      nmmax=nmmaxold
      nfc=nfcold
      ibity=1
      iredo=1
      endpo=endpoold
      icell=icellold
      icelo2=icelo2old
      iend=0
      read (3,1000,end=900)qr,qmol 
********************************************************************** 
*                                                                    *
*     Read in restraint information (optional)                       *
*                                                                    *
********************************************************************** 
      if (qr.eq.'R'.or.qr.eq.'P'.or.qr.eq.'X') then
      qmol=qmol(7:60)
      qmolset(nuge)=qmol
      read (18,1070,end=4,err=4) nrestra
      do i1=1,nrestra
      read (18,1090)irstra(i1,1),irstra(i1,2),rrstra(i1),vkrstr(i1),
     $vkrst2(i1),rrcha(i1)
      end do
    4 continue
      end if
********************************************************************** 
*                                                                    *
*     Read in torsion restraint information (optional)               *
*                                                                    *
********************************************************************** 
      if (qr.eq.'T'.or.qr.eq.'X') then
      if (qr.eq.'T') then
      qmol=qmol(7:60)
      qmolset(nuge)=qmol
      end if
      read (28,1070,end=6,err=6) nrestrat
      do i1=1,nrestrat
      read (28,1091)irstrat(i1,1),irstrat(i1,2),irstrat(i1,3),
     $irstrat(i1,4),trstra(i1),vkrt(i1),vkr2t(i1),rtcha(i1)
      end do
    6 continue
      end if
********************************************************************** 
*                                                                    *
*     Read in valency angle restraint information (optional)         *
*                                                                    *
********************************************************************** 
      if (qr.eq.'V') then
      qmol=qmol(7:60)
      qmolset(nuge)=qmol
      read (38,1070,end=7,err=7) nrestrav
      do i1=1,nrestrav
      read (38,1092)irstrav(i1,1),irstrav(i1,2),irstrav(i1,3),
     $vrstra(i1),vkrv(i1),vkr2v(i1)
      end do
    7 continue
      end if
********************************************************************** 
*                                                                    *
*     Read in geometry                                               *
*                                                                    *
********************************************************************** 
      ibh2=0
      iequi=1
      iexco=0
      if (nequi.gt.0) iequi=0
      axiss(1)=-1.0

      if (qr.eq.'O'.or.qr.eq.'L') stop 'Not xyz-format'

      if (qr.eq.'I') then      !Delphi internal coordinates
      if (nsurp.ge.2) stop 'Int.coordinates only with 1 gemetry' 
      call staint
      goto 20
      end if

      if (qr.eq.'B') then      !Previous geometry with volume reduction
      read (3,*)
      vred=(1.0-0.01*volcha)**(0.33333)
      iexco=1
      na=naold
      do i1=1,3
      qmol=qmol
      axiss(i1)=vred*axis(i1)
      angles(i1)=angle(i1)
      do i2=1,na
      c(i2,i1)=vred*c(i2,i1)
      end do
      end do

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      ibity=2

      goto 20
      end if

      if (qr.eq.'S') then      !Previous geometry with volume expansion
      read (3,*)
      vexp=(1.0+0.01*volcha)**(0.33333)
      na=naold
      iexco=1
      do i1=1,3
      qmol=qmol
      axiss(i1)=vexp*axis(i1)
      angles(i1)=angle(i1)
      do i2=1,na
      c(i2,i1)=vexp*c(i2,i1)
      end do
      end do

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      ibity=2

      goto 20
      end if

      if (qr.eq.'F'.or.qr.eq.'Y'.or.qr.eq.'3'.or.qr.eq.'5'.
     $or.qr.eq.'P') then 
      kx=0
      ky=0
      kz=0
      ibity=2
      read(3,1005)axiss(1),axiss(2),axiss(3)
      read(3,1005)angles(1),angles(2),angles(3)

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)

      end if

      if (qr.eq.'M'.or.qr.eq.'A') then 
      nmmsav=nmm
      nmm=2 
      end if

      if (qr.eq.'A') nmm=1

      if (qr.eq.'D') then
      endpo=endpo/25 
      nmmax=nmmax*5
      qruid='HIGH PRECISION'
      end if

      if (qr.eq.'H') then
      nmmax=nmmax/10
      qruid='LOW PRECISION'
      end if

      if (qr.eq.'1'.or.qr.eq.'5') then
      nmm=1
      nmmax=1
      qruid='SINGLE POINT'
      end if

      if (qr.eq.'Y') then
      icell=0
      qruid='NO CELL OPT'
      end if

   10 read (3,1100,end=20,err=20)ir,qa(na+1),(c(na+1,i2),i2=1,3)
      qlabel(na+1)=qa(na+1)
      qresi1(na+1)='   '
      qresi2(na+1)=' '
      qresi3(na+1)='     '
      qffty(na+1)='     '
      if (ir.eq.0) goto 20
      na=na+1

      if (na.gt.nat) then
      write (*,*)'Maximum number of atom exceeded ',na,nat
      stop 'Maximum number of atoms exceeded'
      end if

      goto 10
   20 continue

      if (imodfile.eq.1) close (3)

      return
  900 iend=1
      return
 1000 format (2x,a1,1x,a60)
 1005 format (3f10.4)
 1070 format (i3)
 1090 format (2i4,2f8.4,f8.6,f10.8)
 1091 format (4i4,2f8.4,3f8.6)
 1092 format (3i4,2f8.4,2f8.6)
 1100 format (i4,1x,a2,3x,3d22.15,1x,a5,1x,i5)
      end
************************************************************************
************************************************************************

      subroutine readbgf(iendf,naold)

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkcha.blk"
#include "cbkcharmol.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkdcell.blk"
#include "cbkdistan.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkh.blk"
#include "cbkimove.blk"
#include "cbkinit.blk"
#include "cbkqa.blk"
#include "cbkrestr.blk"
#include "cbksrtbon1.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
      character*80 qromb
      character*2 qrom
      character*5 quen
      character*5 qlabhulp
      character*25 qfileh
      character*200 qhulp
********************************************************************** 
*                                                                    *
*     Read in BIOGRAF-geometry                                       *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readbgf'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      iendf=0
      ienread=0
      iredo=0
      qremark(1)=' '
      enmol=zero
      formol=zero
c$$$      if (imodfile.eq.1) then
c$$$      open (3,file=qfileh,status='old')
c$$$      end if
      open (3,file='fort.3',status='old')
      read (3,'(a40)',end=900)qromb
      ibity=0
      if (qromb(1:6).eq.'BIOGRF') ibity=1
      if (qromb(1:6).eq.'XTLGRF') ibity=2
      if (ibity.eq.0) then
      write (*,*)qromb(1:6)
      stop 'Unknown Biograf-file'
      end if
      read (qromb,'(6x,i4)')ibgfversion
      if (ibity.eq.1) qr='C'
      if (ibity.eq.2) qr='F'
      iremark=0
      iformat=0
      iline=0
      iexco=0
      iruid=1
      vvol=1.0
      nmcharge=0
      nmmax=nmmaxold
      nfc=nfcold
      ncha=nchaold
      endpo=endpoold
      icell=icellold
      icelo2=icelo2old
      axiss(1)=-1.0

   30 read (3,'(a200)',end=46,err=40)qhulp
      qstrana1(1:200)=qhulp
      iline=iline+1
      irecog=0

      if (qhulp(1:6).eq.'DESCRP') then
      read (qhulp,'(7x,a40)',end=46,err=46)qmol
      irecog=1
      end if

      if (qhulp(1:6).eq.'REMARK') then
      if (iremark.lt.20) iremark=iremark+1
      read (qhulp,'(7x,a40)',end=46,err=46)qremark(iremark)
      irecog=1
      end if

      if (qhulp(1:6).eq.'FORMAT') then
      if (iformat.lt.20) iformat=iformat+1
      read(qhulp,'(7x,a40)',end=46,err=46)qformat(iformat)
      irecog=1
      end if

      if (qhulp(1:7).eq.'VCHANGE') then
      read (qhulp(8:60),*)vvol
      vred=(1.0+(vvol-1.0))**(0.33333333)
      iexco=1
      na=naold
      qmol=qmol
      do i1=1,3
      axiss(i1)=vred*axis(i1)
      angles(i1)=angle(i1)
      do i2=1,na
      cglobal(i2,i1)=vred*cglobal(i2,i1)
      end do
      end do

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      ibity=2
      irecog=1
      end if

      if (qhulp(1:7).eq.'VCHANGX') then
      read (qhulp(8:60),*)vvol
      vred=vvol
      iexco=1
      na=naold
      qmol=qmol
      do i1=1,3
      axiss(i1)=axis(i1)
      angles(i1)=angle(i1)
      do i2=1,na
      cglobal(i2,i1)=cglobal(i2,i1)
      end do
      end do

      axiss(1)=vred*axiss(1)
      do i2=1,na
      cglobal(i2,1)=vred*cglobal(i2,1)
      end do

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      ibity=2
      irecog=1
      end if

      if (qhulp(1:7).eq.'VCHANGY') then
      read (qhulp(8:60),*)vvol
      vred=vvol
      iexco=1
      na=naold
      qmol=qmol
      do i1=1,3
      axiss(i1)=axis(i1)
      angles(i1)=angle(i1)
      do i2=1,na
      cglobal(i2,i1)=cglobal(i2,i1)
      end do
      end do

      axiss(2)=vred*axiss(2)
      do i2=1,na
      cglobal(i2,2)=vred*cglobal(i2,2)
      end do

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      ibity=2
      irecog=1
      end if

      if (qhulp(1:7).eq.'VCHANGZ') then
      read (qhulp(8:60),*)vvol
      vred=vvol
      iexco=1
      na=naold
      qmol=qmol

      do i1=1,3
      axiss(i1)=axis(i1)
      angles(i1)=angle(i1)
      do i2=1,na
      cglobal(i2,i1)=cglobal(i2,i1)
      end do
      end do

      axiss(3)=vred*axiss(3)
      do i2=1,na
      cglobal(i2,3)=vred*cglobal(i2,3)
      end do

      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      ibity=2
      irecog=1
      end if

      if (qhulp(1:6).eq.'CRYSTX') then
      read (qhulp,'(8x,6f11.5)',end=46,err=46)axiss(1),
     $axiss(2),axiss(3),angles(1),angles(2),angles(3)
      kx=0
      ky=0
      kz=0
      halfa=angles(1)*dgrrdn
      hbeta=angles(2)*dgrrdn
      hgamma=angles(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axiss(1)*sinbet*sinphi
      tm21=axiss(1)*sinbet*cosphi
      tm31=axiss(1)*cosbet
      tm22=axiss(2)*sinalf
      tm32=axiss(2)*cosalf
      tm33=axiss(3)
      kx=int(2.0*swb/tm11)
      ky=int(2.0*swb/tm22)
      kz=int(2.0*swb/tm33)
      qr='F'
      if (nmmax.eq.1.and.nmmaxold.gt.1) qr='5'
      if (icell.eq.0.and.icellold.gt.0) qr='Y'
      ibity=2
      irecog=1
      end if

      if (qhulp(1:6).eq.'PERIOD') then
      read (qhulp,'(7x,i3)',end=46,err=46)iperiod
      irecog=1
      end if

      if (qhulp(1:4).eq.'AXES') then
      read (qhulp,'(7x,a3)',end=46,err=46)qbgfaxes
      irecog=1
      end if

      if (qhulp(1:6).eq.'SGNAME') then
      read (qhulp,'(7x,a3)',end=46,err=46)qbgfsgn
      irecog=1
      end if

*     if (qhulp(1:5).eq.'CELLS') then
*     read (qhulp,'(7x,*)',end=40,err=40)kx,ky,kz
*     irecog=1
*     end if

      if (qhulp(1:6).eq.'HETATM') then
      if (ibgfversion.lt.400) then
      read (qhulp,
     $'(7x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)'
     $,end=40,err=40)
     $ir,qlabel(na+1),qresi1(na+1),qresi2(na+1),qresi3(na+1),
     $cglobal(na+1,1),cglobal(na+1,2),
     $cglobal(na+1,3),qffty(na+1),ibgr1(na+1),ibgr2(na+1),
     $chgglobal(na+1)
      else
      stop 'Unsupported Biograf-version'
      end if
      qlabhulp=qlabel(na+1)
      if (qlabhulp(1:1).eq.' ') qlabhulp=qlabhulp(2:5)
      if (qlabhulp(1:1).eq.' ') qlabhulp=qlabhulp(2:4)
      if (qlabhulp(1:1).eq.' ') qlabhulp=qlabhulp(2:3)
      if (qlabhulp(1:1).eq.'C ') qa(na+1)='C '
      if (qlabhulp(1:2).eq.'Ca') qa(na+1)='Ca'
      if (qlabhulp(1:2).eq.'Cl') qa(na+1)='Cl'
      if (qlabhulp(1:2).eq.'Cu') qa(na+1)='Cu'
      if (qlabhulp(1:2).eq.'Co') qa(na+1)='Co'
      if (qlabhulp(1:1).eq.'H ') qa(na+1)='H '
      if (qlabhulp(1:2).eq.'He') qa(na+1)='He'
      if (qlabhulp(1:1).eq.'N ') qa(na+1)='N '
      if (qlabhulp(1:2).eq.'Ni') qa(na+1)='Ni'
      if (qlabhulp(1:1).eq.'O ') qa(na+1)='O '
      if (qlabhulp(1:1).eq.'B ') qa(na+1)='B '
      if (qlabhulp(1:1).eq.'F ') qa(na+1)='F '
      if (qlabhulp(1:2).eq.'Fe') qa(na+1)='Fe'
      if (qlabhulp(1:1).eq.'P ') qa(na+1)='P '
      if (qlabhulp(1:1).eq.'S ') qa(na+1)='S '
      if (qlabhulp(1:1).eq.'Y ') qa(na+1)='Y '
      if (qlabhulp(1:2).eq.'Al ') qa(na+1)='Al'
      if (qlabhulp(1:2).eq.'Au ') qa(na+1)='Au'
      if (qlabhulp(1:2).eq.'Si') qa(na+1)='Si'
      if (qlabhulp(1:2).eq.'Pt') qa(na+1)='Pt'
      if (qlabhulp(1:2).eq.'Mo') qa(na+1)='Mo'
      if (qlabhulp(1:2).eq.'Mg') qa(na+1)='Mg'
      if (qlabhulp(1:2).eq.'Ar') qa(na+1)='Ar'
      if (qlabhulp(1:2).eq.'Zr') qa(na+1)='Zr'
      if (qlabhulp(1:2).eq.'Ti') qa(na+1)='Ti'
      if (qlabhulp(1:2).eq.'Ru') qa(na+1)='Ru'
      if (qlabhulp(1:2).eq.'Ba') qa(na+1)='Ba'
      if (qlabhulp(1:2).eq.'Bi') qa(na+1)='Bi'
      if (qlabhulp(1:2).eq.'Li') qa(na+1)='Li'
      if (qlabhulp(1:2).eq.'V ') qa(na+1)='V '
      if (qlabhulp(1:2).eq.'X ') qa(na+1)='X '
      na=na+1
      if (na.gt.nattot) then
      write (*,*)'Number of atoms:read ',na
      write (*,*)'Maximum number of atoms: ',nattot
      stop 
     $'Maximum number of atoms exceeded; increase nattot in cbka.blk'
      end if
      irecog=1
      end if

      if (qhulp(1:6).eq.'RUTYPE') then          !run-type identifiers
      irecrun=0
      read (qhulp,'(7x,a40)',end=46,err=46)qruid

      if (qruid(1:10).eq.'NORMAL RUN') then
      iruid=0
      irecrun=1
      end if

      if (qruid(1:14).eq.'HIGH PRECISION') then
      endpo=endpo/25 
      nmmax=nmmax*5
      qr='D'
      iruid=1
      irecrun=1
      end if

      if (qruid(1:13).eq.'LOW PRECISION') then
      nmmax=nmmax/10
      qr='H'
      iruid=1
      irecrun=1
      end if

      if (qruid(1:12).eq.'SINGLE POINT') then
      iruid=1
      nmmax=1
      qr='1'
      if (ibity.eq.2) qr='5'
      irecrun=1
      end if

      if (qruid(1:11).eq.'NO CELL OPT') then
      iruid=1
      icell=0
      if (ibity.eq.2) qr='Y'
      irecrun=1
      end if

      if (qruid(1:8).eq.'CELL OPT') then
      iruid=1
      icell=1
      iexco=0   !Override from VCHANGE
      read (qruid,'(8x,i6)',end=46,err=46)ncellopt
      if (ncellopt.eq.2) icell=2 !cell optimisation during energy minimisation 
      if (ncellopt.eq.3) icelo2=4  !c/a optimisation
      if (ncellopt.eq.4) icelo2=1  !only a optimisation
      if (ncellopt.eq.5) icelo2=2  !only b optimisation
      if (ncellopt.eq.6) icelo2=3  !only c optimisation
      if (ncellopt.eq.7) then
      icelo2=4  !c/a optimisation
      icell=2 !cell optimisation during energy minimisation
      end if
      if (ibity.eq.2) qr='F'
      irecrun=1
      end if

      if (qruid(1:6).eq.'MAXMOV') then
      iruid=1
      read (qruid,'(6x,i6)',end=46,err=46)nfc
      irecrun=1
      end if

      if (qruid(1:4).eq.'REDO') then
      iruid=1
      read (qruid,'(4x,i6)',end=46,err=46)iredo
      irecrun=1
      end if

      if (qruid(1:5).eq.'MAXIT') then
      iruid=1
      read (qruid,'(6x,i6)',end=46,err=46)nmmax
      if (qruid(14:18).eq.'ENDPO') then
      read (qruid,'(18x,f6.3)',end=46,err=46)endpo
      end if
      irecrun=1
      end if
      if (qruid(1:5).eq.'ENDPO') then
      iruid=1
      read (qruid,'(6x,f6.3)',end=46,err=46)endpo
      irecrun=1
      end if
      
      if (qruid(1:9).eq.'CHARGEMET') then
      iruid=1
      read (qruid,'(9x,i6)',end=46,err=46)ncha
      irecrun=1
      end if

      if (irecrun.eq.0) then
      write (*,*)'Warning: ignored RUTYPE identifier ',qruid(1:12)
      end if

      irecog=1
      end if

      if (qhulp(1:14).eq.'BOND RESTRAINT') then  
      nrestra=nrestra+1
      istart=15
      call stranal(istart,iend,vout,iout,1)
      irstra(nrestra,1)=iout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      irstra(nrestra,2)=iout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      rrstra(nrestra)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      vkrstr(nrestra)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      vkrst2(nrestra)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      rrcha(nrestra)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      itstart(nrestra)=iout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      itend(nrestra)=iout
      istart=iend
*     read (qhulp,'(15x,2i4,f8.4,f8.2,f8.5,f10.7)',end=46,err=46)
*    $irstra(nrestra,1),irstra(nrestra,2),rrstra(nrestra),
*    $vkrstr(nrestra),vkrst2(nrestra),rrcha(nrestra)
      qr='R'
      irecog=1
      end if

      if (qhulp(1:15).eq.'ANGLE RESTRAINT') then  
      nrestrav=nrestrav+1
      read (qhulp,'(16x,3i4,2f8.2,f8.4,f9.6)',end=46,err=46)
     $irstrav(nrestrav,1),irstrav(nrestrav,2),irstrav(nrestrav,3),
     $vrstra(nrestrav),vkrv(nrestrav),vkr2v(nrestrav),
     $rvcha(nrestrav)
      qr='V'
      irecog=1
      end if

      if (qhulp(1:17).eq.'TORSION RESTRAINT') then  
      nrestrat=nrestrat+1
      read (qhulp,'(18x,4i4,2f8.2,f8.4,f9.6)',end=46,err=46)
     $irstrat(nrestrat,1),irstrat(nrestrat,2),irstrat(nrestrat,3),
     $irstrat(nrestrat,4),trstra(nrestrat),vkrt(nrestrat),
     $vkr2t(nrestrat),rtcha(nrestrat)
      qr='T'
      irecog=1
      end if

      if (qhulp(1:16).eq.'MASCEN RESTRAINT') then  
      nrestram=nrestram+1
      istart=17
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      irstram(nrestram,1)=0
      if (qstrana2.eq.'x') irstram(nrestram,1)=1 
      if (qstrana2.eq.'y') irstram(nrestram,1)=2 
      if (qstrana2.eq.'z') irstram(nrestram,1)=3 
      if (qstrana2.eq.'p') irstram(nrestram,1)=4 !fixed center of mass
      if (irstram(nrestram,1).eq.0) 
     $stop 'Error in mass centre restraint'
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      irstram(nrestram,2)=iout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      irstram(nrestram,3)=iout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      rmstra1(nrestram)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      if (irstram(nrestram,1).le.3) irstram(nrestram,4)=iout
      if (irstram(nrestram,1).eq.4) rmstra2(nrestram)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      if (irstram(nrestram,1).le.3) irstram(nrestram,5)=iout
      if (irstram(nrestram,1).eq.4) rmstra3(nrestram)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      if (irstram(nrestram,1).le.3) rmstra2(nrestram)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      if (irstram(nrestram,1).le.3) rmstra3(nrestram)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      if (irstram(nrestram,1).le.3) rmcha(nrestram)=vout
      irecog=1
      end if

      if (qhulp(1:9).eq.'MOLCHARGE') then
      nmcharge=nmcharge+1
      istart=10
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      iat1mc(nmcharge)=iout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      iat2mc(nmcharge)=iout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      vmcha(nmcharge)=vout
      irecog=1
      end if
      
      if (qhulp(1:8).eq.'FIXATOMS') then
      istart=9
      call stranal(istart,iend,vout,iout,1)
      if1=iout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      if2=iout
      do i12=if1,if2
      imove(i12)=0
      end do
      irecog=1
      end if

      if (qhulp(1:11).eq.'UNIT ENERGY') then
      eenconv=zero
      read (qhulp,'(14x,a5)',end=46,err=46)quen
      if (quen.eq.'eV') eenconv=23.0408
      if (quen.eq.'EV') eenconv=23.0408
      if (quen.eq.'ev') eenconv=23.0408
      if (quen.eq.'h') eenconv=627.5
      if (quen.eq.'H') eenconv=627.5
      if (quen.eq.'kcal') eenconv=1.0
      if (quen.eq.'kCal') eenconv=1.0
      if (quen.eq.'KCAL') eenconv=1.0
      if (eenconv.eq.zero) then 
      write (*,*)quen,': unknown energy unit; assuming kcal/mol'
      eenconv=1.0
      end if
      irecog=1
      end if

      if (qhulp(1:6).eq.'ENERGY') then
      read (qhulp(7:80),*,end=46,err=46)enmol
      ienread=1
      irecog=1
      end if

      if (qhulp(1:6).eq.'GEOUPD') then
      icgeopt(nprob)=0
      icgeo=0
      irecog=1
      end if
 
      if (qhulp(1:9).eq.'NO GEOUPD') then
      icgeopt(nprob)=1
      icgeo=1
      irecog=1
      end if

      if (qhulp(1:9).eq.'FREQUENCY') then
      ifreqset(nprob)=1
      ifreq=1
      irecog=1
      end if

*     if (qhulp(1:5).eq.'FORCE') then
*     read (qhulp(6:80),*,end=46,err=46)formol
*     ienread=1
*     irecog=1
*     end if

      if (qhulp(1:6).eq.'FFIELD') goto 30 
      if (qhulp(1:6).eq.'CONECT') goto 30
      if (qhulp(1:5).eq.'ORDER') goto 30
      if (qhulp(1:1).eq.'#') goto 30
      if (qhulp(1:3).eq.'END') goto 45

      if (irecog.eq.0) then
      write (*,*)'Warning: ignored line starting with: ',qhulp(1:10)
      end if
 
      goto 30
      
   40 write (*,*)'Error on line ',iline+1,' of Biograf-input'
      stop
   45 read (3,*,err=46,end=46)
   46 continue
      if (ienread.eq.1) then 
      if (eenconv.eq.zero) then
      write (*,*)'No energy unit given; assuming kcal/mol' 
      eenconv=1.0
      end if
      enmol=enmol*eenconv                !Convert energies to kcal/mol
      end if
      
      namov=0                            !calculate number of moving atoms
      do i1=1,na
      if (imove(i1).eq.1) namov=namov+1
      end do

      if (imodfile.eq.1) close (3)
      return
  900 iendf=1
      return
      end
************************************************************************
************************************************************************

      subroutine readpdb (iendf)

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkqa.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
#include "cbksrtbon1.blk"
      character*200 qhulp
********************************************************************** 
*                                                                    *
*     Read in .pdb-geometry                                          *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readpdb'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      iendf=1
      qmol='pdb_in'
    5 read (3,'(a200)',end=10,err=900) qhulp
      qstrana1(1:200)=qhulp
      istart=1
      call stranal(istart,iend,vout,iout,1)
      istart=iend

      if (qstrana2(1:6).eq.'HEADER') then
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      qmol=qstrana2(1:20)
      end if

      if (qstrana2(1:6).eq.'HETATM'.or.qstrana2(1:4).eq.'ATOM') then
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      qa(na+1)=qstrana2(1:2)
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      c(na+1,1)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      c(na+1,2)=vout
      call stranal(istart,iend,vout,iout,1)
      istart=iend
      c(na+1,3)=vout
      na=na+1
      end if

      if (qstrana2(1:3).eq.'END'.or.qstrana2(2:4).eq.'END') then
      iendf=0
      goto 10
      end if

      goto 5
   10 continue
      return
  900 write (*,*)'Error reading in .pdb-format'
      stop 'Error reading in .pdb-format'
      end 
************************************************************************
************************************************************************

      subroutine readtreg

************************************************************************
#include "cbka.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "small.blk"
      dimension isumattreg(mtreg)
      character*200 qrom
********************************************************************** 
*                                                                    *
*     Read in temperature regime                                     *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readtreg'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      ntrc=0
      open (19,file='tregime.in',status='old',err=60)
   10 read (19,'(a200)',end=50,err=900)qrom
      qstrana1(1:200)=qrom
      if (qrom(1:1).eq.'#') goto 10
      istart=1
      ntrc=ntrc+1
      if (ntrc.gt.mtreg) then
      write (*,*)'Too many temperature regimes in tregime.in;',
     $' inrease mtreg in cbka.blk'
      stop 'Too many temperature regimes in tregime.in'
      end if
      call stranal(istart,iend,vout,iout,1)
      nittc(ntrc)=iout
      istart=iend

      if (ntrc.gt.1) then
      if (nittc(ntrc).lt.nittc(ntrc-1)) then
      ntrc=ntrc-1
      write (*,*)'Warning: wrong order or empty line in tregime.in'
      write (*,*)'Ignored lines below iteration:',nittc(ntrc)
      goto 50
      end if
      end if

      call stranal(istart,iend,vout,iout,1)
      nntreg(ntrc)=iout
      if (nntreg(ntrc).gt.mtzone) then
      write (*,*)'Too many temperature zones in tregime.in;',
     $' inrease mtzone in cbka.blk'
      stop 'Too many temperature zones in tregime.in'
      end if
      istart=iend
      isumattreg(ntrc)=0
      do i1=1,nntreg(ntrc)
      call stranal(istart,iend,vout,iout,1)
      ia1treg(ntrc,i1)=iout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      ia2treg(ntrc,i1)=iout
      istart=iend
      isumattreg(ntrc)=isumattreg(ntrc)+1+ia2treg(ntrc,i1)-
     $ia1treg(ntrc,i1)
      call stranal(istart,iend,vout,iout,1)
      tsettreg(ntrc,i1)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      tdamptreg(ntrc,i1)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      dttreg(ntrc,i1)=vout
      istart=iend
      end do
      goto 10
   50 continue
      close (19)
   60 continue
**********************************************************************
*                                                                    *
*     Check consistency temperature programs in tregime.in           *
*                                                                    *
**********************************************************************
      if (ntrc.gt.0) then
      do i1=1,ntrc
      if (isumattreg(i1).ne.na) then
      write (*,*)'Inconsistency in temperature regime nr.',i1
      write (*,*)'Number of atoms defined in tregime.in:',
     $isumattreg(i1)
      write (*,*)'Number of atoms in system:',na
      stop 'Inconsistency in tregime.in'
      end if
      end do
      end if

      return
  900 stop 'Error reading tregime.in'
      end
************************************************************************
************************************************************************

      subroutine readvreg

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkvregime.blk"
#include "control.blk"
      character*200 qrom
********************************************************************** 
*                                                                    *
*     Read in volume regime                                          *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readvreg'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      nvrc=0
      open (19,file='vregime.in',status='old',err=60)
   10 read (19,'(a200)',end=50,err=900)qrom
      qstrana1(1:200)=qrom
      if (qrom(1:1).eq.'#') goto 10
      istart=1
      nvrc=nvrc+1
      if (nvrc.gt.mvreg) then
      write (*,*)'Too many volume regimes in vregime.in;',
     $' inrease mvreg in cbka.blk'
      stop 'Too many volume regimes in vregime.in'
      end if

      call stranal(istart,iend,vout,iout,1)
      nitvc(nvrc)=iout
      istart=iend

      if (nvrc.gt.1) then
      if (nitvc(nvrc).lt.nitvc(nvrc-1)) then
      nvrc=nvrc-1
      write (*,*)'Warning: wrong order or empty line in vregime.in'
      write (*,*)'Ignored lines below iteration:',nitvc(nvrc)
      goto 50
      end if
      end if

      call stranal(istart,iend,vout,iout,1)
      nnvreg(nvrc)=iout
      if (nnvreg(nvrc).gt.mvzone) then
      write (*,*)'Too many volume regimes in vregime.in;',
     $' inrease mvzone in cbka.blk'
      stop 'Too many volume zones in vregime.in'
      end if
      istart=iend
      do i1=1,nnvreg(nvrc)
      call stranal(istart,iend,vout,iout,1)
      if (qstrana2(1:1).ne.'a'.and.qstrana2(1:1).ne.'b'.and.
     $qstrana2(1:1).ne.'c'.and.qstrana2(1:4).ne.'alfa'.and.
     $qstrana2(1:4).ne.'beta'.and.qstrana2(1:5).ne.'gamma') then
      write (*,*)qstrana2
      write (*,*)'Invalid cell parameter type in vregime.in ;',
     $' use a,b,c,alfa,beta or gamma'
      stop 'Invalid cell parameter type in vregime.in'
      end if
      qvtype(nvrc,i1)=qstrana2
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      dvvreg(nvrc,i1)=vout
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      ivsca(nvrc,i1)=1
      if (qstrana2(1:1).eq.'n') ivsca(nvrc,i1)=0
      istart=iend
      end do
      goto 10
   50 continue
      close (19)
   60 continue
      return
  900 stop 'Error reading vregime.in'
      end
************************************************************************
************************************************************************

      subroutine readereg

************************************************************************
#include "cbka.blk"
#include "cbkeregime.blk"
#include "control.blk"
      character*200 qrom
********************************************************************** 
*                                                                    *
*     Read in electric field regime                                  *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readereg'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      nerc=0
      open (19,file='eregime.in',status='old',err=60)
   10 read (19,'(a200)',end=50,err=900)qrom
      qstrana1(1:200)=qrom
      if (qrom(1:1).eq.'#') goto 10
      istart=1
      nerc=nerc+1
      if (nerc.gt.mereg) then
      write (*,*)'Too many electric field regimes in eregime.in;',
     $' inrease mereg in cbka.blk'
      stop 'Too many electric field regimes in eregime.in'
      end if
      call stranal(istart,iend,vout,iout,1)
      nitec(nerc)=iout

      if (nerc.gt.1) then
      if (nitec(nerc).lt.nitec(nerc-1)) then
      nerc=nerc-1
      write (*,*)'Warning: wrong order or empty line in eregime.in'
      write (*,*)'Ignored lines below iteration:',nitec(nerc)
      goto 50
      end if
      end if
      
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      nnereg(nerc)=iout
      if (nnereg(nerc).gt.mezone) then
      write (*,*)'Too many electric field zones in eregime.in;',
     $' inrease mezone in cbka.blk'
      stop 'Too many electric field zones in vregime.in'
      end if
      istart=iend
      do i1=1,nnereg(nerc)
      call stranal(istart,iend,vout,iout,1)
      if (qstrana2(1:1).ne.'x'.and.qstrana2(1:1).ne.'y'.and.
     $qstrana2(1:1).ne.'z') then
      write (*,*)qstrana2
      write (*,*)'Invalid field direction in eregime.in ;',
     $' use x,y or z'
      stop 'Invalid field direction in eregime.in'
      end if
      qetype(nerc,i1)=qstrana2
      istart=iend
      call stranal(istart,iend,vout,iout,1)
      ereg(nerc,i1)=vout
      istart=iend
      end do
      goto 10
   50 continue
      close (19)
   60 continue
      return
  900 stop 'Error reading vregime.in'
      end
************************************************************************
************************************************************************

      subroutine readaddmol

************************************************************************
#include "cbka.blk"
#include "cbkatomcoord.blk"
#include "cbkc.blk"
#include "cbkff.blk"
#include "cbkh.blk"
#include "control.blk"
      character*80 qromb
      character*200 qhulp
      character*5 qlabhulp
********************************************************************** 
*                                                                    *
*     Read in molecule coordinates. This molecule will be added to   *
*     the system at regular intervals                                *
*     Accepts only .bgf-format                                       *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readaddmol'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
********************************************************************** 
*                                                                    *
*     Set default values                                             *
*                                                                    *
********************************************************************** 
      iaddfreq=-1   !frequency of molecule addition; <0: no addition
      xadd=-9000.0  !x-coordinate for added molecule; <-5000.0: random
      yadd=-9000.0  !y-coordinate for added molecule; <-5000.0: random
      zadd=-9000.0  !z-coordinate for added molecule; <-5000.0: random
      iveladd=1     !1: random initial velocities; 2: read in velocities
                    !from addmol.vel
      addist=-1.00  !Minimum distance between added molecule and rest
                    !of system. < 0.0: do not check
      nadattempt=10  !Number of attempts at adding the molecule
      taddmol=-1.0  !Temperature added molecule. <0.0: system temperature
      open (19,file='addmol.bgf',status='old',err=60)
      read (19,'(a40)',end=900,err=900)qromb
      if (qromb(1:6).ne.'BIOGRF') then
      write (*,*)'addmol.bgf should start with BIOGRF'
      stop 'addmol.bgf should start with BIOGRF'
      end if
      naa=0
      iline=0
   30 read (19,'(a200)',end=900,err=900)qhulp
      irecog=0
      iline=iline+1

      if (qhulp(1:6).eq.'DESCRP') then
      irecog=1
      end if
 
      if (qhulp(1:6).eq.'FORMAT') then
      irecog=1
      end if

      if (qhulp(1:6).eq.'REMARK') then
      irecog=1
      end if

      if (qhulp(1:6).eq.'HETATM') then
      irecog=1
      read (qhulp,'(7x,i5,1x,a5,1x,3x,1x,1x,1x,5x,3f10.5)'
     $,end=900,err=900)
     $ir,qlabhulp,cadd(naa+1,1),cadd(naa+1,2),cadd(naa+1,3)
      if (qlabhulp(1:1).eq.' ') qlabhulp=qlabhulp(2:5)
      if (qlabhulp(1:1).eq.' ') qlabhulp=qlabhulp(2:4)
      if (qlabhulp(1:1).eq.' ') qlabhulp=qlabhulp(2:3)
      if (qlabhulp(1:1).eq.'C ') qadd(naa+1)='C '
      if (qlabhulp(1:2).eq.'Ca') qadd(naa+1)='Ca'
      if (qlabhulp(1:2).eq.'Cl') qadd(naa+1)='Cl'
      if (qlabhulp(1:2).eq.'Cu') qadd(naa+1)='Cu'
      if (qlabhulp(1:2).eq.'Co') qadd(naa+1)='Co'
      if (qlabhulp(1:1).eq.'H ') qadd(naa+1)='H '
      if (qlabhulp(1:2).eq.'He') qadd(naa+1)='He'
      if (qlabhulp(1:1).eq.'N ') qadd(naa+1)='N '
      if (qlabhulp(1:2).eq.'Ni') qadd(naa+1)='Ni'
      if (qlabhulp(1:1).eq.'O ') qadd(naa+1)='O '
      if (qlabhulp(1:1).eq.'B ') qadd(naa+1)='B '
      if (qlabhulp(1:1).eq.'F ') qadd(naa+1)='F '
      if (qlabhulp(1:2).eq.'Fe') qadd(naa+1)='Fe'
      if (qlabhulp(1:1).eq.'P ') qadd(naa+1)='P '
      if (qlabhulp(1:1).eq.'S ') qadd(naa+1)='S '
      if (qlabhulp(1:1).eq.'Y ') qadd(naa+1)='Y '
      if (qlabhulp(1:2).eq.'Al') qadd(naa+1)='Al'
      if (qlabhulp(1:2).eq.'Au') qadd(naa+1)='Au'
      if (qlabhulp(1:2).eq.'Si') qadd(naa+1)='Si'
      if (qlabhulp(1:2).eq.'Pt') qadd(naa+1)='Pt'
      if (qlabhulp(1:2).eq.'Mo') qadd(naa+1)='Mo'
      if (qlabhulp(1:2).eq.'Mg') qadd(naa+1)='Mg'
      if (qlabhulp(1:2).eq.'Ar') qadd(naa+1)='Ar'
      if (qlabhulp(1:2).eq.'Zr') qadd(naa+1)='Zr'
      if (qlabhulp(1:2).eq.'Ba') qadd(naa+1)='Ba'
      if (qlabhulp(1:2).eq.'X ') qadd(naa+1)='X '
      ityadd(naa+1)=0
      do i1=1,nso  !Find force field type
      if (qadd(naa+1).eq.qas(i1)) ityadd(naa+1)=i1 
      end do
      if (ityadd(naa+1).eq.0) then
      write (*,*) 'Unknown atom type:',qadd(naa+1)
      stop 'Unknown atom type'
      end if
      naa=naa+1
      end if

      if (qhulp(1:7).eq.'FREQADD') then
      irecog=1
      read (qhulp,'(8x,i6)',end=900,err=900) iaddfreq
      end if

      if (qhulp(1:6).eq.'VELADD') then
      irecog=1
      read (qhulp,'(8x,i6)',end=900,err=900) iveladd
      end if

      if (qhulp(1:6).eq.'STARTX') then
      irecog=1
      read (qhulp,'(7x,f8.2)',end=900,err=900) xadd
      end if

      if (qhulp(1:6).eq.'STARTY') then
      irecog=1
      read (qhulp,'(7x,f8.2)',end=900,err=900) yadd
      end if

      if (qhulp(1:6).eq.'STARTZ') then
      irecog=1
      read (qhulp,'(7x,f8.2)',end=900,err=900) zadd
      end if

      if (qhulp(1:6).eq.'ADDIST') then
      irecog=1
      read (qhulp,'(7x,f8.2)',end=900,err=900) addist
      end if

      if (qhulp(1:8).eq.'NATTEMPT') then
      irecog=1
      read (qhulp,'(9x,i6)',end=900,err=900) nadattempt
      end if

      if (qhulp(1:7).eq.'TADDMOL') then
      irecog=1
      read (qhulp,'(8x,f8.2)',end=900,err=900) taddmol
      end if

      if (qhulp(1:6).eq.'FFIELD') goto 30
      if (qhulp(1:6).eq.'CONECT') goto 30
      if (qhulp(1:5).eq.'ORDER') goto 30
      if (qhulp(1:1).eq.'#') goto 30
      if (qhulp(1:3).eq.'END') goto 45

      if (irecog.eq.0) then
      write (*,*)'Warning: ignored line starting with: ',qhulp(1:10)
      end if

      goto 30

   45 continue
      close (19)
      if (iveladd.eq.2) then
      open (19,file='addmol.vel',status='old',err=800)
      read (19,*)
      read (19,'(3d24.15)',err=850,end=850)
     $((veladd(j,i),j=1,3),i=1,naa)       
      close (19)
      end if
************************************************************************
*                                                                      *
*     Place molecule at origin                                         *
*                                                                      *
************************************************************************
      ccx=0.0
      ccy=0.0
      ccz=0.0
      do i1=1,naa
      ccx=ccx+cadd(i1,1)/float(naa)
      ccy=ccy+cadd(i1,2)/float(naa)
      ccz=ccz+cadd(i1,3)/float(naa)
      end do
      do i1=1,naa
      cadd(i1,1)=cadd(i1,1)-ccx
      cadd(i1,2)=cadd(i1,2)-ccy
      cadd(i1,3)=cadd(i1,3)-ccz
      end do

   60 continue
      return
  800 stop 'Error opening addmol.vel'
  850 stop 'Error or end of file reading addmol.vel'
  900 write (*,*)'Error or end-of-file reading addmol.bgf on line:',
     $iline
      return
      end
************************************************************************
**********************************************************************

      subroutine writegeo(nunit1)

**********************************************************************   
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkqa.blk"
#include "cbkrestr.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
#include "cbksrtbon1.blk"
#include "cbkinit.blk"
**********************************************************************   
*                                                                    *
*     Copy new geometries to unit nunit1                             *
*                                                                    *
**********************************************************************   
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In writegeo'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      if (axiss(1).lt.zero) then
      if (nrestra.eq.0.and.nrestrat.eq.0.and.
     $nrestrav.eq.0) 
     $write (nunit1,300)qr,qmol
      if (nrestra.gt.0) write (nunit1,301)qr,
     $rrstra(1),qmol
      if (nrestrav.gt.0) write (nunit1,301)qr,
     $vrstra(1),qmol
      if (nrestrat.gt.0) write (nunit1,301)qr,
     $trstra(1),qmol
      else
      write (nunit1,310)qr,qmol
      write (nunit1,320)axiss(1),axiss(2),axiss(3)
      write (nunit1,320)angles(1),angles(2),angles(3)
      end if
      do i1=1,na
      if (nbiolab.ne.1) write (nunit1,400)i1,qa(i1),(c(i1,i2),i2=1,3)
      if (nbiolab.eq.1) write (nunit1,401)i1,qa(i1),(c(i1,i2),i2=1,3) !Delphi-format
      end do
      if (nbiolab.ne.1) write (nunit1,*)
     
      return

  300 format (2x,a1,1x,a60)
  301 format (2x,a1,1x,f6.2,a60)
  310 format (2x,a1,1x,a60)
  320 format (3f10.4)
  400 format (i4,1x,a2,3x,3(d21.14,1x),1x,a5,1x,i5)
  401 format (i3,2x,a2,3x,3(d21.14,1x),1x,a5,1x,i5)
      end
**********************************************************************
**********************************************************************

      subroutine writebgf(nunit1)

**********************************************************************   
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkcha.blk"
#include "cbkcharmol.blk"
#include "cbkconst.blk"
#include "cbkenergies.blk"
#include "cbkia.blk"
#include "cbkimove.blk"
#include "cbkinit.blk"
#include "cbkqa.blk"
#include "cbkrestr.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "cbksrtbon1.blk"
#include "small.blk"

      dimension qdir(3) 
      character*2 qt
      character*1 qdir
**********************************************************************   
*                                                                    *
*     Copy new Biograf-geometries to unit nunit1                     *
*                                                                    *
**********************************************************************   
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In newbgf'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      irom=1
      qdir(1)='x'
      qdir(2)='y'
      qdir(3)='z'
      ibgfversion=200
      if (ibity.eq.1) write (nunit1,1500)ibgfversion
      if (ibity.eq.2) write (nunit1,1600)ibgfversion
*     if (qr.ne.'F'.and.qr.ne.'5'.and.qr.ne.'Y') 
*    $write (nunit1,1500)ibgfversion
*     if (qr.eq.'F'.or.qr.eq.'5'.or.qr.eq.'Y') 
*    $write (nunit1,1600)ibgfversion
      write (nunit1,1700)qmol
*     write (nunit1,1700)qkeyw(nprob)
      do i1=1,iremark
      write (nunit1,1800)qremark(i1)
      end do
      qruid='NORMAL RUN'
      if (iruid.eq.0) then
      write (nunit1,2000)
      else
      if (abs(endpo-endpoold).gt.1e-5) write (nunit1,2010)endpo
      if (nmmax.ne.nmmaxold) write (nunit1,2020)nmmax
      if (nfc.ne.nfcold) write (nunit1,2030)nfc
      if (ncha.ne.nchaold) write (nunit1,2036)ncha
      if (iredo.gt.1) write (nunit1,2035)iredo
      if (icell.ne.icellold) then
      if (icell.eq.0) write (nunit1,2033)
      if (icell.gt.0) write (nunit1,2034)ncellopt
      end if
      end if
      if (iexco.ne.0.and.nsurp.gt.0) then
      write (nunit1,2040)vvol
      write (nunit1,3500)
      write (nunit1,*)
      return
      end if
      if (nmcharge.gt.0) then
      do i3=1,nmcharge
      write (nunit1,2050)iat1mc(i3),iat2mc(i3),vmcha(i3)
      end do
      end if

      ims=0
      do i1=1,na
      if (ims.eq.0.and.imove(i1).eq.0) then
      if1=i1
      ims=1
      end if
      if (ims.eq.1.and.imove(i1).eq.1) then
      write (nunit1,2060)if1,i1-1
      ims=0
      end if
      end do
      if (ims.eq.1) then
      write (nunit1,2060)if1,na
      end if

*     if (qr.eq.'F'.or.qr.eq.'5'.or.qr.eq.'Y') 
      if (ibity.eq.2)
     $write (nunit1,2100)axiss(1),axiss(2),axiss(3),angles(1),
     $angles(2),angles(3)

      if (nrestra.gt.0) write (nunit1,2300)
      do i2=1,nrestra
      write (nunit1,2400)
     $irstra(i2,1),irstra(i2,2),rrstra(i2),
     $vkrstr(i2),vkrst2(i2),rrcha(i2),itstart(i2),itend(i2)
      end do

      if (nrestrav.gt.0) write (nunit1,2500)
      do i2=1,nrestrav
      write (nunit1,2600)
     $irstrav(i2,1),irstrav(i2,2),irstrav(i2,3),
     $vrstra(i2),vkrv(i2),vkr2v(i2),zero
      end do

      if (nrestrat.gt.0) write (nunit1,2700)
      do i2=1,nrestrat
      write (nunit1,2800)
     $irstrat(i2,1),irstrat(i2,2),irstrat(i2,3),
     $irstrat(i2,4),trstra(i2),vkrt(i2),
     $vkr2t(i2),zero
      end do

      if (nrestram.gt.0) write (nunit1,2810)
      do i2=1,nrestram
      write (nunit1,2820)
     $qdir(irstram(i2,1)),irstram(i2,2),irstram(i2,3),
     $rmstra1(i2),irstram(i2,4),irstram(i2,5),rmstra2(i2),
     $rmstra3(i2),rmcha(i2)
      end do

      if (icgeo.eq.0.and.ingeo.eq.0) write (nunit1,2830)
      if (icgeo.eq.1.and.ingeo.eq.1) write (nunit1,2840)
      if (ifreq.eq.1) write (nunit1,2850)
      write (nunit1,2900)
      do i2=1,na
      write (nunit1,3000)i2,qa(i2),c(i2,1),c(i2,2),c(i2,3),
     $qa(i2),irom,irom,chgbgf(i2)
      end do
      write (nunit1,3100)
      if (nsurp.lt.2) then
      do i1=1,na
      write (nunit1,3200)i1,(iag(i1,2+i2),i2=1,iag(i1,2))
      end do
      write (nunit1,3300)
      write (nunit1,3400)estrc
      end if

      write (nunit1,3500)
      write (nunit1,*)
     
      return
 1500 format ('BIOGRF',i4)
 1600 format ('XTLGRF',i4)
 1700 format ('DESCRP ',a60)
 1800 format ('REMARK ',a60)
 1900 format ('FFIELD ',a40)
 2000 format ('RUTYPE NORMAL RUN')
 2010 format ('RUTYPE ENDPO',f6.3)
 2020 format ('RUTYPE MAXIT',i6)
 2030 format ('RUTYPE MAXMOV',i6)
 2033 format ('RUTYPE NO CELL OPT')
 2034 format ('RUTYPE CELL OPT',i6)
 2035 format ('RUTYPE REDO',i6)
 2036 format ('RUTYPE CHARGEMET',i6)
 2040 format ('VCHANGE',f8.4)
 2050 format ('MOLCHARGE',2i4,f6.2)
 2060 format ('FIXATOMS',2i6)    
 2100 format ('CRYSTX ',6f11.5)
 2200 format ('CELLS ',6i5)
 2300 format ('#              At1 At2   R12    Force1  Force2  ',
     $'dR12/dIter(MD) Start (MD) End (MD)')
 2400 format ('BOND RESTRAINT ',2i4,f8.4,f8.2,f8.4,1x,f10.7,2i8)
 2500 format ('#               At1 At2 At3 Angle   Force1  Force2',
     $'  dAngle/dIteration (MD only)')
 2600 format ('ANGLE RESTRAINT ',3i4,2f8.2,f8.4,f9.6)
 2700 format ('#                 At1 At2 At3 At3 Angle   Force1  ',
     $'Force2  dAngle/dIteration (MD only)')
 2800 format ('TORSION RESTRAINT ',4i4,2f8.2,f8.4,f9.6)
 2810 format ('#              x/y/z At1 At2    R   At3 At4 Force1',
     $'  Force2  dR/dIteration (MD only)')
 2820 format ('MASCEN RESTRAINT ',a1,1x,2i4,f8.2,2i4,2f8.2,f9.6)
 2830 format ('GEOUPD')
 2840 format ('NO GEOUPD')
 2850 format ('FREQUENCY')
 2900 format ('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,',
     $'3f10.5,1x,a5,i3,i2,1x,f8.5)')
 3000 format ('HETATM',1x,i5,1x,a2,3x,1x,3x,1x,1x,1x,5x,3f10.5,1x,
     $a5,i3,i2,1x,f8.5)
 3100 format ('FORMAT CONECT (a6,12i6)')
 3200 format ('CONECT',12i6)
 3300 format ('UNIT ENERGY   kcal')
 3400 format ('ENERGY',5x,f14.6)
 3500 format ('END')

      end
**********************************************************************
**********************************************************************

      subroutine writeen(tottime,sum1,sdev,sdeva,sum12,sumt,sump,
     $sumtt,tmax,eaver,eav2,eav3,etot2,ediff)
**********************************************************************   
#include "cbka.blk"
#include "cbkcha.blk"
#include "cbkenergies.blk"
#include "cbkrestr.blk"
#include "cbktorang.blk"
#include "cbktorsion.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "small.blk"

      dimension disres(mrestra)
**********************************************************************
*                                                                    *
*     Write out MD statistics to units 71,73 and 76                  *
*                                                                    *
**********************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In writeen'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      if (nrep1.gt.1)
     $sdev=sqrt((sum12-sum1*sum1/float(nrep1))/float(nrep1-1))
      eavn=eaver/float(mdstep)
      if (mdstep.gt.1)
     $sdeva=sqrt((eav3-eav2*eav2/float(mdstep))/float(mdstep-1))
C      open (71,file='fort.71',status='unknown',access='append')
C      open (73,file='fort.73',status='unknown',access='append')
      write (71,'(i8,2i4,1x,19(f10.2,1x))')mdstep+nprevrun,nmolo,
     $nmolo5,estrc,ekin,estrc+ekin,tempmd,sum1/float(nrep1),eavn,
     $sumt/float(nrep1),tmax,sump/float(nrep1),sdev,sdeva,tset,
     $tstep*1d+15,rmsg,tottime
      write (73,'(i8,1x,14(f10.2,1x))')mdstep+nprevrun,eb,ea,elp,
     $emol,ev,ecoa,ehb,et,eco,ew,ep,ech,efi
      close (71)
      close (73)
 
      if ((sumt/float(nrep1)).gt.tset) then
      if (invt.eq.0) write (*,*)'Switched to NVT in iteration',mdstep
      invt=1
      end if
 
C      if (nrestra.gt.0.or.nrestrat.gt.0)
C     $open (76,file='fort.76',status='unknown',access='append')
 
      if (nrestra.gt.0) then
      do i2=1,nrestra
      call dista2(irstra(i2,1),irstra(i2,2),disres(i2),dx,dy,dz)
      end do
C      open (76,file='fort.76',status='unknown',access='append')
      write (76,'(i8,1x,40f12.4)')mdstep,eres,estrc,
     $(rrstra(i2),disres(i2),i2=1,nrestra)
      end if
 
      if (nrestrat.gt.0) then
C      open (76,file='fort.76',status='unknown',access='append')
      do i2=1,nrestrat
      do i3=1,ntor
      ih1=irstrat(i2,1)
      ih2=irstrat(i2,2)
      ih3=irstrat(i2,3)
      ih4=irstrat(i2,4)
      if (ih1.eq.it(i3,2).and.ih2.eq.it(i3,3).and.ih3.eq.it(i3,4)
     $.and.ih4.eq.it(i3,5)) ittr=i3
      end do
      write (76,'(i8,1x,40f12.4)')mdstep,eres,
     $trstra(i2),thg(ittr)
      end do
      end if
 
      if (nrestra.gt.0.or.nrestrat.gt.0) close(76)
 
      if (nrestram.gt.0) then
C      open (76,file='fort.76',status='unknown',access='append')
      do i2=1,nrestram
      write (76,'(2i8,1x,20f12.4)')mdstep,i2,eres,rmstra1(i2),
     $dismacen(i2)
      end do
      close (76)
      end if

      return
      end
**********************************************************************
************************************************************************

      subroutine molanal

************************************************************************
#include "cbka.blk"
#include "cbkbo.blk"
#include "cbkconst.blk"
#include "cbkdcell.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkrbo.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
#include "cbksrtbon1.blk"
      dimension iam(nat,mbond+3),nmolata(nmolmax,nat)
      dimension molfra(nmolmax,nsort),ndup(nmolmax)
      character*40 qmolan1
      character*100 qmolan
      logical found
************************************************************************
*                                                                      *
*     Analyse and output molecular fragments                           *
*                                                                      *
************************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In molanal'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      do i1=1,nmolmax
      do i2=1,nsort
      molfra(i1,i2)=0
      end do
      ndup(i1)=1
      end do

      do i1=1,na
      do i2=1,mbond+3
      iam(i1,i2)=0
      end do
      end do
************************************************************************
*                                                                      *
*     Create connection table based on corrected bond orders           *
*                                                                      *
************************************************************************
      do i1=1,nbon
      if (bo(i1).gt.cutof3) then
      j1=ib(i1,2)
      j2=ib(i1,3)
      iam(j1,2)=iam(j1,2)+1
      iam(j1,2+iam(j1,2))=j2
      iam(j2,2)=iam(j2,2)+1
      iam(j2,2+iam(j2,2))=j1
      end if
      end do
**********************************************************************
*                                                                    *
*     Find molecules                                                 *
*                                                                    *
**********************************************************************
      nmolo6=0
      found=.FALSE.
      DO 61 k1=1,na
      IF (iam(K1,3+mbond).EQ.0) found=.TRUE.
   61 IF (iam(K1,3+mbond).GT.nmolo6) nmolo6=iam(K1,3+mbond)
      IF (.NOT.FOUND) GOTO 62
************************************************************************
*                                                                      *
*     Molecule numbers are assigned. No restrictions are made for the  *
*     sequence of the numbers in the connection table.                 *
*                                                                      *
************************************************************************
      N3=1
   64 N2=N3
      nmolo6=nmolo6+1
      if (nmolo6.gt.nmolmax) stop 'Too many molecules in system'
      iam(N2,3+mbond)=nmolo6
   67 FOUND=.FALSE.
      DO 66 N1=N2+1,na
      IF (iam(N1,3+mbond).NE.0) GOTO 66
      DO 65 L=1,mbond
      IF (iam(N1,l+2).EQ.0) GOTO 66
      IF (iam(iam(N1,l+2),3+mbond).EQ.nmolo6) THEN
      FOUND=.TRUE.
      iam(N1,3+mbond)=nmolo6
      GOTO 66
      ENDIF
   65 CONTINUE
   66 CONTINUE
      IF (FOUND) GOTO 67
      DO 63 N3=N2+1,NA
   63 if (iam(N3,3+mbond).eq.0) goto 64
************************************************************************
*                                                                      *
*     The assigned or input molecule numbers are checked for their     *
*     consistency.                                                     *
*                                                                      *
************************************************************************
   62 FOUND=.FALSE.
      DO 72 N1=1,NA
      DO 71 L=1,mbond
      IF (iam(N1,L+2).EQ.0) GOTO 72
      IF (iam(iam(N1,L+2),3+mbond).NE.iam(N1,3+mbond)) THEN
      FOUND=.TRUE.
      ENDIF
   71 CONTINUE
   72 CONTINUE
      IF (FOUND) THEN
      write (7,'(i4,a40)')na,qmol
      do i1=1,na
      write (7,'(40i4)')i1,iam(i1,1),(iam(i1,2+i2),i2=1,nsbmax),
     $iam(i1,3+mbond)
      end do
      STOP' Mol.nrs. not consistent; maybe wrong cell parameters'
      ENDIF

      do i1=1,nmolo6
      natmol=0
      do i2=1,na
      if (iam(i2,3+mbond).eq.i1) then
      natmol=natmol+1
      nmolata(i1,natmol+1)=i2
      end if
      end do
      nmolata(i1,1)=natmol
      end do
************************************************************************
*                                                                      *
*     Analyze molecules                                                *
*                                                                      *
************************************************************************
      do i1=1,nmolo6
      do i2=1,nmolata(i1,1)
      i3=nmolata(i1,1+i2)
      ityp=ia(i3,1)
      molfra(i1,ityp)=molfra(i1,ityp)+1
      end do
      end do

      do i1=1,nmolo6
      isee=0
      do i2=1,nmolo6
      isee2=1
      do i3=1,nso
      if (molfra(i1,i3).ne.molfra(i2,i3)) isee2=0
      end do
      if (isee2.eq.1.and.i1.gt.i2.and.isee.eq.0) then  !molecule type already exists
      ndup(i2)=ndup(i2)+1
      ndup(i1)=0
      isee=1
      end if

      end do
      end do

C      open (45,file='molfra.out',status='unknown',access='append')
      if (mdstep.eq.0) write (45,100)cutof3
      write (45,110)
      ntotmol=0
      ntotat=0
      vtotmass=zero
      do i1=1,nmolo6
      if (ndup(i1).gt.0) then
*     write (45,110)i1,(molfra(i1,i2),i2=1,nso),ndup(i1)
      ntotmol=ntotmol+ndup(i1)
      qmolan=' '
      qmolan1=' '
      istart=-4
      ihulp=0
      vmass=zero
      do i2=1,nso
      vmass=vmass+molfra(i1,i2)*amas(i2)
      ntotat=ntotat+molfra(i1,i2)*ndup(i1)
      if (molfra(i1,i2).gt.0) then
      istart=istart+6
      iend=istart+5
      if (molfra(i1,i2).gt.1) then
      write (qmolan(istart:iend),'(a2,i3)')qas(i2),molfra(i1,i2)
      else
      write (qmolan(istart:iend-2),'(a2)')qas(i2)
      end if
      end if
      end do
      ihulp=1
      do i2=1,iend
      if (qmolan(i2:i2).ne.' ') then
      qmolan1(ihulp:ihulp)=qmolan(i2:i2)
      ihulp=ihulp+1
      end if
      end do

*     write (45,120)ndup(i1),qmolan(1:iend),vmass
      write (45,120)mdstep,ndup(i1),qmolan1,vmass
      vtotmass=vtotmass+ndup(i1)*vmass
      end if
      end do
      write (45,*)'Total number of molecules:',ntotmol
      write (45,*)'Total number of atoms:',ntotat
      write (45,*)'Total system mass:',vtotmass
      close (45)
      return
  100 format('Bond order cutoff:',f6.4)
  110 format('Iteration Freq. Molecular formula',15x,'Molecular mass')
  120 format(i8,i4,' x  ',a35,f10.4)
      end
************************************************************************
************************************************************************

      subroutine stranal(istart,iend,vout,iout,icheck)

************************************************************************
#include "cbka.blk"
#include "cbkconst.blk"
#include "opt.blk"

      character*1 qchar
      dimension qchar(5)
**********************************************************************
*                                                                    *
*     Analyze string for special characters; find words in string    *
*                                                                    *
**********************************************************************
      qchar(1)=' '
      qchar(2)='/'
      
      ifound1=0
      do i1=istart,200
      ifound2=0
      do i2=1,icheck

      if (qstrana1(i1:i1).eq.qchar(i2)) then
      ifound2=1
      if (ifound1.eq.1) then     !End of word
      iend=i1
      goto 10
      end if

      end if

      end do

      if (ifound2.eq.0.and.ifound1.eq.0) then     !Start of word
      istart2=i1
      ifound1=1
      end if

      end do

   10 continue
      qstrana2=' '
      vout=zero
      iout=0

      if (ifound1.eq.1) then
      qstrana2=qstrana1(istart2:iend-1)
      istart=istart2
      vout=zero
      read (qstrana2,*,end=20,err=20) vout
   20 iout=int(vout)
      end if

      return
      end
************************************************************************
********************************************************************** 

      subroutine dipmom(naold,dpmm,xdip,ydip,zdip,xdir,ydir,zdir)

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkch.blk"
#include "cbkconst.blk"
#include "control.blk"
#include "small.blk"
********************************************************************** 
*                                                                    *
*     Calculate and output dipole moment                             *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In dipmom'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
************************************************************************
*                                                                      *
*     CONVERSION FACTOR TO DEBYE UNITS IS CALCULATED                   *
*     THE CALCULATION IS INITIALIZED                                   *
*                                                                      *
************************************************************************
 
      ELCHG=1.60217733D-19       ! [C]       = [As]
      CLIGHT=2.99792458D8        !             [m/s]
      DBCONV=ONE/(CLIGHT*ELCHG*1.0D11)
 
      CHCPX=ZERO
      CHCPY=ZERO
      CHCPZ=ZERO
      CHCMX=ZERO
      CHCMY=ZERO
      CHCMZ=ZERO
      XDIP=ZERO
      YDIP=ZERO
      ZDIP=ZERO
      XGRD=ZERO
      YGRD=ZERO
      ZGRD=ZERO
************************************************************************
*                                                                      *
*     CALCULATION OF MAGNITUDE AND CENTRES OF POSITIVE AND NEGATIVE    *
*     CHARGES                                                          *
*                                                                      *
************************************************************************
 
      if (na.eq.0) na=naold
      CHRG=ZERO
      DO 4 K1=1,NA
      CHK1=CH(K1)
      IF (CHK1.EQ.ZERO) GOTO 4
      IF (CHK1.LT.ZERO) GOTO 3
      CHRG=CHRG+CHK1
      CHCPX=CHCPX+CHK1*C(K1,1)
      CHCPY=CHCPY+CHK1*C(K1,2)
      CHCPZ=CHCPZ+CHK1*C(K1,3)
      GOTO 4
    3 CHCMX=CHCMX-CHK1*C(K1,1)
      CHCMY=CHCMY-CHK1*C(K1,2)
      CHCMZ=CHCMZ-CHK1*C(K1,3)
    4 CONTINUE
 
************************************************************************
*                                                                      *
*     CALCULATION OF DISTANCE BETWEEN CENTRES AND OF DIPOLE MOMENT     *
*     IN DEBIJE UNITS                                                  *
*                                                                      *
************************************************************************
 
      CHDSTX=CHCPX-CHCMX
      CHDSTY=CHCPY-CHCMY
      CHDSTZ=CHCPZ-CHCMZ
      DPMM=SQRT(CHDSTX*CHDSTX+CHDSTY*CHDSTY+CHDSTZ*CHDSTZ)/DBCONV
      IF(DPMM.LT.1.0D-4)RETURN
      XDIP=HALF*(CHCPX+CHCMX)/CHRG
      YDIP=HALF*(CHCPY+CHCMY)/CHRG
      ZDIP=HALF*(CHCPZ+CHCMZ)/CHRG
      GRTST=MAX(CHDSTX,CHDSTY,CHDSTZ)
      XDIR=-CHDSTX/GRTST
      YDIR=-CHDSTY/GRTST
      ZDIR=-CHDSTZ/GRTST
      open (64,file='dipole.out',status='unknown')
      write (64,100)dpmm,xdip,ydip,zdip,xdir,ydir,zdir
      close (64)
 
  100 format ('Dipole moment (Debye):',f12.4,' Location:',3f12.4,
     $' Direction (-side):',3f12.4) 
      return
      end
************************************************************************
**********************************************************************

      subroutine readtraj(ivels)

**********************************************************************
#include "cbka.blk"
#include "cbkatomcoord.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkdistan.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"
#include "cbkinit.blk"
**********************************************************************
*                                                                    *
*     Read in trajectory file                                        *
*                                                                    *
**********************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In readtraj'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      open(unit=66,file='moldyn.vel',status='old',err=10)
      ivels=1
      read (66,*)
      read (66,100)aaxis,baxis,caxis
      read (66,100)angles(1),angles(2),angles(3)
      if (qr.eq.'F'.or.qr.eq.'P'.or.ngeofor.eq.1) then
      axis(1)=aaxis
      axis(2)=baxis
      axis(3)=caxis
      axiss(1)=axis(1)
      axiss(2)=axis(2)
      axiss(3)=axis(3)
      angle(1)=angles(1)
      angle(2)=angles(2)
      angle(3)=angles(3)
      halfa=angle(1)*dgrrdn
      hbeta=angle(2)*dgrrdn
      hgamma=angle(3)*dgrrdn
      sinalf=sin(halfa)
      cosalf=cos(halfa)
      sinbet=sin(hbeta)
      cosbet=cos(hbeta)
      cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
      if (cosphi.gt.1.0) cosphi=1.0
      sinphi=sqrt(one-cosphi*cosphi)
      tm11=axis(1)*sinbet*sinphi
      tm21=axis(1)*sinbet*cosphi
      tm31=axis(1)*cosbet
      tm22=axis(2)*sinalf
      tm32=axis(2)*cosalf
      tm33=axis(3)
      end if
      if (aaxis.ne.axis(1).or.baxis.ne.axis(2).or.caxis.ne.axis(3))
     $stop 'Wrong cell parameters in moldyn.vel'
      read (66,200)nan
      if (nan.ne.na) stop 'Wrong number of atoms in moldyn.vel-file'
      if (nbiolab.eq.1) write (*,*)'Warning: using labels in vels-file'
      read (66,250)((c(i,j),j=1,3),qlabel(i),i=1,na)
      read (66,*)
      read (66,300)((vel(j,i),j=1,3),i=1,na)
      read (66,*)
      read (66,300)((accel(j,i),j=1,3),i=1,na)
      read (66,*)
      read (66,300,end=10,err=10)((aold(j,i),j=1,3),i=1,na)
      read (66,*)
      read (66,300,end=10,err=10)tempmd
      read (66,*)
      read (66,350,end=10,err=10)nsbma2
   10 continue
**********************************************************************
*                                                                    *
*     Format part                                                    *
*                                                                    *
**********************************************************************
  100 format(3d15.8)
  200 format(i4)
  250 format(3d24.15,1x,a5)
  300 format(3d24.15)
  350 format(i3)
  400 format (8i3,8f8.4)
      return
      end
**********************************************************************

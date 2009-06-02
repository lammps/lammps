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

      subroutine encalc

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkcha.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkdcell.blk"
#include "cbkenergies.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"
********************************************************************** 
*                                                                    *
*     Calculate energy and first derivatives                         *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In encalc'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      estrc=0.0
      do i1=1,na
      do i2=1,3
      d(i2,i1)=0.0
      estrain(i1)=0.0
      end do
      end do
      eb=zero
      ea=zero
      elp=zero
      emol=zero
      ev=zero
      ehb=zero
      ecoa=zero
      epen=zero
      et=zero
      eco=zero
      eres=zero
      eradbo=zero
      efi=zero

      if(Lvirial.eq.1) then
         do k1 = 1,6
            virial(k1) = zero
         end do
      endif

      if (Latomvirial.eq.1) then
      do i1=1,na
      do i2=1,6
      atomvirial(i2,i1)=0.0
      end do
      end do
      endif

      call boncor
      call lonpar
      call covbon
      call ovcor

      call srtang   !Determine valency angles
      call srttor   !Determine torsion angles
*     call srtoop   !Determine out of plane angles
      call srthb    !Determine hydrogen bonds
  
      call calval
      call valang
  
*     call oopang

      call torang
      call hbond
      !print *, 'called hbond'
      !print *, nchaud

      call nonbon
      call efield

      call restraint

c
c     Use this to print out fort.73-style energies
c     It only works correctly in serial mode
c
c      write (6,'(i8,1x,14(f21.10,1x))')mdstep+nprevrun,eb,ea,elp,
c     $emol,ev+epen,ecoa,ehb,et,eco,ew,ep,ech,efi

      estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+ncha2*ech+efi
      if (estrc.gt.zero) return
      if (estrc.le.zero) then
      goto 10
      else
      write (*,*)mdstep
      write (92,*)eb,ea,elp,ev,ecoa,emol,epen,eoop,et,eco,ew,
     $ep,ech,eres,eradbo
      stop 'Energy not a number'
      end if

   10 continue
      return
      end
********************************************************************** 
********************************************************************** 

      subroutine reaxinit

********************************************************************** 
#include "cbka.blk"
#include "cbkatomcoord.blk"
#include "cbkc.blk"
#include "cbkcha.blk"
#include "cbkconst.blk"
#include "cbkdcell.blk"
#include "cbkdistan.blk"
#include "cbkenergies.blk"
#include "cbkia.blk"
#include "cbkimove.blk"
#include "cbkinit.blk"
#include "cbktregime.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In init'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
********************************************************************** 
*                                                                    *
*     Initialize variables                                           *
*                                                                    *
********************************************************************** 
      convmd=4.184*1.0d26
      pi=3.14159265
      avognr=6.0221367d23
      rdndgr=180.0/pi
      dgrrdn=1.0/rdndgr
      rgasc=8.314510
      caljou=4.184
      xjouca=1.0/caljou
      ech=zero
      zero=0.0
      one=1.0
      two=2.0
      three=3.0
      half=one/two
      nzero=0
      none=1
      ntwo=2
      nthree=3
      invt=0
      ndata2=0
      iheatf=0
      nradcount=0
      itemp=1
      xinh=zero
      ifieldx=0
      ifieldy=0
      ifieldz=0
      mdstep=0
      kx=0
      ky=0
      kz=0
      nit=0
      nbon=0
      angle(1)=90.0
      angle(2)=90.0
      angle(3)=90.0
      axiss(1)=zero
      axiss(2)=zero
      axiss(3)=zero
      do i1=1,nat
      id(i1,1)=0
      id(i1,2)=0
      id(i1,3)=0
      end do
      icgeo=0
      sumhe=zero
      ustime=zero
      systime=zero
      vpmax=zero
      vpmin=zero
      dseed=0
      iagain=0
      do i1=1,nat
      do i2=1,mbond+3
      ia(i1,i2)=0
      iag(i1,i2)=0
      end do
      end do

      ioldchg=0
      na=0
      nrestra=0
      nrestrav=0
      nrestrat=0
      nrestram=0
      tset=tsetor
      tm11=axis(1)
      tm21=zero
      tm31=zero
      tm22=axis(2)
      tm32=zero
      tm33=axis(3)
      qruid='NORMAL RUN'
c$$$      do i1=1,nat
c$$$      imove(i1)=1
c$$$      end do

********************************************************************** 
*                                                                    *
*     Write file headers                                             *
*                                                                    *
********************************************************************** 
Cc      open (71,file='fort.71',status='unknown',access='append')
c      write (71,10)
c      close (71)
Cc      open (73,file='fort.73',status='unknown',access='append')
c      write (73,20)
c      close (73)
c      if (ntrc.gt.0) then
Cc      open (75,file='fort.75',status='unknown',access='append')
c      write (75,30)
c      close (75)
c      end if
c      if (nmethod.eq.4) then
Cc      open (59,file='fort.59',status='unknown',access='append')
c      write (59,40)
c      close (59)
c      end if

      return
********************************************************************** 
*                                                                    *
*     Format part                                                    *
*                                                                    *
********************************************************************** 
   10 format ('   Iter. Nmol    Epot         Ekin      Etot ',
     $'       T(K)  Eaver(block) Eaver(total) Taver      Tmax   ',
     $' Pres(MPa)   sdev(Epot)  sdev(Eaver)    Tset      Timestep',
     $'    RMSG     Totaltime')
   20 format ('  Iter.      Ebond       Eatom       Elp        Emol',
     $'       Eval       Ecoa       Ehbo       Etors      Econj',
     $'      Evdw      Ecoul    Echarge  Efield')
   30 format ('   Iter.  Tsys    Tzone1  Tset1   Tzone2  Tset2')
   40 format ('    Iter.     a             b          c        px',
     $'(MPa)    py(MPa)      pz(MPa)     pset(MPa)  Volume ')
      end
********************************************************************** 
************************************************************************
 
c      subroutine timer(nunit)
 
************************************************************************
c#include "cbka.blk"
c#include "cbkinit.blk"
c      real timear
c      real tarray(2)
c#ifdef _IBM
c      call dtime_(tarray,timear)
c#else
c      call dtime(tarray,timear)
c#endif

c      ustime=ustime+tarray(1)
c      systime=systime+tarray(2)
c      write (nunit,100)ustime,systime,ustime+systime
c      return
c  100 format ('User time:',f20.4,' System time:',f20.4,
c     $' Total time:',f20.4)
c      end
************************************************************************
************************************************************************

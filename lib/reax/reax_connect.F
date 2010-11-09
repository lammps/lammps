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

      subroutine srtatom

********************************************************************** 
#include "cbka.blk"
#include "cbkatomcoord.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkqa.blk"
#include "control.blk"
#include "opt.blk"
#include "small.blk"
********************************************************************** 
*                                                                    *
*     Determine atom types in system                                 *
*                                                                    *
********************************************************************** 
*     Requires the following variables
*     ndebug - opt.blk; determines whether to debug or not; everywhere
*     xmasmd - cbka.blk; some sort of atmoic mass?; srtatom, reac.f
*     molin - opt.blk; keeps info on?; srtatom
*     nso - cbka.blk; number of atoms?; srtatom, inout.f
*     nprob - cbka.blk; does?; connect.f, inout.f, reac.f
*     nasort - cbka.blk; a sorting array; srtatom
*     ia - cbka.blk; atom numbers?; poten.f, inout.f, connect.f, charges.f
*     iag - cbka.blk; ; connect.f, inout.f, poten.f, reac.f
*     xmasat - cbka.blk; does?; srtatom, reac.f
*     amas - cbka.blk; ? ; srtatom, ffinpt, molanal, ovcor
*     qa - cbka.blk; some sort of error statement variable?; srtatom, srtbon1, inout.f, radbo
*     
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In srtatom'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      xmasmd=0.0
      do i1=1,nso
      molin(nprob,i1)=0
      nasort(i1)=0
      end do
      do i1=1,na
      ia(i1,1)=0
      iag(i1,1)=0
      do i2=1,nso
      if (qa(i1).eq.qas(i2)) then
      ia(i1,1)=i2
      iag(i1,1)=i2
      molin(nprob,i2)=molin(nprob,i2)+1
      xmasat(i1)=amas(i2)
      xmasmd=xmasmd+amas(i2)
      nasort(i2)=nasort(i2)+1
      end if
      end do
      if (ia(i1,1).eq.0) then
      write (*,*)'Unknown atom type: ',qa(i1)
      stop 'Unknown atom type'
      end if
      end do

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine molec

********************************************************************** 
#include "cbka.blk"
#include "cbkdcell.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkmolec.blk"
#include "cbknmolat.blk"
#include "control.blk"
#include "small.blk"
      dimension nmolo2(nat),iseen(nmolmax),isee2(nmolmax)
********************************************************************** 
*                                                                    *
*     Determine changes in molecules                                 *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In molec'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      npreac=0

      do i1=1,nmolo
      natmol=0
      do i2=1,na
      if (ia(i2,3+mbond).eq.i1) then
      natmol=natmol+1
      nmolat(i1,natmol+1)=i2
      end if
      end do
      nmolat(i1,1)=natmol
      end do

      if (nmolo5.lt.nmolo5o) nradcount=0     !reset reaction counter
      do i1=1,nmolo5
      natmol=0
      do i2=1,na
      if (iag(i2,3+mbond).eq.i1) then
      natmol=natmol+1
      nmolat2(i1,natmol+1)=i2
      end if
      end do
      nmolat2(i1,1)=natmol
      end do
      nmolo5o=nmolo5

      do i1=nmolo+1,nmoloold
      do i2=1,nmolat(i1,1)
      nmolat(i1,1+i2)=0
      end do
      nmolat(i1,1)=0
      end do

      do i1=1,nmolo
      elmol(i1)=0.0
      do i2=1,nmolat(i1,1)
      ihu=nmolat(i1,i2+1)
      ity=ia(ihu,1)
      elmol(i1)=elmol(i1)+stlp(ity)
      end do
      end do
 
      do i1=1,nmolo5
      elmol2(i1)=0.0
      do i2=1,nmolat2(i1,1)
      ihu=nmolat2(i1,i2+1)
      ity=iag(ihu,1)
      elmol2(i1)=elmol2(i1)+stlp(ity)
      end do
      end do
 
      return
      end
********************************************************************** 
********************************************************************** 

      subroutine dista2 (n1,n2,dista,dx,dy,dz)

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
********************************************************************** 
*                                                                    *
*     Determine interatomic distances                                *
*                                                                    *
********************************************************************** 
c$$$*     if (ndebug.eq.1) then
c$$$C*     open (65,file='fort.65',status='unknown',access='append')
c$$$*     write (65,*) 'In dista2'
c$$$*     call timer(65)
c$$$*     close (65)
c$$$*     end if
      
      dx=c(n1,1)-c(n2,1)
      dy=c(n1,2)-c(n2,2)
      dz=c(n1,3)-c(n2,3)
      dista=sqrt(dx*dx+dy*dy+dz*dz)

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine srtbon1(lprune,lhb,hbcut_in,lhbnew_in,ltripstaball_in)

********************************************************************** 
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkbo.blk"
#include "cbkbosi.blk"
#include "cbkbopi.blk"
#include "cbkbopi2.blk"
#include "cbkc.blk"
#include "cbkch.blk"
#include "cbkconst.blk"
#include "cbkdbopidc.blk"
#include "cbkdrdc.blk"
#include "cbkia.blk"
#include "cbknubon2.blk"
#include "cbknvlbo.blk"
#include "cbkpairs.blk"
#include "cbknvlown.blk"
#include "cbkqa.blk"
#include "cbkrbo.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"
#include "cbkdbodc.blk"
#include "cbksrtbon1.blk"
#include "cbkff.blk"
#include "cbksrthb.blk"
#include "cbkcovbon.blk"
      logical found
      integer nboncol(nboallmax)
      integer iball(nboallmax,3)

********************************************************************** 
*                                                                    *
*     Determine connections within the molecule                      *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In srtbon1'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

c Transfer hbcut, lhbnew, and ltripstaball from C++ calling function
      hbcut = hbcut_in
      lhbnew = lhbnew_in
      ltripstaball = ltripstaball_in

      do i1=1,na
      abo(i1)=0.0d0
      end do
      nbonall=0
      nbon2=0
      nsbmax=0
      nsbma2=0

      if (imolde.eq.0) then

      nmolo=0
      nmolo5=0
      end if
      if (imolde.eq.0) then
      do i1=1,na
      do i2=2,mbond+3
      ia(i1,i2)=0
      iag(i1,i2)=0
      end do
      end do
     
      else

      do i1=1,na
      do i2=2,mbond+2
      ia(i1,i2)=0
      iag(i1,i2)=0
      end do
      end do

      end if

      do i1=1,na
      do i2=1,mbond
      nubon1(i1,i2)=0
      nubon2(i1,i2)=0
      end do
      end do

* First detect all bonds and create preliminary list
 
      do 11 ivl=1,nvpair
      if (nvlbo(ivl).eq.0) goto 11  !not in bond order range
      i1=nvl1(ivl)
      i2=nvl2(ivl)
      call dista2(i1,i2,dis,dxm,dym,dzm)
      ih1=ia(i1,1)
      ih2=ia(i2,1)
      disdx=dxm/dis
      disdy=dym/dis
      disdz=dzm/dis
      itype=0
      if (ih1.gt.ih2) then
      ih1=ia(i2,1)
      ih2=ia(i1,1)
      end if
      do i3=1,nboty2
      if (ih1.eq.nbs(i3,1).and.ih2.eq.nbs(i3,2)) itype=i3
      end do
      if (itype.eq.0.and.rat(ih1).gt.zero.and.rat(ih2).gt.zero) then
c$$$      call mdsav(1,qfile(nprob))
      write (*,*)qa(i1),'-',qa(i2),'Fatal: Unknown bond in molecule'
      stop 
      end if

      rhulp=dis/rob1(ih1,ih2)
 
********************************************************************** 
*                                                                    *
*     Determine bond orders                                          *
*                                                                    *
********************************************************************** 
      rh2=zero
      rh2p=zero
      rh2pp=zero
      ehulp=zero
      ehulpp=zero
      ehulppp=zero
      if (rapt(ih1).gt.zero.and.rapt(ih2).gt.zero) then
      rhulp2=dis/rob2(ih1,ih2)
      rh2p=rhulp2**ptp(itype)
      ehulpp=exp(pdp(itype)*rh2p)
      end if
      if (vnq(ih1).gt.zero.and.vnq(ih2).gt.zero) then
      rhulp3=dis/rob3(ih1,ih2)
      rh2pp=rhulp3**popi(itype)
      ehulppp=exp(pdo(itype)*rh2pp)
      end if

      if (rat(ih1).gt.zero.and.rat(ih2).gt.zero) then
      rh2=rhulp**bop2(itype)
      ehulp=(1.0+cutoff)*exp(bop1(itype)*rh2)
      end if

      bor=ehulp+ehulpp+ehulppp

      j1=i1
      j2=i2

********************************************************************** 
*                                                                    *
*     Determine bond orders                                          *
*                                                                    *
********************************************************************** 
      if (bor.gt.cutoff) then
      nbonall=nbonall+1
      if (nbonall.gt.nboallmax) then 
      write (6,*)'nbonall = ',nbonall,
     $ ' reax_defs.h::NBOALLMAXDEF = ',NBOALLMAXDEF,
     $     ' after',ivl, ' of ',nvpair,' pairs completed.' 
      stop 'Too many bonds; maybe wrong cell parameters.'
      end if
      iball(nbonall,1)=itype
      iball(nbonall,2)=j1
      iball(nbonall,3)=j2

      ia(i1,2)=ia(i1,2)+1
      if (ia(i1,2).gt.mbond) then 
         write (6,*)'ia(i1,2) = ',ia(i1,2),
     $        ' reax_defs.h::MBONDDEF = ',MBONDDEF,
     $        ' after',ivl, ' of ',nvpair,' pairs completed.' 
         stop 'Too many bonds on atom. Increase MBONDDEF'
      end if

      if (i1.ne.i2) then
         ia(i2,2)=ia(i2,2)+1
         if (ia(i2,2).gt.mbond) then 
            write (6,*)'ia(i1,2) = ',ia(i1,2),
     $           ' reax_defs.h::MBONDDEF = ',MBONDDEF,
     $           ' after',ivl, ' of ',nvpair,' pairs completed.' 
            stop 'Too many bonds on atom. Increase MBONDDEF'
         end if
      endif

      ia(i1,ia(i1,2)+2)=i2 
      ia(i2,ia(i2,2)+2)=i1 
      if (abs(de1(iball(nbonall,1))).gt.-0.01) then
      nubon2(i1,ia(i1,2))=nbonall
      nubon2(i2,ia(i2,2))=nbonall
      else
      nbonall=nbonall-1      !Inorganics
      end if
      end if
 11   continue

********************************************************************** 
*                                                                    *
*     lprune controls level of bond-pruning performed to increase    *
*     performance. For correct results, it should be set to 4.       *
*     However, making it smaller can speed up                        *
*     force calculation and may not have a big effect on forces.     *
*     Setting it to 0 turns off pruning, useful for debugging.       *
*                                                                    *
********************************************************************** 
********************************************************************* 
*                                                                    *
*     lhb controls whether or not to unprune ghost bonds that        *
*     may possibly form ghost hydrogen bonds.                        *
*     Setting it to 1 causes unpruning, and so is the safe option.   *
*     If lprune = 0, then pruning is not used, results are exact     *
*     and lhb has no effect.                                         *
*                                                                    *
********************************************************************** 
      if (lprune .gt. 0) then
********************************************************************** 
*                                                                    *
*     Eliminate bonds that are not in 1-6 interaction                *
*     with local atom, or closer.                                    *
*     Need additional sweep to catch possible hydrogen bonds         *
*                                                                    *
********************************************************************** 

      ntmp0 = 0
      ntmp1 = 0
      ntmp2 = 0
      ntmp3 = 0
      ntmp4 = 0
      ntmp5 = 0
      ntmp6 = 0
      ntmphb = 0

* color 1 are bonds with two local atoms
* color 2 are bonds with one local atom
* color 3 are bonds adjacent to bond with one local atom

      do i1 = 1,nbonall
         if (iball(i1,2).le.na_local) then
            if (iball(i1,3).le.na_local) then
               nboncol(i1) = 1
               ntmp1 = ntmp1+1
            else
               nboncol(i1) = 2
               ntmp2 = ntmp2+1
            endif
         else if (iball(i1,3).le.na_local) then
            nboncol(i1) = 2
            ntmp2 = ntmp2+1
         else
            nboncol(i1) = 0
         endif
      end do

      if (lprune .ge. 3) then
      do i1 = 1,nbonall
         if (nboncol(i1).eq.2) then
            if (iball(i1,2).le.na_local) then
               i3=iball(i1,3)
            else
               i3=iball(i1,2)
            endif
           
            do i4 = 1,ia(i3,2)
               i5=nubon2(i3,i4)
               if (nboncol(i5).eq.0) then
                  nboncol(i5)=3
                  ntmp3 = ntmp3+1
               endif
            end do
         endif
      end do
      endif
* color 4 bonds are part of a 1-4 interaction with local atom

      if (lprune .ge. 4) then
      do i1 = 1,nbonall
         if (nboncol(i1).eq.3) then
* One end definitely has a bond of color 2
* Find it and color bonds on other end 4
            i3=iball(i1,2)
            i3b=0
            do i4 = 1,ia(i3,2)
               i5=nubon2(i3,i4)
               if (nboncol(i5).eq.2) then
                  i3b=iball(i1,3)
               endif
            end do

            if (i3b.eq.0) then
               i3=iball(i1,3)
               i3b=0
               do i4 = 1,ia(i3,2)
                  i5=nubon2(i3,i4)
                  if (nboncol(i5).eq.2) then
                     i3b=iball(i1,2)
                  endif
               end do
            endif

            if (i3b.eq.0) then
               stop 'Could not find color 2 from color 3 bond'
            endif

            do i4 = 1,ia(i3b,2)
               i5=nubon2(i3b,i4)
               if (nboncol(i5).eq.0) then
                  nboncol(i5)=4
                  ntmp4 = ntmp4+1
               endif
            end do
           
         endif
      end do
      endif

* color 5 bonds are part of a 1-5 interaction with local atom

      if (lprune .ge. 5) then
      do i1 = 1,nbonall
         if (nboncol(i1).eq.4) then
* One end definitely has a bond of color 3
* Find it and color bonds on other end 5
            i3=iball(i1,2)
            i3b=0
            do i4 = 1,ia(i3,2)
               i5=nubon2(i3,i4)
               if (nboncol(i5).eq.3) then
                  i3b=iball(i1,3)
               endif
            end do

            if (i3b.eq.0) then
               i3=iball(i1,3)
               i3b=0
               do i4 = 1,ia(i3,2)
                  i5=nubon2(i3,i4)
                  if (nboncol(i5).eq.3) then
                     i3b=iball(i1,2)
                  endif
               end do
            endif

            if (i3b.eq.0) then
               stop 'Could not find color 3 from color 4 bond'
            endif

            do i4 = 1,ia(i3b,2)
               i5=nubon2(i3b,i4)
               if (nboncol(i5).eq.0) then
                  nboncol(i5)=5
                  ntmp5 = ntmp5+1
               endif
            end do
           
         endif
      end do
      endif

* color 6 bonds are part of a 1-6 interaction with local atom

      if (lprune .ge. 6) then
      do i1 = 1,nbonall
         if (nboncol(i1).eq.5) then
* One end definitely has a bond of color 4
* Find it and color bonds on other end 6
            i3=iball(i1,2)
            i3b=0
            do i4 = 1,ia(i3,2)
               i5=nubon2(i3,i4)
               if (nboncol(i5).eq.4) then
                  i3b=iball(i1,3)
               endif
            end do

            if (i3b.eq.0) then
               i3=iball(i1,3)
               i3b=0
               do i4 = 1,ia(i3,2)
                  i5=nubon2(i3,i4)
                  if (nboncol(i5).eq.4) then
                     i3b=iball(i1,2)
                  endif
               end do
            endif

            if (i3b.eq.0) then
               stop 'Could not find color 4 from color 5 bond'
            endif

            do i4 = 1,ia(i3b,2)
               i5=nubon2(i3b,i4)
               if (nboncol(i5).eq.0) then
                  nboncol(i5)=6
                  ntmp6 = ntmp6+1
               endif
            end do
           
         endif
      end do
      endif

* Catch all the possible hydrogen bonds
* This section replicates the logic used in srthb()
      if (lhb .eq. 1) then
c  Outer loop must be Verlet list, because ia() does not store Verlet entries,
c  but it does store bond entries in nubon2()
      do ivl=1,nvpair    !Use Verlet-list to find donor-acceptor pairs

      j1=nvl1(ivl)
      j2=nvl2(ivl)
      ihhb1=nphb(ia(j1,1))
      ihhb2=nphb(ia(j2,1))

      if (ihhb1.gt.ihhb2) then        !Make j1 donor(H) atom and j2 acceptor(O) atom
      j2=nvl1(ivl)
      j1=nvl2(ivl)
      ihhb1=nphb(ia(j1,1))
      ihhb2=nphb(ia(j2,1))
      end if

*     Only need to compute bonds where j1 is local
      if (j1 .le. na_local) then

      if (ihhb1.eq.1.and.ihhb2.eq.2) then
         call dista2(j1,j2,dishb,dxm,dym,dzm)
         if (dishb.lt.hbcut) then
            do i23=1,ia(j1,2)   !Search for acceptor atoms bound to donor atom
               if (nboncol(nubon2(j1,i23)).eq.0) then
                  j3=ia(j1,2+i23)
                  if (nphb(ia(j3,1)).eq.2.and.j3.ne.j2) then
                     nboncol(nubon2(j1,i23))=-1
                     ntmphb = ntmphb+1
                  endif
               endif
            end do
         end if
      end if
      end if
      end do
      end if

* Compact the list, removing all uncolored bonds

      nbon = 0
      do i1 = 1,nbonall
         if (nboncol(i1).eq.0) then
            ntmp0=ntmp0+1
         else
            nbon = nbon+1

            if (nbon.gt.nbomax) then 
               write (6,*)nbon,nbomax
               write (6,*)'nbon = ',nbon,' reax_defs.h::NBOMAXDEF = ',
     $              NBOMAXDEF,' after',i1, ' of ',nbonall,
     $              ' initial bonds completed.' 
               stop 'Too many pruned bonds; increase NBOMAXDEF'
            end if


            ib(nbon,1) = iball(i1,1)
            ib(nbon,2) = iball(i1,2)
            ib(nbon,3) = iball(i1,3)
         endif
      end do

********************************************************************** 
*                                                                    *
*     Do not perform ghost-bond pruning                              *
*                                                                    *
********************************************************************** 

      else

      nbon = 0
      do i1 = 1,nbonall
         nbon = nbon+1

         if (nbon.gt.nbomax) then 
            write (6,*)nbon,nbomax
            write (6,*)'nbon = ',nbon,' reax_defs.h::NBOMAXDEF = ',
     $           NBOMAXDEF,' after',i1, ' of ',nbonall,
     $           ' initial bonds completed.' 
            stop 'Too many pruned bonds; increase NBOMAXDEF'
         end if

         ib(nbon,1) = iball(i1,1)
         ib(nbon,2) = iball(i1,2)
         ib(nbon,3) = iball(i1,3)
      end do

      endif

      do i1=1,na
      do i2=2,mbond+2
      ia(i1,i2)=0
      iag(i1,i2)=0
      end do
      end do

* Generate full set of bond data structures

      do 10 i0 = 1,nbon
      i1 = ib(i0,2)
      i2 = ib(i0,3)
      call dista2(i1,i2,dis,dxm,dym,dzm)
*     do 10 i1=1,na-1
*     do 10 i2=i1+1,na
*     call dista2(i1,i2,dis,dxm,dym,dzm)
      ih1=ia(i1,1)
      ih2=ia(i2,1)
*     if (dis.gt.5.0*rob) goto 10
      disdx=dxm/dis
      disdy=dym/dis
      disdz=dzm/dis
      itype=0
      if (ih1.gt.ih2) then
      ih1=ia(i2,1)
      ih2=ia(i1,1)
      end if
      do i3=1,nboty2
      if (ih1.eq.nbs(i3,1).and.ih2.eq.nbs(i3,2)) itype=i3
      end do
      if (itype.eq.0.and.rat(ih1).gt.zero.and.rat(ih2).gt.zero) then
c$$$      call mdsav(1,qfile(nprob))
      write (*,*)qa(i1),'-',qa(i2),'Fatal: Unknown bond in molecule'
      stop 
      end if

      rhulp=dis/rob1(ih1,ih2)
 
********************************************************************** 
*                                                                    *
*     Determine bond orders                                          *
*                                                                    *
********************************************************************** 
      rh2=zero
      rh2p=zero
      rh2pp=zero
      ehulp=zero
      ehulpp=zero
      ehulppp=zero
      if (rapt(ih1).gt.zero.and.rapt(ih2).gt.zero) then
      rhulp2=dis/rob2(ih1,ih2)
      rh2p=rhulp2**ptp(itype)
      ehulpp=exp(pdp(itype)*rh2p)
      end if
      if (vnq(ih1).gt.zero.and.vnq(ih2).gt.zero) then
      rhulp3=dis/rob3(ih1,ih2)
      rh2pp=rhulp3**popi(itype)
      ehulppp=exp(pdo(itype)*rh2pp)
      end if

      if (rat(ih1).gt.zero.and.rat(ih2).gt.zero) then
      rh2=rhulp**bop2(itype)
      ehulp=(1.0+cutoff)*exp(bop1(itype)*rh2)
      end if

      bor=ehulp+ehulpp+ehulppp
      borsi=ehulp
      borpi=ehulpp
      borpi2=ehulppp
      dbordrob=bop2(itype)*bop1(itype)*rh2*(1.0/dis)*ehulp+
     $ptp(itype)*pdp(itype)*rh2p*(1.0/dis)*ehulpp+
     $popi(itype)*pdo(itype)*rh2pp*(1.0/dis)*ehulppp
      dborsidrob=bop2(itype)*bop1(itype)*rh2*(1.0/dis)*ehulp
      dborpidrob=ptp(itype)*pdp(itype)*rh2p*(1.0/dis)*ehulpp
      dborpi2drob=popi(itype)*pdo(itype)*rh2pp*(1.0/dis)*ehulppp
     
      nbon2=nbon2+1
      j1=i1
      j2=i2

********************************************************************** 
*                                                                    *
*     Determine bond orders                                          *
*                                                                    *
********************************************************************** 
      ib(i0,1)=itype
      ib(i0,2)=j1
      ib(i0,3)=j2
      ibsym(i0)=ivl
      drdc(i0,1,1)=disdx
      drdc(i0,2,1)=disdy
      drdc(i0,3,1)=disdz
      drdc(i0,1,2)=-disdx
      drdc(i0,2,2)=-disdy
      drdc(i0,3,2)=-disdz
      abo(i1)=abo(i1)+bor-cutoff
      if (i1.ne.i2) abo(i2)=abo(i2)+bor-cutoff
      bo(i0)=bor-cutoff
      bos(i0)=bor-cutoff
      bosi(i0)=borsi-cutoff
      bopi(i0)=borpi
      bopi2(i0)=borpi2
      rbo(i0)=dis
      dbodr(i0)=dbordrob
*     dbosidr(i0)=dborsidrob
      dbopidr(i0)=dborpidrob
      dbopi2dr(i0)=dborpi2drob
      dbodc(i0,1,1)=dbodr(i0)*drdc(i0,1,1)
      dbodc(i0,2,1)=dbodr(i0)*drdc(i0,2,1)
      dbodc(i0,3,1)=dbodr(i0)*drdc(i0,3,1)
      dbodc(i0,1,2)=dbodr(i0)*drdc(i0,1,2)
      dbodc(i0,2,2)=dbodr(i0)*drdc(i0,2,2)
      dbodc(i0,3,2)=dbodr(i0)*drdc(i0,3,2)
*     dbosidc(i0,1,1)=dbosidr(i0)*drdc(i0,1,1)
*     dbosidc(i0,2,1)=dbosidr(i0)*drdc(i0,2,1)
*     dbosidc(i0,3,1)=dbosidr(i0)*drdc(i0,3,1)
*     dbosidc(i0,1,2)=dbosidr(i0)*drdc(i0,1,2)
*     dbosidc(i0,2,2)=dbosidr(i0)*drdc(i0,2,2)
*     dbosidc(i0,3,2)=dbosidr(i0)*drdc(i0,3,2)
      dbopidc(i0,1,1)=dbopidr(i0)*drdc(i0,1,1)
      dbopidc(i0,2,1)=dbopidr(i0)*drdc(i0,2,1)
      dbopidc(i0,3,1)=dbopidr(i0)*drdc(i0,3,1)
      dbopidc(i0,1,2)=dbopidr(i0)*drdc(i0,1,2)
      dbopidc(i0,2,2)=dbopidr(i0)*drdc(i0,2,2)
      dbopidc(i0,3,2)=dbopidr(i0)*drdc(i0,3,2)
      dbopi2dc(i0,1,1)=dbopi2dr(i0)*drdc(i0,1,1)
      dbopi2dc(i0,2,1)=dbopi2dr(i0)*drdc(i0,2,1)
      dbopi2dc(i0,3,1)=dbopi2dr(i0)*drdc(i0,3,1)
      dbopi2dc(i0,1,2)=dbopi2dr(i0)*drdc(i0,1,2)
      dbopi2dc(i0,2,2)=dbopi2dr(i0)*drdc(i0,2,2)
      dbopi2dc(i0,3,2)=dbopi2dr(i0)*drdc(i0,3,2)
      ia(i1,2)=ia(i1,2)+1
      if (i1.ne.i2) ia(i2,2)=ia(i2,2)+1
      ia(i1,ia(i1,2)+2)=i2 
      ia(i2,ia(i2,2)+2)=i1 
      if (ia(i1,2).gt.nsbma2) nsbma2=ia(i1,2)
      if (ia(i2,2).gt.nsbma2) nsbma2=ia(i2,2)
      if (bor.gt.cutof3) then
      iag(i1,2)=iag(i1,2)+1
      iag(i2,2)=iag(i2,2)+1
      iag(i1,iag(i1,2)+2)=i2 
      iag(i2,iag(i2,2)+2)=i1 
      nubon1(i1,iag(i1,2))=i0
      nubon1(i2,iag(i2,2))=i0
      if (iag(i1,2).gt.nsbmax) nsbmax=iag(i1,2)
      if (iag(i2,2).gt.nsbmax) nsbmax=iag(i2,2)
      end if
      nubon2(i1,ia(i1,2))=i0
      nubon2(i2,ia(i2,2))=i0

   10 continue

********************************************************************** 
*                                                                    *
*     Sort molecules                                                 *
*                                                                    *
********************************************************************** 
      imolde = 1
      if (imolde.eq.1) return    !fixed molecular definitions

      FOUND=.FALSE.
      DO 31 K1=1,NA
      IF (IA(K1,3+mbond).EQ.0) FOUND=.TRUE.
   31 IF (IA(K1,3+mbond).GT.NMOLO) NMOLO=IA(K1,3+mbond)
      IF (.NOT.FOUND) GOTO 32
************************************************************************
*                                                                      *
*     Molecule numbers are assigned. No restrictions are made for the  *
*     sequence of the numbers in the connection table.                 *
*                                                                      *
************************************************************************
      N3=1
   34 N2=N3
      NMOLO=NMOLO+1
      if (nmolo.gt.nmolmax) then
      write (*,*)nmolmax
      write (*,*)'Too many molecules in system; increase nmolmax'
      write (*,*)'nmolmax = ',nmolmax
      write (*,*)'nmolo = ',nmolo
      write (*,*)'n2 = ',n2
      stop 'Too many molecules in system'
      end if
      IA(N2,3+mbond)=NMOLO
   37 FOUND=.FALSE.
      DO 36 N1=N2+1,NA
      IF (IA(N1,3+mbond).NE.0) GOTO 36
      DO 35 L=1,mbond
      IF (IA(N1,l+2).EQ.0) GOTO 36
      IF (IA(IA(N1,l+2),3+mbond).EQ.NMOLO) THEN
      FOUND=.TRUE.
      IA(N1,3+mbond)=NMOLO
      GOTO 36
      ENDIF
   35 CONTINUE
   36 CONTINUE
      IF (FOUND) GOTO 37
      DO 33 N3=N2+1,NA
   33 IF (IA(N3,3+mbond).EQ.0) GOTO 34
************************************************************************
*                                                                      *
*     The assigned or input molecule numbers are checked for their     *
*     consistency.                                                     *
*                                                                      *
************************************************************************
   32 FOUND=.FALSE.
      DO 42 N1=1,NA
      DO 41 L=1,mbond
      IF (IA(N1,L+2).EQ.0) GOTO 42
      IF (IA(IA(N1,L+2),3+mbond).NE.IA(N1,3+mbond)) THEN
      FOUND=.TRUE.
      ENDIF
   41 CONTINUE
   42 CONTINUE
      IF (FOUND) THEN
      write (7,1000)NA,qmol
      do i1=1,NA
      write (7,1100)i1,ia(i1,1),(ia(i1,2+i2),i2=1,nsbmax),
     $ia(i1,3+mbond)
      end do
      write (7,*)tm11,tm22,tm33,angle(1),angle(2),angle(3)
      STOP' Mol.nrs. not consistent; maybe wrong cell parameters'
      end if
********************************************************************** 
*                                                                    *
*     Sort molecules again                                           *
*     This sort is on iag, enforces bond order cutoff                *
*                                                                    *
********************************************************************** 
      FOUND=.FALSE.
      DO 61 K1=1,NA
      IF (IAG(K1,3+mbond).EQ.0) FOUND=.TRUE.
   61 IF (IAG(K1,3+mbond).GT.NMOLO5) NMOLO5=IAG(K1,3+mbond)
      IF (.NOT.FOUND) GOTO 62
************************************************************************
*                                                                      *
*     Molecule numbers are assigned. No restrictions are made for the  *
*     sequence of the numbers in the connection table.                 *
*                                                                      *
************************************************************************
      N3=1
   64 N2=N3
      NMOLO5=NMOLO5+1
      if (nmolo5.gt.nmolmax) stop 'Too many molecules in system'
      IAG(N2,3+mbond)=NMOLO5
   67 FOUND=.FALSE.
      DO 66 N1=N2+1,NA
      IF (IAG(N1,3+mbond).NE.0) GOTO 66
      DO 65 L=1,mbond
      IF (IAG(N1,l+2).EQ.0) GOTO 66
      IF (IAG(IAG(N1,l+2),3+mbond).EQ.NMOLO5) THEN
      FOUND=.TRUE.
      IAG(N1,3+mbond)=NMOLO5
      GOTO 66
      ENDIF
   65 CONTINUE
   66 CONTINUE
      IF (FOUND) GOTO 67
      DO 63 N3=N2+1,NA
   63 IF (IAG(N3,3+mbond).EQ.0) GOTO 64
************************************************************************
*                                                                      *
*     The assigned or input molecule numbers are checked for their     *
*     consistency.                                                     *
*                                                                      *
************************************************************************
   62 FOUND=.FALSE.
      DO 72 N1=1,NA
      DO 71 L=1,mbond
      IF (IAG(N1,L+2).EQ.0) GOTO 72
      IF (IAG(IAG(N1,L+2),3+mbond).NE.IAG(N1,3+mbond)) THEN
      FOUND=.TRUE.
      ENDIF
   71 CONTINUE
   72 CONTINUE
      IF (FOUND) THEN
      write (7,1000)NA,qmol
      do i1=1,NA
      write (7,1100)i1,iag(i1,1),(iag(i1,2+i2),i2=1,nsbmax),
     $iag(i1,3+mbond)
      end do
      write (7,*)tm11,tm22,tm33,angle(1),angle(2),angle(3)
      STOP' Mol.nrs. not consistent; maybe wrong cell parameters'
      ENDIF

********************************************************************** 
*                                                                    *
*     Format part                                                    *
*                                                                    *
********************************************************************** 
 1000 format (i3,2x,a60)
 1100 format (8i3)
      end
********************************************************************** 
********************************************************************** 

      subroutine srtang

********************************************************************** 
#include "cbka.blk"
#include "cbkbo.blk"
#include "cbknubon2.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkrbo.blk"
#include "cbkvalence.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"

      dimension a(3),b(3),j(3)
      dimension ityva(100)

********************************************************************** 
*                                                                    *
*     Find valency angles in molecule                                *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In srtang'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      nval=0
      if (nvaty.eq.0) return
      do iindexatom=1,na
      inumbonds=ia(iindexatom,2)
      do jindexbond=1,inumbonds-1
      jindexbondlist = nubon2(iindexatom,jindexbond)
      if (bo(jindexbondlist).lt.cutof2) goto 51
      k4=ib(jindexbondlist,2)
      k5=ib(jindexbondlist,3)
      do kindexbond=jindexbond+1,inumbonds
      kindexbondlist = nubon2(iindexatom,kindexbond)
      iju=0
      if (bo(kindexbondlist).lt.cutof2) goto 50
      if (bo(jindexbondlist)*bo(kindexbondlist).lt.0.001) goto 50
      k7=ib(kindexbondlist,2)
      k8=ib(kindexbondlist,3)

*     Exclude angles that have no local atoms.
*     Angles with non-local center atom are not needed for angle
*     energies, but are needed to construct torsions. 
      if ( k4 .le. na_local .or. 
     $     k5 .le. na_local .or. 
     $     k7 .le. na_local .or. 
     $     k8 .le. na_local) then 

      if (k4.eq.k7.and.k5.eq.k8.and.k4.ne.k8.and.k5.ne.k7) then
      nval=nval+1
      iv(nval,2)=k5
      iv(nval,3)=k4
      iv(nval,4)=k8
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      nval=nval+1
      iv(nval,2)=k4
      iv(nval,3)=k5
      iv(nval,4)=k7
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      iju=2
      write(6,*) 'Aaaah!'
      end if
      if (iju.eq.2) goto 50

      if (k4.eq.k8.and.k5.eq.k7.and.k4.ne.k7.and.k5.ne.k8) then
      nval=nval+1
      iv(nval,2)=k5
      iv(nval,3)=k4
      iv(nval,4)=k7
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      nval=nval+1
      iv(nval,2)=k4
      iv(nval,3)=k5
      iv(nval,4)=k8
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      iju=2
      write(6,*) 'Aaaah!'
      end if
      if (iju.eq.2) goto 50

      if (k4.eq.k7) then
      nval=nval+1
      iv(nval,2)=k5
      iv(nval,3)=k4
      iv(nval,4)=k8
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      iju=1
      end if
      if (iju.eq.1) goto 50

      if (k4.eq.k8) then
      nval=nval+1
      iv(nval,2)=k5
      iv(nval,3)=k4
      iv(nval,4)=k7
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      iju=1
      end if
      if (iju.eq.1) goto 50

      if (k5.eq.k7) then
      nval=nval+1
      iv(nval,2)=k4
      iv(nval,3)=k5
      iv(nval,4)=k8
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      iju=1
      end if
      if (iju.eq.1) goto 50

      if (k5.eq.k8) then
      nval=nval+1
      iv(nval,2)=k4
      iv(nval,3)=k5
      iv(nval,4)=k7
      iv(nval,5)=jindexbondlist
      iv(nval,6)=kindexbondlist
      iju=1
      end if
      if (iju.eq.1) goto 50

      write (6,*)'nval = ',nval,
     $        ' after',iindexatom, ' of ',na,' atoms completed.'
      stop 'Adjacent bonds did not make an angle'

      endif

   50 continue

      if (nval.gt.nvamax) then
         write (6,*)'nval = ',nval,' reax_defs.h::NVAMAXDEF = ',
     $        NVAMAXDEF,
     $        ' after',iindexatom, ' of ',na,' atoms completed.'
         stop 'Too many valency angles. Increase NVAMAXDEF'
      endif

      if (iju.gt.0) then
**********************************************************************
*                                                                    *
*     Determine force field types of angles                          *
*                                                                    *
**********************************************************************
      ityva(1)=0
      ih1=ia(iv(nval,2),1)
      ih2=ia(iv(nval,3),1)
      ih3=ia(iv(nval,4),1)
      if (ih3.lt.ih1) then
      ih3=ia(iv(nval,2),1)
      ih2=ia(iv(nval,3),1)
      ih1=ia(iv(nval,4),1)
      end if

      nfound=0
      do i3=1,nvaty
      if (ih1.eq.nvs(i3,1).and.ih2.eq.nvs(i3,2).and.
     $ih3.eq.nvs(i3,3)) then
      nfound=nfound+1
      ityva(nfound)=i3
      end if
      end do

      if (ityva(1).eq.0.or.abs(vka(ityva(1))).lt.0.001) then   !Valence angle does not exist in force field;ignore
      nval=nval-1
      ihul=0
      else
      iv(nval,1)=ityva(1)
      ihul=1

      do i3=1,nfound-1           !Found multiple angles of the same type
      nval=nval+1
      iv(nval,1)=ityva(i3+1)
      do i4=2,6
      iv(nval,i4)=iv(nval-1,i4)
      end do

      end do

      end if

      if (iju.eq.2) then
      ityva(1)=0
      ih1=ia(iv(nval-ihul,2),1)
      ih2=ia(iv(nval-ihul,3),1)
      ih3=ia(iv(nval-ihul,4),1)
      if (ih3.lt.ih1) then
      ih3=ia(iv(nval-ihul,2),1)
      ih2=ia(iv(nval-ihul,3),1)
      ih1=ia(iv(nval-ihul,4),1)
      end if

      nfound=0
      do i3=1,nvaty
      if (ih1.eq.nvs(i3,1).and.ih2.eq.nvs(i3,2).and.
     $ih3.eq.nvs(i3,3)) then
      nfound=nfound+1
      ityva(nfound)=i3
      end if
      end do

      if (ityva(1).eq.0.or.abs(vka(ityva(1))).lt.0.001) then   !Valence angle does not exist in force field;ignore
      if (ihul.eq.1) then
      do i3=1,6
      iv(nval-1,i3)=iv(nval,i3)
      end do
      end if
      nval=nval-1
      else
      iv(nval-ihul,1)=ityva(1)

      do i3=1,nfound-1           !Found multiple angles of the same type
      nval=nval+1
      iv(nval,1)=ityva(i3+1)
      do i4=2,6
      iv(nval,i4)=iv(nval-1,i4)
      end do

      end do

      end if

      end if

      end if

      end do
   51 continue
      end do
      end do

      nbonop=0
     
      return
      end
********************************************************************** 
********************************************************************** 

      subroutine srttor

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkbo.blk"
#include "cbkrbo.blk"
#include "cbkia.blk"
#include "cbktorsion.blk"
#include "cbkvalence.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"
#include "cbknubon2.blk"
********************************************************************** 
*                                                                    *
*     Find torsion angles in molecule                                *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In srttor'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      ntor=0
      if (ntoty.eq.0) return
      do 61 i1=1,nbon
      k2=ib(i1,2)
      k3=ib(i1,3)
c     Only compute interaction if both atoms local
c     are local or else flip a coin 
      if (k2 .gt. na_local) go to 61
      if (k3 .gt. na_local) then
         if (itag(k2) .lt. itag(k3)) go to 61
         if (itag(k2) .eq. itag(k3)) then
            if(c(k2,3) .gt. c(k3,3)) go to 61
            if(c(k2,3) .eq. c(k3,3) .and. 
     $           c(k2,2) .gt. c(k3,2)) go to 61
            if(c(k2,3) .eq. c(k3,3) .and. 
     $           c(k2,2) .eq. c(k3,2) .and. 
     $           c(k2,1) .gt. c(k3,1)) go to 61
         endif
      endif

      iob1=ia(k2,2)
      iob2=ia(k3,2)
      do 60 i2=1,iob1        !Atoms connected to k2
      k4=ia(k2,2+i2)
      ibo2=nubon2(k2,i2)
      do 60 i3=1,iob2        !Atoms connected to k3
      k5=ia(k3,2+i3)
      ibo3=nubon2(k3,i3)
      bopr=bo(i1)*bo(ibo2)*bo(ibo3)
      if (bopr.gt.cutof2.and.k2.ne.k5.and.k3.ne.k4.and.k4.ne.k5) then

      ntor=ntor+1
      it(ntor,2)=k4
      it(ntor,3)=k2
      it(ntor,4)=k3
      it(ntor,5)=k5
      it(ntor,6)=ibo2
      it(ntor,7)=i1
      it(ntor,8)=ibo3

**********************************************************************
*                                                                    *
*     Determine force field types of torsion angles                  *
*                                                                    *
**********************************************************************
      ity=0
      ih1=ia(it(ntor,2),1)
      ih2=ia(it(ntor,3),1)
      ih3=ia(it(ntor,4),1)
      ih4=ia(it(ntor,5),1)

      if (ih2.gt.ih3) then
      ih1=ia(it(ntor,5),1)
      ih2=ia(it(ntor,4),1)
      ih3=ia(it(ntor,3),1)
      ih4=ia(it(ntor,2),1)
      end if

      if (ih2.eq.ih3.and.ih4.lt.ih1) then
      ih1=ia(it(ntor,5),1)
      ih2=ia(it(ntor,4),1)
      ih3=ia(it(ntor,3),1)
      ih4=ia(it(ntor,2),1)
      end if

      do i4=1,ntoty
      if (ih1.eq.nts(i4,1).and.ih2.eq.nts(i4,2).and.ih3.eq.nts(i4,3)
     $.and.ih4.eq.nts(i4,4)) ity=i4
      end do

      if (ity.eq.0) then
      do i4=1,ntoty
      if (nts(i4,1).eq.0.and.ih2.eq.nts(i4,2).and.ih3.eq.nts(i4,3)
     $.and.nts(i4,4).eq.0) ity=i4
      end do
      end if

      if (ity.eq.0) then
      ntor=ntor-1           !Torsion angle does not exist in force field: ignore
      else
      it(ntor,1)=ity
      end if

      end if

   60 continue
 61   continue

      if (ntor.gt.ntomax) stop 'Too many torsion angles'
*     do i1=1,ntor
*     write (41,'(20i4)')i1,it(i1,1),it(i1,2),it(i1,3),
*    $it(i1,4),it(i1,5),it(i1,6),it(i1,7),it(i1,8)
*     end do

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine srtoop

********************************************************************** 
#include "cbka.blk"
#include "cbkbo.blk"
#include "cbkrbo.blk"
#include "cbkvalence.blk"
#include "control.blk"
#include "small.blk"
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In srtoop'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
**********************************************************************
*                                                                    *
*     Find out of plane angles in molecule                           *
*                                                                    *
**********************************************************************
      noop=0
      do i1=1,nval
      k2=iv(i1,2)
      k3=iv(i1,3)
      k4=iv(i1,4)
      k5=iv(i1,5)
      k6=iv(i1,6)
      do i2=1,nbon
      k7=ib(i2,2)
      k8=ib(i2,3)
      if (bo(i2).gt.cutof2) then
      if (k7.eq.k3.and.k8.ne.k4.and.k8.ne.k2) then
      noop=noop+1 
      ioop(noop,2)=k8
      ioop(noop,3)=k3
      ioop(noop,4)=k2
      ioop(noop,5)=k4
      ioop(noop,6)=i2
      ioop(noop,7)=iv(i1,5)
      ioop(noop,8)=iv(i1,6)
      ioop(noop,9)=i1
      end if
      if (k8.eq.k3.and.k7.ne.k4.and.k7.ne.k2) then
      noop=noop+1 
      ioop(noop,2)=k7
      ioop(noop,3)=k3
      ioop(noop,4)=k2
      ioop(noop,5)=k4
      ioop(noop,6)=i2
      ioop(noop,7)=iv(i1,5)
      ioop(noop,8)=iv(i1,6)
      ioop(noop,9)=i1
      end if
      end if
      end do
      end do
      
      do i1=1,noop
      call caltor(ioop(i1,2),ioop(i1,3),ioop(i1,4),ioop(i1,5),hoop)
      end do
      
********************************************************************** 
      return
      end
**********************************************************************

**********************************************************************

      subroutine srthb

**********************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkbo.blk"
#include "cbkconst.blk"
#include "cbkia.blk"
#include "cbkrbo.blk"
#include "cbksrthb.blk"
#include "control.blk"
#include "small.blk"
#include "cbkpairs.blk"
#include "cbknvlown.blk"
#include "cbknubon2.blk"
**********************************************************************
*                                                                    *
*     Find hydrogen bonds in molecule                                *
*                                                                    *
**********************************************************************
c$$$      if (ndebug.eq.1) then
c$$$      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In srthb'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      nhb=0
**********************************************************************
*                                                                    *
*     Locate donor/acceptor bonds                                    *
*                                                                    *
**********************************************************************
c  Outer loop must be Verlet list, because ia() does not store Verlet entries,
c  but it does store bond entries in nubon2()
c
c The problem with using the nvlown ownership criterion
c is that it would require that we unprune every bond that is within 
c certain distance, as well as its first and second neighbor bonds.
c
c For the ownership criterion based on H atom location no unpruning is required.
c Apparently lprune=4 is sufficient here, implying that we need to capture first and
c second neighbor bonds of the O-H bond, and of course we need to include all hydrogen
c bond partners within hbcut. 
c 

      do 20 ivl=1,nvpair    !Use Verlet-list to find donor-acceptor pairs

      j1=nvl1(ivl)
      j2=nvl2(ivl)
      ity1=ia(j1,1)
      ity2=ia(j2,1)
      ihhb1=nphb(ia(j1,1))
      ihhb2=nphb(ia(j2,1))

      if (ihhb1.gt.ihhb2) then        !Make j1 donor(H) atom and j2 acceptor(O) atom
      j2=nvl1(ivl)
      j1=nvl2(ivl)
      ity1=ia(j1,1)
      ity2=ia(j2,1)
      ihhb1=nphb(ia(j1,1))
      ihhb2=nphb(ia(j2,1))
      end if

* Only need to compute bonds where j1 is local
      if (j1 .le. na_local) then

      if (ihhb1.eq.1.and.ihhb2.eq.2) then
      call dista2(j1,j2,dishb,dxm,dym,dzm)
      if (dishb.lt.hbcut) then
      do 10 i23=1,ia(j1,2)                !Search for acceptor atoms bound to donor atom
      j3=ia(j1,2+i23)
      ity3=ia(j3,1)
      nbohb=nubon2(j1,i23)
      if (nphb(ity3).eq.2.and.j3.ne.j2.and.bo(nbohb).gt.0.01) then
**********************************************************************
*                                                                    *
*     Accept hydrogen bond and find hydrogen bond type               *
*                                                                    *
**********************************************************************
      nhb=nhb+1

      if (nhb.gt.nhbmax) then
      write (*,*)nhb,nhbmax
      write (*,*)'Maximum number of hydrogen bonds exceeded'
      stop 'Maximum number of hydrogen bonds exceeded'
      end if

      ihb(nhb,1)=0

      do i3=1,nhbty
      if (ity3.eq.nhbs(i3,1).and.ity1.eq.nhbs(i3,2).and.ity2.eq.
     $nhbs(i3,3)) ihb(nhb,1)=i3
      end do

      if (ihb(nhb,1).eq.0) then    !Hydrogen bond not in force field
      nhb=nhb-1 
*     write (*,*)'Warning: added hydrogen bond ',ity3,ity1,ity2
*     nhbty=nhbty+1
*     nhbs(nhbty,1)=ity3
*     nhbs(nhbty,2)=ity1
*     nhbs(nhbty,3)=ity2
*     rhb(nhbty)=2.70
*     dehb(nhbty)=zero
*     vhb1(nhbty)=5.0
*     vhb2(nhbty)=20.0
*     ihb(nhb,1)=nhbty
      end if

      ihb(nhb,2)=j3
      ihb(nhb,3)=j1
      ihb(nhb,4)=j2
      ihb(nhb,5)=nbohb
      ihb(nhb,6)=k1
      ihb(nhb,7)=k2
      ihb(nhb,8)=k3
*     write (64,*)nhb,ihb(nhb,1),j3,j1,j2,nbohb,k1,k2,k3,bo(nbohb),
*    $dishb

      end if

   10 continue

      end if
      end if
      end if
   20 end do

*     stop 'end in srthb'
      return
      end
**********************************************************************

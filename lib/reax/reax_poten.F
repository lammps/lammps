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
******************************************************************** 

      subroutine calval

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkdhdc.blk"
#include "cbkdrdc.blk"
#include "cbkh.blk"
#include "cbkrbo.blk"
#include "cbkvalence.blk"
#include "cellcoord.blk"
#include "control.blk"
      dimension a(3),b(3),j(3),dradc(3,3),drbdc(3,3),dtdc(3,3),
     $dargdc(3,3),dndc(3,3),dadc(3),dbdc(3)
********************************************************************** 
*                                                                    *
*     Calculate valency angles and their derivatives to cartesian    *
*     coordinates                                                    *
*     Valency angle energies are calculated in valang                *
*                                                                    *
********************************************************************** 
********************************************************************** 
*     Description of variables used in this routine.       
*                                                          
*     ndebug: stored in cbka.blk; control-parameter
*     third: local variable                  
*     twothird: local variable    
*     dadc(3): local array; stores derivative distance to cartesians    
*     dbdc(3): local array; stores derivative distance to cartesians
*     i1: local do-loop counter
*     i2: local do-loop counter
*     k1: local do-loop counter
*     k2: local do-loop counter
*     dradc(3,3): local array; stores derivatives bond lengths to
*                 cartesians
*     drbdc(3,3): local array; stores derivatives bond lengths to
*                 cartesians
*     nval: stored in cbka.blk; number of valence angles
*     ity: local integer; atom type                        
*     iv(nvalmax,6): stored in cbka.blk; valence angle identifiers
*     j(3): local integer array; stores valence angle atom numbers
*     la: local integer: stores bond numbers in valence angle
*     lb: local integer: stores bond numbers in valence angle
*     ivl1: local integer; stores symmetric copy number of bond
*     ivl2: local integer; stores symmetric copy number of bond
*     ibsym(nbomax): stored in cbka.blk; symmetric copy number of bond
*     isign1: local integer; -1 or 1                       
*     isign2: local integer; -1 or 1                       
*     rla: local variable; stores bond length for bond la
*     rlb: local variable; stores bond length for bond lb
*     rbo(nbomax): stored in cbka.blk; stores bond lengths
*     ix1,iy1,iz1,ix2,iy2,iz2: local integers; periodic cell shifts
*     a(3): local variable; distance in x,y and z-direction between atoms
*     b(3): local variable; distance in x,y and z-direction between atoms
*     c(nat,3): stored in cbka.blk; cartesian coordinate array
*     tm11,tm21,tm22,tm31,tm32,tm33: stored in cbka.blk; periodic cell
*                  matrix                                  
*     poem: local variable; product of bond lengths                                 
*     tel: local variable; cross-product of x,y and z-interatomic
*                  distances                    
*     arg: local variable; cosine of angle between bonds a and b
*     arg2: local variable; square of arg
*     s1ma22: local variable; used to check whether angle gets to 180
*                  degrees
*     s1ma2: local variable; square root of s1ma22                         
*     hl: local variable; angle (in radians) between bonds a and b
*     h(nvamax): stored in cbka.blk; angle (in radians) between bonds a
*                  and b
*     ib(nbomax,3): stored in cbka.blk: bond distance identifiers                   
*     drdc(nbomax,3,2): stored in cbka.blk; derivatives bond distances
*                  to cartesian coordinates
*     dndc(3,3): local variable; temporary storage for calculating
*                  derivatives of valence angle to cartesians
*     dtdc(3,3): local variable; temporary storage for calculating
*                  derivatives of valence angle to cartesians
*     dargdc(3,3): local variable; temporary storage for calculating
*                  derivatives of valence angle to cartesians
*     dhdc(nvamax,3,3): stored in cbka.blk; derivatives of valence angle 
*                  to cartesians
*                                                          
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In calval'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if

      third=1.0/3.0
      twothird=2.0/3.0
      dadc(1)=-1.0
      dadc(2)=1.0
      dadc(3)=0.0
      dbdc(1)=0.0
      dbdc(2)=1.0
      dbdc(3)=-1.0
      do k1=1,3
      do k2=1,3
      dradc(k1,k2)=0.0
      drbdc(k1,k2)=0.0
      end do
      end do
      if (nval.eq.0) return
 
      do 10 i1=1,nval
      ity=iv(i1,1)
      j(1)=iv(i1,2)
      j(2)=iv(i1,3)
      j(3)=iv(i1,4)
**********************************************************************
*                                                                    *
*     Determine valency angle                                        *
*                                                                    *
**********************************************************************
      la=iv(i1,5)
      lb=iv(i1,6)
      ivl1=ibsym(la)
      ivl2=ibsym(lb)
      isign1=1
      isign2=1
      rla=rbo(la)
      rlb=rbo(lb)
 
      call dista2(j(2),j(1),dis,a(1),a(2),a(3))
      call dista2(j(2),j(3),dis,b(1),b(2),b(3))

      poem=rla*rlb
      tel=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      arg=tel/poem
      arg2=arg*arg
      s1ma22=1.0-arg2
      if (s1ma22.lt.1.0d-10) s1ma22=1.0d-10
      s1ma2=sqrt(s1ma22)
      if (arg.gt.1.0) arg=1.0
      if (arg.lt.-1.0) arg=-1.0
      hl=acos(arg)
      h(i1)=hl
**********************************************************************
*                                                                    *
*     Calculate derivative valency angle to cartesian coordinates    *
*                                                                    *
**********************************************************************
      if (j(1).eq.ib(la,2)) then
      do k1=1,3
      dradc(k1,1)=drdc(la,k1,1)
      dradc(k1,2)=drdc(la,k1,2)
      end do
      else
      do k1=1,3
      dradc(k1,1)=drdc(la,k1,2)
      dradc(k1,2)=drdc(la,k1,1)
      end do
      end if
      if (j(2).eq.ib(lb,2)) then
      do k1=1,3
      drbdc(k1,2)=drdc(lb,k1,1)
      drbdc(k1,3)=drdc(lb,k1,2)
      end do
      else
      do k1=1,3
      drbdc(k1,2)=drdc(lb,k1,2)
      drbdc(k1,3)=drdc(lb,k1,1)
      end do
      end if
      do k1=1,3
      do k2=1,3
      dndc(k1,k2)=rla*drbdc(k1,k2)+rlb*dradc(k1,k2)
      dtdc(k1,k2)=a(k1)*dbdc(k2)+b(k1)*dadc(k2)
      dargdc(k1,k2)=(dtdc(k1,k2)-arg*dndc(k1,k2))/poem
      dhdc(i1,k1,k2)=-dargdc(k1,k2)/s1ma2
      end do
      end do
      
   10 continue

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine boncor

********************************************************************** 
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkc.blk"
#include "cbkbo.blk"
#include "cbkboncor.blk"
#include "cbkbosi.blk"
#include "cbkbopi.blk"
#include "cbkbopi2.blk"
#include "cbkconst.blk"
#include "cbkdbopi2ndc.blk"
#include "cbkdbopidc.blk"
#include "cbkdbopindc.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbknubon2.blk"
#include "cbkrbo.blk"
#include "control.blk"
#include "small.blk"
#include "cbkdbodc.blk"
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In boncor'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
********************************************************************** 
*                                                                    *
*     Correction for overcoordination and 1-3 bond orders            *
*                                                                    *
********************************************************************** 
********************************************************************** 
*     Description of variables used in this routine.       
*                                                          
*     ndebug: stored in cbka.blk; control-parameter
*     i1: local do-loop counter
*     i2: local do-loop counter
*     k1: local do-loop counter
*     k2: local do-loop counter
*     nbon: stored in cbka.blk; number of bonds in system
*     ibt: local integer; stores bond type
*     ib(nbomax,3): stored in cbka.blk: bond distance identifiers                   
*     j1: local integer; stores atom number 1st atom in bond
*     j2: local integer; stores atom number 2nd atom in bond
*     ovc(nbotym): stored in cbka.blk: force field parameter for
*                  overcoordination correction
*     v13cor(nbotym): stored in cbka.blk: force field parameter for
*                  1-3 bond order correction
*     idbo1(nbomax): stored in cbka.blk; number of atoms in the
*                  derivative of the bond order
*     idbo(nbomax,2*mbond): stored in cbka.blk; atom numbers of the
*                  atoms in the derivative of the bond order
*     dbondc(nbomax,3,2*mbond): stored in cbka.blk; derivative of
*                  corrected total bond orders to cartesians
*     dbosindc(nbomax,3,2*mbond): stored in cbka.blk; derivative of
*                  corrected sigma bond orders to cartesians
*     dbopindc(nbomax,3,2*mbond): stored in cbka.blk; derivative of
*                  corrected pi bond orders to cartesians
*     dbopi2ndc(nbomax,3,2*mbond): stored in cbka.blk; derivative of
*                  corrected double pi bond orders to cartesians
*     dbodc(nbomax,3,2): stored in cbka.blk; derivative of
*                  uncorrected total bond orders to cartesians
*     dbosidc(nbomax,3,2): stored in cbka.blk; derivative of
*                  uncorrected sigma bond orders to cartesians
*     dbopidc(nbomax,3,2): stored in cbka.blk; derivative of
*                  uncorrected pi bond orders to cartesians
*     dbopi2dc(nbomax,3,2): stored in cbka.blk; derivative of
*                  uncorrected double pi bond orders to cartesians
*     boo: local variable; storage of uncorrected total bond order
*     bo(nbomax): stored in cbka.blk; total bond order          
*     bopi(nbomax): stored in cbka.blk; pi bond order          
*     bopi2(nbomax): stored in cbka.blk; double pi bond order          
*     bopio: local variable; storage of uncorrected pi bond order
*     bopi2o: local variable; storage of uncorrected double pi bond order
*     iti: local integer; atom type first atom in bond     
*     itj: local integer; atom type second atom in bond     
*     ia(nat,mbond+3): stored in cbka.blk; connection table without bond
*                  order cutoff
*     aboi: local variable: total bond order around atom i
*     aboj: local variable: total bond order around atom j
*     abo(nat): stored in cbka.blk; total bond order around atoms
*     vp131: local variable; force field cross-term        
*     vp132: local variable; force field cross-term        
*     vp133: local variable; force field cross-term        
*     bo131(nsort): stored in cbka.blk; force field parameter for 1-3
*                  bond order correction
*     bo132(nsort): stored in cbka.blk; force field parameter for 1-3
*                  bond order correction
*     bo133(nsort): stored in cbka.blk; force field parameter for 1-3
*                  bond order correction
*     corrtot:local variable; total correction on bond order
*     dbodsboi1: local variable; derivative of bond order to sum of bond
*                  orders around atom i
*     dbodsboj1: local variable; derivative of bond order to sum of bond
*                  orders around atom j
*     ovi: local variable; overcoordination on atom i      
*     ovj: local variable; overcoordination on atom j      
*     aval(nat): stored in cbka.blk; nr. of valence electrons on atom       
*     exphu1: local variable; stores exponential
*     exphu2: local variable; stores exponential
*     exp11: local variable; stores exponential
*     exp21: local variable; stores exponential
*     vpar(npamax): stored in cbka.blk: general parameters                                  
*     exphu12: local variable; stores sum of exponential
*     ovcor: local variable; temporary storage for BO/ovcor corr.
*     huli: local variable; temporary storage for BO/ovcor corr.
*     hulj: local variable; temporary storage for BO/ovcor corr.
*     corr1: local variable; temporary storage for BO/ovcor corr.
*     corr2: local variable; temporary storage for BO/ovcor corr.
*     dbodsboi2: local variable; derivative of 1-3 BO correction to sum 
*                  of bond orders around atom i 
*     dbodsboj2: local variable; derivative of 1-3 BO correction to sum 
*                  of bond orders around atom i 
*     bocor1: local variable; 1-3 bond order correction                            
*     bocor2: local variable; 1-3 bond order correction                            
*     ovi2: local variable; overcoordination on atom i with reference to
*                  total number of electrons on atom i, including lone
*                  pairs     
*     ovj2: local variable; overcoordination on atom j with reference to
*                  total number of electrons on atom j, including lone
*                  pairs     
*     valf(nsort): stored in cbka.blk; total number of electrons on 
*                  atom, including lone pairs
*     cor1: local variable; temporary storage for BO/1-3 bond corr.
*     cor2: local variable; temporary storage for BO/1-3 bond corr.
*     exphu3: local variable; storage exponential
*     exphu4: local variable; storage exponential
*     corrtot2: local variable; square of corrtot                                  
*     dbodboo: local variable; derivative of corrected total bond order to 
*                   uncorrected bond order
*     dbopidbopio: local variable; derivative of corrected pi bond order
*                   to uncorrected pi bond order
*     dbopidboo: local variable; derivative of corrected pi bond order
*                   to uncorrected total bond order
*     dbopi2dbopi2o: local variable; derivative of corrected double pi bond order
*                   to uncorrected double pi bond order
*     dbopi2dboo: local variable; derivative of corrected double pi bond order
*                   to uncorrected total bond order
*     dbodsboit: local variable; derivative of total bond order to sum
*                   of bond orders around atom i
*     dbodsbojt: local variable; derivative of total bond order to sum
*                   of bond orders around atom j
*     vhui: local variable; temporary storage              
*     vhuj: local variable; temporary storage              
*     dbopidsboit: local variable; derivative of pi bond order to sum
*                   of bond orders around atom i
*     dbopidsbojt: local variable; derivative of pi bond order to sum
*                   of bond orders around atom j
*     dbopi2dsboit: local variable; derivative of pi bond order to sum
*                   of bond orders around atom i
*     dbopi2dsbojt: local variable; derivative of pi bond order to sum
*                   of bond orders around atom j
*     nco: local integer; counter for number of atoms in derivative
*     ihl: local integer; helps to access right dbodc-term
*     nubon2(nat,mbond): stored in cbka.blk; stored bond number as a
*                   function of atom number and connection number
*     iob: local integer; atom number of second atom in bond
*     ncubo: local integer; stores number of current bond  
*     na: stored in cbka.blk: number of atoms in system 
*     zero: stored in cbka.blk: value 0.00
*                                                          
********************************************************************** 
      do 10 i1=1,nbon
      ibt=ib(i1,1)
      j1=ib(i1,2)
      j2=ib(i1,3)
      if (ovc(ibt).lt.0.001.and.v13cor(ibt).lt.0.001) then
      idbo1(i1)=2
      idbo(i1,1)=j1
      idbo(i1,2)=j2
      do k1=1,3
      dbondc(i1,k1,1)=dbodc(i1,k1,1)
      dbondc(i1,k1,2)=dbodc(i1,k1,2)
      dbosindc(i1,k1,1)=dbosidc(i1,k1,1)
      dbosindc(i1,k1,2)=dbosidc(i1,k1,2)
      dbopindc(i1,k1,1)=dbopidc(i1,k1,1)
      dbopindc(i1,k1,2)=dbopidc(i1,k1,2)
      dbopi2ndc(i1,k1,1)=dbopi2dc(i1,k1,1)
      dbopi2ndc(i1,k1,2)=dbopi2dc(i1,k1,2)
      end do
      goto 10
      end if
      boo=bo(i1)
      bopio=bopi(i1)
      bopi2o=bopi2(i1)
      iti=ia(j1,1)
      itj=ia(j2,1)
      aboi=abo(j1)
      aboj=abo(j2)
      vp131=sqrt(bo131(iti)*bo131(itj))
      vp132=sqrt(bo132(iti)*bo132(itj))
      vp133=sqrt(bo133(iti)*bo133(itj))
      corrtot=1.0
      dbodsboi1=zero
      dbodsboj1=zero
      if (ovc(ibt).gt.0.001) then
      ovi=aboi-aval(iti)
      ovj=aboj-aval(itj)

********************************************************************** 
*                                                                    *
*     Correction for overcoordination                                *
*                                                                    *
********************************************************************** 
      exphu1=exp(-vpar(2)*ovi)
      exphu2=exp(-vpar(2)*ovj)
      exp11=exp(-vpar(1)*ovi)
      exp21=exp(-vpar(1)*ovj)
      exphu12=(exphu1+exphu2)
      ovcor=-(1.0/vpar(2))*log(0.50*exphu12)
*     huli=((1.0/ovc(ibt))*aval(iti)+exp11+exp21)
*     hulj=((1.0/ovc(ibt))*aval(itj)+exp11+exp21)
      huli=aval(iti)+exp11+exp21
      hulj=aval(itj)+exp11+exp21
      corr1=huli/(huli+ovcor)
      corr2=hulj/(hulj+ovcor)
      corrtot=0.50*(corr1+corr2)

      dbodsboi1=0.50*(-vpar(1)*exp11/(huli+ovcor)-
     $(corr1/(huli+ovcor))*
     $(-vpar(1)*exp11+exphu1/exphu12)-vpar(1)*exp11/(hulj+ovcor)-
     $(corr2/(hulj+ovcor))*(-vpar(1)*exp11+exphu1/exphu12))
      dbodsboj1=0.50*(-vpar(1)*exp21/(huli+ovcor)-
     $(corr1/(huli+ovcor))*
     $(-vpar(1)*exp21+exphu2/exphu12)-vpar(1)*exp21/(hulj+ovcor)-
     $(corr2/(hulj+ovcor))*(-vpar(1)*exp21+exphu2/exphu12))
      end if
********************************************************************** 
*                                                                    *
*     Correction for 1-3 bond orders                                 *
*                                                                    *
********************************************************************** 
      dbodsboi2=zero
      dbodsboj2=zero
      bocor1=1.0
      bocor2=1.0
      if (v13cor(ibt).gt.0.001) then
      ovi2=aboi-vval3(iti)                  !Modification for metal surfaces
      ovj2=aboj-vval3(itj)
*     ovi2=aboi-valf(iti)
*     ovj2=aboj-valf(itj)
*     ovi2=aboi-aval(iti)
*     ovj2=aboj-aval(itj)
      cor1=vp131*boo*boo-ovi2
      cor2=vp131*boo*boo-ovj2
*     exphu3=v13cor(ibt)*exp(-vp132*cor1+vp133)
*     exphu4=v13cor(ibt)*exp(-vp132*cor2+vp133)
      exphu3=exp(-vp132*cor1+vp133)
      exphu4=exp(-vp132*cor2+vp133)
      bocor1=1.0/(1.0+exphu3)
      bocor2=1.0/(1.0+exphu4)
      dbodsboi2=-bocor1*bocor1*bocor2*vp132*exphu3
      dbodsboj2=-bocor1*bocor2*bocor2*vp132*exphu4
      end if

      bo(i1)=boo*corrtot*bocor1*bocor2
      if (bo(i1).lt.1e-10) bo(i1)=zero
      corrtot2=corrtot*corrtot
      bopi(i1)=bopio*corrtot2*bocor1*bocor2
      bopi2(i1)=bopi2o*corrtot2*bocor1*bocor2
      if (bopi(i1).lt.1e-10) bopi(i1)=zero
      if (bopi2(i1).lt.1e-10) bopi2(i1)=zero

      dbodboo=corrtot*bocor1*bocor2+corrtot*
     $bocor1*bocor1*bocor2*boo*vp132*vp131*2.0*boo*exphu3+
     $corrtot*bocor1*bocor2*bocor2*boo*
     $vp132*vp131*exphu4*2.0*boo

      dbopidbopio=corrtot2*bocor1*bocor2

      dbopidboo=corrtot2*
     $bocor1*bocor1*bocor2*boo*vp132*vp131*2.0*bopio*exphu3+
     $corrtot2*bocor1*bocor2*bocor2*boo*
     $vp132*vp131*exphu4*2.0*bopio

      dbopi2dbopi2o=corrtot2*bocor1*bocor2

      dbopi2dboo=corrtot2*
     $bocor1*bocor1*bocor2*boo*vp132*vp131*2.0*bopi2o*exphu3+
     $corrtot2*bocor1*bocor2*bocor2*boo*
     $vp132*vp131*exphu4*2.0*bopi2o

      dbodsboit=boo*dbodsboi1*bocor1*bocor2+boo*corrtot*dbodsboi2
      dbodsbojt=boo*dbodsboj1*bocor1*bocor2+boo*corrtot*dbodsboj2

      vhui=2.0*corrtot*dbodsboi1*bocor1*bocor2+corrtot2*dbodsboi2
      vhuj=2.0*corrtot*dbodsboj1*bocor1*bocor2+corrtot2*dbodsboj2
      dbopidsboit=bopio*vhui
      dbopidsbojt=bopio*vhuj

      dbopi2dsboit=bopi2o*vhui
      dbopi2dsbojt=bopi2o*vhuj

********************************************************************** 
*                                                                    *
*     Calculate bond order derivatives                               *
*                                                                    *
********************************************************************** 
      idbo1(i1)=2+ia(j1,2)+ia(j2,2)
      idbo(i1,1)=j1
      idbo(i1,2)=j2
      nco=0
      do k1=1,3
      dbondc(i1,k1,1)=dbodc(i1,k1,1)*dbodboo
      dbondc(i1,k1,2)=dbodc(i1,k1,2)*dbodboo
*     dbosindc(i1,k1,1)=dbosidc(i1,k1,1)*dbosidboo
*     dbosindc(i1,k1,2)=dbosidc(i1,k1,2)*dbosidboo
      dbopindc(i1,k1,1)=dbopidc(i1,k1,1)*dbopidbopio+
     $dbodc(i1,k1,1)*dbopidboo
      dbopindc(i1,k1,2)=dbopidc(i1,k1,2)*dbopidbopio+
     $dbodc(i1,k1,2)*dbopidboo
      dbopi2ndc(i1,k1,1)=dbopi2dc(i1,k1,1)*dbopi2dbopi2o+
     $dbodc(i1,k1,1)*dbopi2dboo
      dbopi2ndc(i1,k1,2)=dbopi2dc(i1,k1,2)*dbopi2dbopi2o+
     $dbodc(i1,k1,2)*dbopi2dboo
      end do
      do i2=1,ia(j1,2)
      ihl=0
      iob=ia(j1,2+i2)
      if (iob.lt.j1) ihl=1
      ncubo=nubon2(j1,i2)
      idbo(i1,2+nco+1)=iob
      do k1=1,3
      dbondc(i1,k1,1)=dbondc(i1,k1,1)+dbodc(ncubo,k1,1+ihl)*dbodsboit
      dbondc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbodsboit

*     dbosindc(i1,k1,1)=dbosindc(i1,k1,1)+
*    $dbodc(ncubo,k1,1+ihl)*dbosidsboit
*     dbosindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbosidsboit

      dbopindc(i1,k1,1)=dbopindc(i1,k1,1)+
     $dbodc(ncubo,k1,1+ihl)*dbopidsboit
      dbopindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopidsboit

      dbopi2ndc(i1,k1,1)=dbopi2ndc(i1,k1,1)+
     $dbodc(ncubo,k1,1+ihl)*dbopi2dsboit
      dbopi2ndc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopi2dsboit

      end do
      nco=nco+1
      end do
      do i2=1,ia(j2,2)
      ihl=0
      iob=ia(j2,2+i2)
      if (iob.lt.j2) ihl=1
      ncubo=nubon2(j2,i2)
      idbo(i1,2+nco+1)=iob
      do k1=1,3

      dbondc(i1,k1,2)=dbondc(i1,k1,2)+dbodc(ncubo,k1,1+ihl)*dbodsbojt
      dbondc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbodsbojt

*     dbosindc(i1,k1,2)=dbosindc(i1,k1,2)+
*    $dbodc(ncubo,k1,1+ihl)*dbosidsbojt
*     dbosindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbosidsbojt

      dbopindc(i1,k1,2)=dbopindc(i1,k1,2)+
     $dbodc(ncubo,k1,1+ihl)*dbopidsbojt
      dbopindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopidsbojt

      dbopi2ndc(i1,k1,2)=dbopi2ndc(i1,k1,2)+
     $dbodc(ncubo,k1,1+ihl)*dbopi2dsbojt
      dbopi2ndc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopi2dsbojt

      end do
      nco=nco+1
      end do

   10 continue

      do i1=1,na
      abo(i1)=zero
      end do
*     do i1=1,na
*     do i2=1,ia(i1,2)
*     iob=ia(i1,2+i2)
*     ncubo=nubon2(i1,i2)
*     abo(i1)=abo(i1)+bo(ncubo)
*     end do
*     end do
      do i1=1,nbon
      j1=ib(i1,2)
      j2=ib(i1,3)
      abo(j1)=abo(j1)+bo(i1)
      if (j1.ne.j2) abo(j2)=abo(j2)+bo(i1)
      end do

   15 continue
      return
      end
********************************************************************** 
********************************************************************** 

      subroutine lonpar

********************************************************************** 
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkconst.blk"
#include "cbkc.blk"
#include "cbkd.blk"
#include "cbkdcell.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbklonpar.blk"
#include "cbknubon2.blk"
#include "control.blk"
#include "small.blk"
      dimension virial_tmp(3,3),virialsym(6)
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In lonpar'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
**********************************************************************
*                                                                    *
*     Calculate lone pair energy and first derivatives               *
*                                                                    *
**********************************************************************
      elp=zero
      do i1=1,na
**********************************************************************
*                                                                    *
*     Determine number of lone pairs on atoms
*                                                                    *
**********************************************************************
      ity=ia(i1,1)
      voptlp=0.50*(stlp(ity)-aval(ity))
      vlp(i1)=zero
      vund=abo(i1)-stlp(ity)
      vlph=2.0*int(vund/2.0)
      vlpex=vund-vlph
      vp16h=vpar(16)-1.0

      expvlp=exp(-vpar(16)*(2.0+vlpex)*(2.0+vlpex))
      dvlpdsbo(i1)=-vpar(16)*2.0*(2.0+vlpex)*expvlp
      vlp(i1)=expvlp-int(vund/2.0)
*     expvlp=exp(-vpar(16)*(2.0+vlpex))
*     dvlpdsbo(i1)=-vpar(16)*expvlp
*     expvlp=exp(-6.0*((-0.50*vlpex)**vpar(16)))
*     vlp(i1)=(1.0-expvlp)-int(vund/2.0)
*     dvlpdsbo(i1)=-0.5*6.0*vpar(16)*((-0.5*vlpex)**vp16h)*
*    $expvlp
**********************************************************************
*                                                                    *
*     Calculate lone pair energy                                     *
*                                                                    *
**********************************************************************
      if (i1 .le. na_local) then

      diffvlp=voptlp-vlp(i1)
      exphu1=exp(-75.0*diffvlp)
      hulp1=1.0/(1.0+exphu1)
      elph=vlp1(ity)*diffvlp*hulp1
*     elph=vlp1(ity)*diffvlp
      delpdvlp=-vlp1(ity)*hulp1-vlp1(ity)*diffvlp*hulp1*hulp1*
     $75.0*exphu1

      elp=elp+elph
      estrain(i1)=estrain(i1)+elph !atom energy

      delpdsbo=delpdvlp*dvlpdsbo(i1)
**********************************************************************
*                                                                    *
*     Calculate first derivative of lone pair energy to              *
*     cartesian coordinates                                          *
*                                                                    *
**********************************************************************
      do i3=1,ia(i1,2)
      iob=ia(i1,2+i3)
      ncubo=nubon2(i1,i3)

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = delpdsbo*dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
         do k1p=1,3
            virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
         end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif
      endif

      end do
      endif

      end do

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine covbon

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkabo.blk"
#include "cbkbo.blk"
#include "cbkbosi.blk"
#include "cbkbopi.blk"
#include "cbkbopi2.blk"
#include "cbkconst.blk"
#include "cbkcovbon.blk"
#include "cbkd.blk"
#include "cbkdbopi2ndc.blk"
#include "cbkdbopindc.blk"
#include "cbkdcell.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbknubon2.blk"
#include "cbkqa.blk"
#include "cbkrbo.blk"
#include "control.blk"
#include "small.blk"
      dimension virial_tmp(3,3),virialsym(6)
********************************************************************** 
*                                                                    *
*     Calculate bond energy and first derivatives                    *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In covbon'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      eb=0.0d0
      if (nbon.eq.0) return
********************************************************************** 
*                                                                    *
*     Calculate bond energies                                        *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write(65,*) 'Bond forces'
c$$$      write(65,*) 'nbon = ',nbon
c$$$      endif

      do 20 i1=1,nbon

      boa=bo(i1)
*     if (boa.lt.cutof2) goto 20
      j1=ib(i1,2)
      j2=ib(i1,3)

c     Only compute interaction if both atoms
c     are local or else flip a coin 
      if (j1 .gt. na_local) go to 20
      if (j2 .gt. na_local) then
         if (itag(j1) .lt. itag(j2)) go to 20
         if (itag(j1) .eq. itag(j2)) then
            if(c(j1,3) .gt. c(j2,3)) go to 20
            if(c(j1,3) .eq. c(j2,3) .and. 
     $           c(j1,2) .gt. c(j2,2)) go to 20
            if(c(j1,3) .eq. c(j2,3) .and. 
     $           c(j1,2) .eq. c(j2,2) .and. 
     $           c(j1,1) .gt. c(j2,1)) go to 20
         endif
      endif
      vsymm=1.0
      if (j1.eq.j2) vsymm=0.5

      bopia=bopi(i1)
      bopi2a=bopi2(i1)
      bosia=boa-bopia-bopi2a
      if (bosia.lt.zero) bosia=zero
      it1=ia(j1,1)
      it2=ia(j2,1)
      ibt=ib(i1,1)
      de1h=vsymm*de1(ibt)
      de2h=vsymm*de2(ibt)
      de3h=vsymm*de3(ibt)

      bopo1=bosia**psp(ibt)
      exphu1=exp(psi(ibt)*(1.0-bopo1))
      ebh=-de1h*bosia*exphu1-de2h*bopia-de3h*bopi2a

      debdbo=-de1h*exphu1+de1h*exphu1*psp(ibt)*psi(ibt)*bopo1
      debdbopi=-de2h
      debdbopi2=-de3h

      eb=eb+ebh
      estrain(j1)=estrain(j1)+0.50*ebh !1st atom energy
      estrain(j2)=estrain(j2)+0.50*ebh !2nd atom energy

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(i1)
      ihu=idbo(i1,i2)
      do k1=1,3
      ftmp = debdbo*(dbondc(i1,k1,i2)-dbopindc(i1,k1,i2)-
     $dbopi2ndc(i1,k1,i2))+
     $debdbopi*dbopindc(i1,k1,i2)+
     $debdbopi2*dbopi2ndc(i1,k1,i2)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(i1)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(i1)
                  ihu=idbo(i1,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

**********************************************************************
*                                                                    *
*     Stabilisation terminal triple bond in CO                       *
*                                                                    *
**********************************************************************
      if (boa.lt.1.00) goto 20
* Stabilization for all triple bonds (not just for CO) in ReaxFF combustion FF
      if (ltripstaball.eq.1 .or.
     $  (qa(j1).eq.'C '.and.qa(j2).eq.'O ').or.
     $  (qa(j1).eq.'O '.and.qa(j2).eq.'C ')) then

      ba=(boa-2.50)*(boa-2.50)
      exphu=exp(-vpar(8)*ba)
      oboa=abo(j1)-boa
      obob=abo(j2)-boa
      exphua1=exp(-vpar(4)*oboa)
      exphub1=exp(-vpar(4)*obob)
      ovoab=abo(j1)-aval(it1)+abo(j2)-aval(it2)
      exphuov=exp(vpar(5)*ovoab)
      hulpov=1.0/(1.0+25.0*exphuov)
      
      estriph=vpar(11)*exphu*hulpov*(exphua1+exphub1)

      eb=eb+estriph
      estrain(j1)=estrain(j1)+0.50*estriph !1st atom energy
      estrain(j2)=estrain(j2)+0.50*estriph !2nd atom energy

      decobdbo=vpar(4)*vpar(11)*exphu*hulpov*(exphua1+exphub1)
     $-2.0*vpar(11)*vpar(8)*(boa-2.50)*hulpov*exphu*
     $(exphua1+exphub1)
      decobdboua=-25.0*vpar(5)*vpar(11)*exphu*exphuov*hulpov*hulpov*
     $(exphua1+exphub1)-vpar(11)*exphu*vpar(4)*hulpov*exphua1
      decobdboub=-25.0*vpar(5)*vpar(11)*exphu*exphuov*hulpov*hulpov*
     $(exphua1+exphub1)-vpar(11)*exphu*vpar(4)*hulpov*exphub1

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(i1)
      ihu=idbo(i1,i2)
      do k1=1,3
      ftmp = decobdbo*dbondc(i1,k1,i2)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(i1)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(i1)
                  ihu=idbo(i1,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      do i3=1,ia(j1,2)
      iob=ia(j1,2+i3)
      ncubo=nubon2(j1,i3)

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = decobdboua*dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      end do

      do i3=1,ia(j2,2)
      iob=ia(j2,2+i3)
      ncubo=nubon2(j2,i3)
      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif
      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = decobdboub*dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      end do

      endif

   20 continue

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine ovcor

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkabo.blk"
#include "cbkbo.blk"
#include "cbkbopi.blk"
#include "cbkbopi2.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkdbopi2ndc.blk"
#include "cbkdbopindc.blk"
#include "cbkdcell.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbklonpar.blk"
#include "cbknubon2.blk"
#include "cbkrbo.blk"
#include "control.blk"
#include "small.blk"
********************************************************************** 
*                                                                    *
*     Calculate atom energy                                          *
*     Correction for over- and undercoordinated atoms                *
*                                                                    *
********************************************************************** 
      dimension vlptemp(nat)
      dimension virial_tmp(3,3),virialsym(6)
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In ovcor'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      do i1=1,na
      ity1=ia(i1,1)
      vlptemp(i1)=vlp(i1)
      if (amas(ity1).gt.21.0) vlptemp(i1)=0.50*(stlp(ity1)-aval(ity1))  !Only for 1st-row elements
      end do
   25 ea=zero
      eaot=zero
      eaut=zero
      epen=0.0

      do 30 i1=1,na_local
      ity1=ia(i1,1)
      dfvl=1.0                         
      if (amas(ity1).gt.21.0) dfvl=0.0  !Only for 1st-row elements
********************************************************************** 
*                                                                    *
*     Calculate overcoordination energy                              *
*     Valency is corrected for lone pairs                            *
*                                                                    *
********************************************************************** 

      voptlp=0.50*(stlp(ity1)-aval(ity1))
      diffvlph=dfvl*(voptlp-vlptemp(i1))
********************************************************************** 
*                                                                    *
*     Determine coordination neighboring atoms                       *
*                                                                    *
********************************************************************** 
      sumov=0.0
      sumov2=0.0
      do i3=1,ia(i1,2)
      iat2=ia(i1,2+i3)
      ity2=ia(iat2,1)
      ncubo=nubon2(i1,i3)
      if (bo(ncubo).gt.0.0) then
      ibt=ib(ncubo,1)
      voptlp2=0.50*(stlp(ity2)-aval(ity2))
      diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
      sumov=sumov+(bopi(ncubo)+bopi2(ncubo))*
     $(abo(iat2)-aval(ity2)-diffvlp2)
      sumov2=sumov2+vover(ibt)*de1(ibt)*bo(ncubo)
      endif
      end do

      exphu1=exp(vpar(32)*sumov)
      vho=1.0/(1.0+vpar(33)*exphu1)
      diffvlp=diffvlph*vho
      
      vov1=abo(i1)-aval(ity1)-diffvlp
      dvov1dsumov=diffvlph*vpar(32)*vpar(33)*vho*vho*exphu1
      exphuo=exp(vovun(ity1)*vov1)
      hulpo=1.0/(1.0+exphuo)

      hulpp=(1.0/(vov1+aval(ity1)+1e-8))

      eah=sumov2*hulpp*hulpo*vov1
      deadvov1=-sumov2*hulpp*hulpp*vov1*hulpo+
     $sumov2*hulpp*hulpo-sumov2*hulpp*vov1*vovun(ity1)*
     $hulpo*hulpo*exphuo

      ea=ea+eah
      estrain(i1)=estrain(i1)+eah !atom energy
********************************************************************** 
*                                                                    *
*     Calculate first derivative of overcoordination energy to       *
*     cartesian coordinates                                          *
*                                                                    *
********************************************************************** 
      do i3=1,ia(i1,2)
      iob=ia(i1,2+i3)
      ncubo=nubon2(i1,i3)
      if (bo(ncubo).gt.0.0) then
      ibt=ib(ncubo,1)
      deadbo=vover(ibt)*de1(ibt)*hulpp*hulpo*vov1

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = deadvov1*(1.0+dfvl*vho*dvlpdsbo(i1))*
     $dbondc(ncubo,k1,i4)+deadbo*dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do
 
      do i2=1,ia(i1,2)
 
      iat2=ia(i1,2+i2)
      ity2=ia(iat2,1)
      nbosa=nubon2(i1,i2)
      if (bo(nbosa).gt.0.0) then
      deadvov2=deadvov1*dvov1dsumov*(bopi(nbosa)+bopi2(nbosa))
 
      voptlp2=0.50*(stlp(ity2)-aval(ity2))
      diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
      deadpibo=deadvov1*dvov1dsumov*(abo(iat2)-aval(ity2)-diffvlp2)
 
      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(nbosa)
      ihu=idbo(nbosa,i4)
      do k1=1,3
      ftmp = deadpibo*(dbopindc(nbosa,k1,i4)+
     $dbopi2ndc(nbosa,k1,i4))
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(nbosa)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(nbosa)
                  ihu=idbo(nbosa,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif
      endif
 
      do i3=1,ia(iat2,2)
      iob=ia(iat2,2+i3)
      ncubo=nubon2(iat2,i3)
      if (bo(ncubo).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = deadvov2*(1.0+dfvl*dvlpdsbo(iat2))*
     $dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp
 
      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do

      endif

      end do

********************************************************************** 
*                                                                    *
*     Calculate undercoordination energy                             *
*                                                                    *
********************************************************************** 
      if (valp1(ity1).lt.zero) goto 30  !skip undercoordination
      exphu2=exp(vpar(10)*sumov)
      vuhu1=1.0+vpar(9)*exphu2
      hulpu2=1.0/vuhu1

      exphu3=-exp(vpar(7)*vov1)
      hulpu3=-(1.0+exphu3)

      dise2=valp1(ity1)
      exphuu=exp(-vovun(ity1)*vov1)
      hulpu=1.0/(1.0+exphuu)
      eahu=dise2*hulpu*hulpu2*hulpu3
      deaudvov1=dise2*hulpu2*vovun(ity1)*hulpu*hulpu*exphuu*hulpu3-
     $dise2*hulpu*hulpu2*vpar(7)*exphu3

      ea=ea+eahu
      estrain(i1)=estrain(i1)+eahu !atom energy

      deaudsumov=-dise2*hulpu*vpar(9)*vpar(10)*hulpu3*exphu2*
     $hulpu2*hulpu2

********************************************************************** 
*                                                                    *
*     Calculate first derivative of atom energy to cartesian         *
*     coordinates                                                    *
*                                                                    *
********************************************************************** 

      do i3=1,ia(i1,2)
      iob=ia(i1,2+i3)
      ncubo=nubon2(i1,i3)
      if (bo(ncubo).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = deaudvov1*(1.0+dfvl*vho*dvlpdsbo(i1))*
     $dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do
 
      do i2=1,ia(i1,2)
 
      iat2=ia(i1,2+i2)
      ity2=ia(iat2,1)
      nbosa=nubon2(i1,i2)
      if (bo(nbosa).gt.0.0) then
      deadvov2=(deaudsumov+dvov1dsumov*deaudvov1)*
     $(bopi(nbosa)+bopi2(nbosa))
 
      voptlp2=0.50*(stlp(ity2)-aval(ity2))
      diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
      deadpibo1=(dvov1dsumov*deaudvov1+deaudsumov)*
     $(abo(iat2)-aval(ity2)-diffvlp2)
 
      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(nbosa)
      ihu=idbo(nbosa,i4)
      do k1=1,3
      ftmp = deadpibo1*
     $(dbopindc(nbosa,k1,i4)+dbopi2ndc(nbosa,k1,i4))
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(nbosa)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(nbosa)
                  ihu=idbo(nbosa,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      do i3=1,ia(iat2,2)
      iob=ia(iat2,2+i3)
      ncubo=nubon2(iat2,i3)
      if (bo(ncubo).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = deadvov2*(1.0+dfvl*dvlpdsbo(iat2))*
     $dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do
 
      endif

      end do


   30 continue
 
**********************************************************************
*                                                                    *
*     Calculate correction for C2                                    *
*                                                                    *
**********************************************************************
      if (abs(vpar(6)).gt.0.001) then
      do 40 i1=1,na_local
      ity1=ia(i1,1)
      vov4=abo(i1)-aval(ity1)
 
      do i2=1,ia(i1,2)
      iat2=ia(i1,2+i2)
      nbohu=nubon2(i1,i2)
      if (bo(nbohu).gt.0.0) then
 
      ibt=ib(nbohu,1)
      elph=zero
      deahu2dbo=zero
      deahu2dsbo=zero
      vov3=bo(nbohu)-vov4-0.040*(vov4**4)
      if (vov3.gt.3.0) then
      elph=vpar(6)*(vov3-3.0)*(vov3-3.0)
      deahu2dbo=2.0*vpar(6)*(vov3-3.0)
      deahu2dsbo=2.0*vpar(6)*(vov3-3.0)*(-1.0-
     $0.16*(vov4**3))
      end if
 
      elp=elp+elph
      estrain(i1)=estrain(i1)+elph !atom energy

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(nbohu)
      ihu=idbo(nbohu,i3)
      do k1=1,3
      ftmp = deahu2dbo*dbondc(nbohu,k1,i3)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(nbohu)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(nbohu)
                  ihu=idbo(nbohu,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      do i3=1,ia(i1,2)
      iob=ia(i1,2+i3)
      ncubo=nubon2(i1,i3)
      if (bo(ncubo).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = deahu2dsbo*dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      end if
      end do

      end if
      end do
 
   40 continue
      end if

      return
      end
**********************************************************************
********************************************************************** 

      subroutine molen

********************************************************************** 
#include "cbka.blk"
#include "cbkbo.blk"
#include "cbkconst.blk"
#include "cbkc.blk"
#include "cbkd.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbknmolat.blk"
#include "cbknubon2.blk"
#include "control.blk"
#include "small.blk"
      dimension virial_tmp(3,3),virialsym(6)
********************************************************************** 
*                                                                    *
*     Calculate molecular energy and first derivatives               *
*     Only used to prevent creating virtual electrons                *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In molen'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      emol=zero
      return
      do i1=1,nmolo

      enelm=0.0
      do i2=1,na
      if (ia(i2,3+mbond).eq.i1) then
      it1=ia(i2,1)
      enelm=enelm+aval(it1)
      end if
      end do

      na1m=nmolat(i1,1)

      enelm=2*int(enelm*0.50)
*     enelm=elmol(i1)
      bomsum=zero
      do i2=1,na1m
      ihu=nmolat(i1,i2+1)
      do i3=1,ia(ihu,2)
      ihu2=nubon2(ihu,i3)
      bomsum=bomsum+bo(ihu2)
      end do
      end do
      diff=(bomsum-enelm)
      exphu=exp(-vpar(37)*diff)
      exphu2=1.0/(1.0+15.0*exphu)
      emolh=zero
      demoldsbo=zero
      emolh=vpar(38)*exphu2
      emol=emol+emolh
      demoldsbo=vpar(38)*vpar(37)*15.0*exphu2*exphu2*exphu

      do i2=1,na1m
      ihu1=nmolat(i1,i2+1)
      do i3=1,ia(ihu1,2)
      iob=ia(ihu1,2+i3)
      ncubo=nubon2(ihu1,i3)

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i4=1,idbo1(ncubo)
      ihu=idbo(ncubo,i4)
      do k1=1,3
      ftmp = demoldsbo*dbondc(ncubo,k1,i4)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i4=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i4)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      end do
      end do


      end do


      return
      end
********************************************************************** 
********************************************************************** 

      subroutine valang

********************************************************************** 
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkbo.blk"
#include "cbkbopi.blk"
#include "cbkbopi2.blk"
#include "cbkconst.blk"
#include "cbkc.blk"
#include "cbkd.blk"
#include "cbkdbopi2ndc.blk"
#include "cbkdbopindc.blk"
#include "cbkdcell.blk"
#include "cbkdhdc.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkh.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbklonpar.blk"
#include "cbknubon2.blk"
#include "cbkvalence.blk"
#include "control.blk"
#include "valang.blk"
#include "small.blk"
      dimension j(3)
      dimension virial_tmp(3,3),virialsym(6)
********************************************************************** 
*                                                                    *
*     Calculate valency angle energies and first derivatives         *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In valang'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
*     eco=0.0
      ev=0.0
      ecoa=0.0
      if (nval.eq.0) return

      do 10 i1=1,nval
      ity=iv(i1,1)
      j(1)=iv(i1,2)
      j(2)=iv(i1,3)
      j(3)=iv(i1,4)

      if (j(2) .le. na_local) then

      la=iv(i1,5)
      lb=iv(i1,6)
      boa=bo(la)-cutof2
      bob=bo(lb)-cutof2
      if (boa.lt.zero.or.bob.lt.zero) goto 10

      hl=h(i1)     ! Calculated earlier in routine calval
********************************************************************** 
*                                                                    *
*     Calculate valency angle energy                                 *
*                                                                    *
********************************************************************** 
      nbocen=ia(j(2),2)
      sbo2=0.0
      vmbo=1.0
     
      do i2=1,nbocen
      ibv=nubon2(j(2),i2)
      if (bo(ibv).gt.0.0) then
      vmbo=vmbo*exp(-bo(ibv)**8)
      sbo2=sbo2+bopi(ibv)+bopi2(ibv)
      endif
      end do

      ity2=ia(j(2),1)
*     exbo=abo(j(2))-stlp(ia(j(2),1))
      exbo=abo(j(2))-valf(ity2)
*     if (exbo.gt.zero) exbo=zero
*     expov=exp(vka8(ity)*exbo)
*     expov2=exp(-vpar(13)*exbo)
*     htov1=2.0+expov2
*     htov2=1.0+expov+expov2
*     evboadj=htov1/htov2
      evboadj=1.0
      expun=exp(-vkac(ity)*exbo)
      expun2=exp(vpar(15)*exbo)
      htun1=2.0+expun2
      htun2=1.0+expun+expun2
      evboadj2=vval4(ity2)-(vval4(ity2)-1.0)*htun1/htun2
********************************************************************** 
*                                                                    *
*     Calculate number of lone pairs                                 *
*                                                                    *
********************************************************************** 
      dsbo2dvlp=(1.0-vmbo)
      vlpadj=zero
      exlp1=abo(j(2))-stlp(ia(j(2),1))
      exlp2=2.0*int(exlp1/2.0)
      exlp=exlp1-exlp2
      if (exlp.lt.zero) then
*     expvlp=exp(-vpar(16)*(2.0+exlp)*(2.0+exlp))
*     vlpadj=expvlp-int(exlp1/2.0)
*     dsbo2dvlp=(1.0-vmbo)*(1.0-vpar(34)*2.0*
*     $(2.0+exlp)*vpar(16)*expvlp)
      vlpadj=vlp(j(2))
      dsbo2dvlp=(1.0-vmbo)*(1.0+vpar(34)*dvlpdsbo(j(2)))
      end if

      sbo2=sbo2+(1.0-vmbo)*(-exbo-vpar(34)*vlpadj)
      dsbo2dvmbo=exbo+vpar(34)*vlpadj

      sbo2h=sbo2
      powv=vpar(17)
      if (sbo2.le.0.0) sbo2h=0.0
      if (sbo2.gt.0.0.and.sbo2.le.1.0) sbo2h=sbo2**powv
      if (sbo2.gt.1.0.and.sbo2.lt.2.0) sbo2h=2.0-(2.0-sbo2)**powv
      if (sbo2.gt.2.0) sbo2h=2.0
      thba=th0(ity)
      expsbo=exp(-vpar(18)*(2.0-sbo2h))
      thetao=180.0-thba*(1.0-expsbo)

      thetao=thetao*dgrrdn
      thdif=(thetao-hl)
      thdi2=thdif*thdif
      dthsbo=dgrrdn*thba*vpar(18)*expsbo
      if (sbo2.lt.0.0) dthsbo=zero
      if (sbo2.gt.0.0.and.sbo2.le.1.0) 
     $dthsbo=powv*(sbo2**(powv-1.0))*dgrrdn*thba*vpar(18)*expsbo
      if (sbo2.gt.1.0.and.sbo2.lt.2.0) 
     $dthsbo=powv*((2.0-sbo2)**(powv-1.0))*dgrrdn*thba*vpar(18)*expsbo
      if (sbo2.gt.2.0) dthsbo=zero

      exphu=vka(ity)*exp(-vka3(ity)*thdi2)
      exphu2=vka(ity)-exphu
      if (vka(ity).lt.zero) exphu2=exphu2-vka(ity)             !To avoid linear Me-H-Me angles (6/6/06)
      boap=boa**vval2(ity)
      boap2=boa**(vval2(ity)-1.0)
      bobp=bob**vval2(ity)
      bobp2=bob**(vval2(ity)-1.0)
      exa=exp(-vval1(ity2)*boap)
      exb=exp(-vval1(ity2)*bobp)
      dexadboa=vval2(ity)*vval1(ity2)*exa*boap2
      dexbdbob=vval2(ity)*vval1(ity2)*exb*bobp2
      exa2=(1.0-exa)
      exb2=(1.0-exb)

      evh=evboadj2*evboadj*exa2*exb2*exphu2
      devdlb=evboadj2*evboadj*dexbdbob*exa2*exphu2
      devdla=evboadj2*evboadj*dexadboa*exb2*exphu2
      devdsbo=2.0*evboadj2*evboadj*dthsbo*exa2*exb2*
     $vka3(ity)*thdif*exphu
      devdh=-2.0*evboadj2*evboadj*exa2*exb2*vka3(ity)*thdif*exphu

      devdsbo2=
     $evboadj*exa2*exb2*exphu2*(vval4(ity2)-1.0)*(-vpar(15)*expun2/htun2
     $+htun1*(vpar(15)*expun2-vkac(ity)*expun)/(htun2*htun2))

*     devdsbo2=-evboadj2*exa2*exb2*exphu2*(vpar(13)*expov2/htov2+
*    $htov1*(vka8(ity)*expov-vpar(13)*expov2)/(htov2*htov2))+
*    $evboadj*exa2*exb2*exphu2*(vpar(14)-1.0)*(-vpar(15)*expun2/htun2
*    $+htun1*(vpar(15)*expun2-vkac(ity)*expun)/(htun2*htun2))

      if (j(2) .le. na_local) then
         ev=ev+evh
         estrain(j(2))=estrain(j(2))+evh         !central atom energy
      endif

*     write (64,'(4i8,18f8.2)')mdstep,j(1),j(2),j(3),sbo2,sbo2h,
*    $thetao*rdndgr,hl*rdndgr,bo(la),bo(lb),bopi(la),
*    $vlp(j(2)),exbo,vlpadj,vmbo,evh,ev,vka(ity)
********************************************************************** 
*                                                                    *
*     Calculate penalty for two double bonds in valency angle        *
*                                                                    *
********************************************************************** 
      exbo=abo(j(2))-aval(ia(j(2),1))
      expov=exp(vpar(22)*exbo)
      expov2=exp(-vpar(21)*exbo)
      htov1=2.0+expov2
      htov2=1.0+expov+expov2
      ecsboadj=htov1/htov2
      exphu1=exp(-vpar(20)*(boa-2.0)*(boa-2.0))
      exphu2=exp(-vpar(20)*(bob-2.0)*(bob-2.0))

      epenh=vkap(ity)*ecsboadj*exphu1*exphu2
      estrain(j(2))=estrain(j(2))+epenh
      epen=epen+epenh
      decoadboa=-2.0*vpar(20)*epenh*(boa-2.0)
      decoadbob=-2.0*vpar(20)*epenh*(bob-2.0)

      decdsbo2=-vkap(ity)*exphu1*exphu2*(vpar(21)*expov2/htov2+htov1*
     $(vpar(22)*expov-vpar(21)*expov2)/(htov2*htov2))
**********************************************************************
*                                                                    *
*     Calculate valency angle conjugation energy                     *
*                                                                    *
**********************************************************************
      unda=abo(j(1))-boa
*     ovb=abo(j(2))-valf(ia(j(2),1))
      ovb=abo(j(2))-vval3(ia(j(2),1))    !Modification for Ru  7/6/2004

      undc=abo(j(3))-bob
      ba=(boa-1.50)*(boa-1.50)
      bb=(bob-1.50)*(bob-1.50)
      exphua=exp(-vpar(31)*ba)
      exphub=exp(-vpar(31)*bb)
      exphuua=exp(-vpar(39)*unda*unda)
      exphuob=exp(vpar(3)*ovb)
      exphuuc=exp(-vpar(39)*undc*undc)
      hulpob=1.0/(1.0+exphuob)
      ecoah=vka8(ity)*exphua*exphub*exphuua*exphuuc*hulpob
      decodbola=-2.0*vka8(ity)*(boa-1.50)*vpar(31)*exphua*exphub
     $*exphuua*exphuuc*hulpob+vpar(39)*vka8(ity)*exphua*exphub*
     $exphuua*exphuuc*hulpob*2.0*unda
      decodbolb=-2.0*vka8(ity)*(bob-1.50)*vpar(31)*exphua*exphub
     $*exphuua*exphuuc*hulpob+vpar(39)*vka8(ity)*exphua*exphub*
     $exphuua*exphuuc*hulpob*2.0*undc
      decodboua=-2.0*unda*vka8(ity)*vpar(39)*exphua*exphub
     $*exphuua*exphuuc*hulpob
      decodbouc=-2.0*undc*vka8(ity)*vpar(39)*exphua*exphub
     $*exphuua*exphuuc*hulpob
      decodboob=-vka8(ity)*exphua*exphub*exphuua*exphuuc*hulpob*
     $hulpob*vpar(3)*exphuob
*     decodboob=zero
*     decodboua=zero
*     decodbouc=zero

      ecoa=ecoa+ecoah
      estrain(j(2))=estrain(j(2))+ecoah !central atom energy

********************************************************************** 
*                                                                    *
*     Calculate derivative valency energy to cartesian coordinates   *
*                                                                    *
********************************************************************** 
      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do k1=1,3
      do k2=1,3
      ftmp = devdh*dhdc(i1,k1,k2)
      d(k1,j(k2))=d(k1,j(k2))+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(j(k2),k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/3
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do k2=1,3
                  ihu=j(k2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(la)
      ihu=idbo(la,i2)
      do k1=1,3
      ftmp = (devdla+decoadboa+decodbola)*
     $dbondc(la,k1,i2)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(la)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(la)
                  ihu=idbo(la,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(lb)
      ihu=idbo(lb,i2)
      do k1=1,3
      ftmp = (devdlb+decoadbob+decodbolb)*
     $dbondc(lb,k1,i2)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(lb)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(lb)
                  ihu=idbo(lb,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      do i2=1,nbocen
      j5=ia(j(2),2+i2)
      ibv=nubon2(j(2),i2)
      if (bo(ibv).gt.0.0) then
      dvmbodbo=-vmbo*8.0*bo(ibv)**7

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(ibv)
      ihu=idbo(ibv,i3)
      do k1=1,3
      ftmp = (-dsbo2dvlp*devdsbo+devdsbo2+decdsbo2
     $+dvmbodbo*dsbo2dvmbo*devdsbo)*
     $dbondc(ibv,k1,i3)+devdsbo*(dbopindc(ibv,k1,i3)+
     $dbopi2ndc(ibv,k1,i3))
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ibv)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(ibv)
                  ihu=idbo(ibv,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do
 
      do i2=1,ia(j(1),2)
      j5=ia(j(1),2+i2)
      ibv=nubon2(j(1),i2)
      if (bo(ibv).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(ibv)
      ihu=idbo(ibv,i3)
      do k1=1,3
      ftmp = decodboua*dbondc(ibv,k1,i3)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ibv)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(ibv)
                  ihu=idbo(ibv,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do
 
      do i2=1,ia(j(2),2)
      j5=ia(j(2),2+i2)
      ibv=nubon2(j(2),i2)
      if (bo(ibv).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(ibv)
      ihu=idbo(ibv,i3)
      do k1=1,3
      ftmp = decodboob*dbondc(ibv,k1,i3)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ibv)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(ibv)
                  ihu=idbo(ibv,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do

      do i2=1,ia(j(3),2)
      j5=ia(j(3),2+i2)
      ibv=nubon2(j(3),i2)
      if (bo(ibv).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(ibv)
      ihu=idbo(ibv,i3)
      do k1=1,3
      ftmp = decodbouc*dbondc(ibv,k1,i3)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ibv)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(ibv)
                  ihu=idbo(ibv,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do

      endif

   10 continue

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine hbond

********************************************************************** 
#include "cbka.blk"
#include "cbkbo.blk"
#include "cbkconst.blk"
#include "cbkc.blk"
#include "cbkd.blk"
#include "cbkdcell.blk"
#include "cbkenergies.blk"
#include "cbkidbo.blk"
#include "cbksrthb.blk"
#include "control.blk"
#include "cbkhbond.blk"
#include "small.blk"
      dimension drda(3),j(3),dvdc(3,3),dargdc(3,3)
      dimension virial_tmp(3,3),virialsym(6)
********************************************************************** 
*                                                                    *
*     Calculate hydrogen bond energies and first derivatives         *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In hbond'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      ehb=zero
      do 10 i1=1,nhb
      ityhb=ihb(i1,1)
      j(1)=ihb(i1,2)
      j(2)=ihb(i1,3)
      j(3)=ihb(i1,4)
      la=ihb(i1,5)
      boa=bo(la)
      call dista2(j(2),j(3),rda,dxm,dym,dzm)
      drda(1)=dxm/rda
      drda(2)=dym/rda
      drda(3)=dzm/rda
      call calvalhb(j(1),j(2),j(3),ix,iy,iz,arg,hhb(i1),dvdc,dargdc)
      rhu1=rhb(ityhb)/rda
      rhu2=rda/rhb(ityhb)
      sinhu=sin(hhb(i1)/2.0)
      sin2=sinhu*sinhu
      exphu1=exp(-vhb1(ityhb)*boa)
      exphu2=exp(-vhb2(ityhb)*(rhu1+rhu2-2.0))
      if (lhbnew .eq. 0) then
         ehbh=(1.0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*sin2*sin2
      else
         ehbh=(1.0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2
      endif
      ehb=ehb+ehbh
      estrain(j(2))=estrain(j(2))+ehbh !2nd atom energy

********************************************************************** 
*                                                                    *
*     Calculate first derivatives                                    *
*                                                                    *
********************************************************************** 
      if (lhbnew .eq. 0) then
         dehbdbo=vhb1(ityhb)*exphu1*dehb(ityhb)*exphu2*sin2*sin2*
     $        sin2*sin2
         dehbdv=(1.0-exphu1)*dehb(ityhb)*exphu2*
     $        4.0*sin2*sin2*sin2*sinhu*cos(hhb(i1)/2.0)
         dehbdrda=(1.0-exphu1)*dehb(ityhb)*sin2*sin2*sin2*sin2*
     $        vhb2(ityhb)*(rhb(ityhb)/(rda*rda)-1.0/rhb(ityhb))*exphu2
      else
         dehbdbo=vhb1(ityhb)*exphu1*dehb(ityhb)*exphu2*sin2*sin2
         dehbdv=(1.0-exphu1)*dehb(ityhb)*exphu2*
     $        2.0*sin2*sinhu*cos(hhb(i1)/2.0)
         dehbdrda=(1.0-exphu1)*dehb(ityhb)*sin2*sin2*
     $        vhb2(ityhb)*(rhb(ityhb)/(rda*rda)-1.0/rhb(ityhb))*exphu2
      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do k1=1,3
      ftmp = dehbdrda*drda(k1)
      d(k1,j(2))=d(k1,j(2))+ftmp
      d(k1,j(3))=d(k1,j(3))-ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+
     $        ftmp*c(j(2),k1p)-ftmp*c(j(3),k1p)
      end do
      endif

      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/2
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               ihu = j(2)
               atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               ihu = j(3)
               atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do k1=1,3
      do k2=1,3
      ftmp = dehbdv*dvdc(k1,k2)
      d(k1,j(k2))=d(k1,j(k2))+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(j(k2),k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/3
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do k2=1,3
                  ihu=j(k2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(la)
      ihu=idbo(la,i2)
      do k1=1,3
      ftmp = dehbdbo*dbondc(la,k1,i2)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(la)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(la)
                  ihu=idbo(la,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

   10 continue
      return
      end 

********************************************************************** 
********************************************************************** 

      subroutine torang

********************************************************************** 
#include "cbka.blk"
#include "cbkabo.blk"
#include "cbkbo.blk"
#include "cbkbopi.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkdbopindc.blk"
#include "cbkdcell.blk"
#include "cbkdhdc.blk"
#include "cbkdrdc.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkfftorang.blk"
#include "cbkh.blk"
#include "cbkia.blk"
#include "cbkidbo.blk"
#include "cbkinit.blk"
#include "cbknubon2.blk"
#include "cbkrbo.blk"
#include "cbktorang.blk"
#include "cbktorsion.blk"
#include "cbktregime.blk"
#include "cbkvalence.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"

      DIMENSION  A(3),DRDA(3),DADC(4),DRADC(3,4),DRBDC(3,4),
     $DRCDC(3,4),DHDDC(3,4),DHEDC(3,4),DRVDC(3,4),DTDC(3,4),
     $DNDC(3,4)
      dimension j(4),dh1rdc(3,3),dh2rdc(3,3),dargdc(3,3)
      dimension virial_tmp(3,3),virialsym(6)
********************************************************************** 
*                                                                    *
*     Calculate torsion angle energies and first derivatives         *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In torang'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      do k1=1,3
      do k2=1,4
      dhddc(k1,k2)=0.0
      dhedc(k1,k2)=0.0
      dradc(k1,k2)=0.0
      drbdc(k1,k2)=0.0
      drcdc(k1,k2)=0.0
      end do
      end do
      et=0.0
      eth12=0.0
      eco=0.0
      dadc(1)=1.0
      dadc(2)=0.0
      dadc(3)=0.0
      dadc(4)=-1.0
      if (ntor.eq.0) return

      do 10 i1=1,ntor
      j(1)=it(i1,2)
      j(2)=it(i1,3)
      j(3)=it(i1,4)
      j(4)=it(i1,5)

      ity=it(i1,1)  
      la=it(i1,6)
      lb=it(i1,7)
      lc=it(i1,8)
      call calvalres(j(1),j(2),j(3),arg1,ht1,dh1rdc,dargdc)
      call calvalres(j(2),j(3),j(4),arg2,ht2,dh2rdc,dargdc)
      boa=bo(la)-cutof2
      bob=bo(lb)-cutof2
      boc=bo(lc)-cutof2
      if (boa.lt.zero.or.bob.lt.zero.or.boc.lt.zero) 
     $goto 10
      r42=0.0
      ivl1=ibsym(la)
      ivl2=ibsym(lb)
      ivl3=ibsym(lc)
      isign1=1
      isign2=1
      isign3=1
      rla=rbo(la)
      rlb=rbo(lb)
 
      call dista2(j(1),j(4),r4,a(1),a(2),a(3))
********************************************************************** 
*                                                                    *
*     Determine torsion angle                                        *
*                                                                    *
********************************************************************** 
      d142=r4*r4
      rla=rbo(la)
      rlb=rbo(lb)
      rlc=rbo(lc)
      coshd=cos(ht1)
      coshe=cos(ht2)
      sinhd=sin(ht1)
      sinhe=sin(ht2)
      poem=2.0*rla*rlc*sinhd*sinhe
      poem2=poem*poem
      tel=rla*rla+rlb*rlb+rlc*rlc-d142-2.0*(rla*rlb*coshd-rla*rlc*
     $coshd*coshe+rlb*rlc*coshe)
      if (poem.lt.1e-20) poem=1e-20
      arg=tel/poem
      if (arg.gt.1.0) arg=1.0
      if (arg.lt.-1.0) arg=-1.0
      arg2=arg*arg
      thg(i1)=acos(arg)*rdndgr
      k1=j(1)
      k2=j(2)
      k3=j(3)
      k4=j(4)
      call dista2(k3,k2,dis,x3,y3,z3)
      y32z32=y3*y3+z3*z3
      wort1=sqrt(y32z32)+1e-6
      wort2=sqrt(y32z32+x3*x3)+1e-6
*     if (wort1.lt.1e-6) wort1=1e-6
*     if (wort2.lt.1e-6) wort2=1e-6
      sinalf=y3/wort1
      cosalf=z3/wort1
      sinbet=x3/wort2
      cosbet=wort1/wort2
      call dista2(k1,k2,dis,x1,y1,z1)
      x1=x1*cosbet-y1*sinalf*sinbet-z1*cosalf*sinbet
      y1=y1*cosalf-z1*sinalf
      wort3=sqrt(x1*x1+y1*y1)+1e-6
*     if (wort3.lt.1e-6) wort3=1e-6
      singam=y1/wort3
      cosgam=x1/wort3
      call dista2(k4,k2,dis,x4,y4,z4)
      x4=x4*cosbet-y4*sinalf*sinbet-z4*cosalf*sinbet
      y4=y4*cosalf-z4*sinalf
      y4=x4*singam-y4*cosgam
      if (y4.gt.0.0) thg(i1)=-thg(i1)
      if (thg(i1).lt.-179.999999d0) thg(i1)=-179.999999d0
      if (thg(i1).gt.179.999999d0) thg(i1)=179.999999d0
      th2=thg(i1)*dgrrdn
********************************************************************** 
*                                                                    *
*     Calculate torsion angle energy                                 *
*                                                                    *
********************************************************************** 
      exbo1=abo(j(2))-valf(ia(j(2),1))
      exbo2=abo(j(3))-valf(ia(j(3),1))
      htovt=exbo1+exbo2
      expov=exp(vpar(26)*htovt)
      expov2=exp(-vpar(25)*(htovt))
      htov1=2.0+expov2
      htov2=1.0+expov+expov2
      etboadj=htov1/htov2

      btb2=bopi(lb)-1.0+etboadj
      bo2t=1.0-btb2
      bo2p=bo2t*bo2t
      bocor2=exp(v4(ity)*bo2p)

      hsin=sinhd*sinhe
      ethhulp=0.50*v1(ity)*(1.0+arg)+v2(ity)*bocor2*(1.0-arg2)+
     $v3(ity)*(0.50+2.0*arg2*arg-1.50*arg)

      exphua=exp(-vpar(24)*boa)
      exphub=exp(-vpar(24)*bob)
      exphuc=exp(-vpar(24)*boc)
      bocor4=(1.0-exphua)*(1.0-exphub)*(1.0-exphuc)
      eth=hsin*ethhulp*bocor4

      detdar=hsin*bocor4*(0.50*v1(ity)-2.0*v2(ity)*bocor2*arg+
     $v3(ity)*(6.0*arg2-1.5d0))
      detdhd=coshd*sinhe*bocor4*ethhulp
      detdhe=sinhd*coshe*bocor4*ethhulp

      detdboa=vpar(24)*exphua*(1.0-exphub)*(1.0-exphuc)*ethhulp*hsin
      detdbopib=-bocor4*2.0*v4(ity)*v2(ity)*
     $bo2t*bocor2*(1.0-arg2)*hsin
      detdbob=vpar(24)*exphub*(1.0-exphua)*
     $(1.0-exphuc)*ethhulp*hsin
      detdboc=vpar(24)*exphuc*(1.0-exphua)*
     $(1.0-exphub)*ethhulp*hsin

      detdsbo1=-(detdbopib)*
     $(vpar(25)*expov2/htov2+htov1*
     $(vpar(26)*expov-vpar(25)*expov2)/(htov2*htov2))

      et=et+eth
      estrain(j(2))=estrain(j(2))+0.50*eth !2nd atom energy
      estrain(j(3))=estrain(j(3))+0.50*eth !3rd atom energy

********************************************************************** 
*                                                                    *
*     Calculate conjugation energy                                   *
*                                                                    *
********************************************************************** 
      ba=(boa-1.50)*(boa-1.50)
      bb=(bob-1.50)*(bob-1.50)
      bc=(boc-1.50)*(boc-1.50)
      exphua1=exp(-vpar(28)*ba)
      exphub1=exp(-vpar(28)*bb)
      exphuc1=exp(-vpar(28)*bc)
      sbo=exphua1*exphub1*exphuc1
      dbohua=-2.0*(boa-1.50)*vpar(28)*exphua1*exphub1*exphuc1
      dbohub=-2.0*(bob-1.50)*vpar(28)*exphua1*exphub1*exphuc1
      dbohuc=-2.0*(boc-1.50)*vpar(28)*exphua1*exphub1*exphuc1
      arghu0=(arg2-1.0)*sinhd*sinhe
      ehulp=vconj(ity)*(arghu0+1.0)
      ecoh=ehulp*sbo
      decodar=sbo*vconj(ity)*2.0*arg*sinhd*sinhe
      decodbola=dbohua*ehulp
      decodbolb=dbohub*ehulp
      decodbolc=dbohuc*ehulp
      decodhd=coshd*sinhe*vconj(ity)*sbo*(arg2-1.0)
      decodhe=coshe*sinhd*vconj(ity)*sbo*(arg2-1.0)

      eco=eco+ecoh
      estrain(j(2))=estrain(j(2))+0.50*ecoh !2nd atom energy
      estrain(j(3))=estrain(j(3))+0.50*ecoh !3rd atom energy

    1 continue
********************************************************************** 
*                                                                    *
*     Calculate derivative torsion angle and conjugation energy      *
*     to cartesian coordinates                                       *
*                                                                    *
********************************************************************** 
      SINTH=SIN(THG(i1)*DGRRDN)
      IF (SINTH.GE.0.0.AND.SINTH.LT.1.0D-10) SINTH=1.0D-10
      IF (SINTH.LT.0.0.AND.SINTH.GT.-1.0D-10) SINTH=-1.0D-10
      IF (j(1).EQ.IB(LA,2)) THEN
      DO  K1=1,3
      DRADC(K1,1)=DRDC(LA,K1,1)
      DRADC(K1,2)=DRDC(LA,K1,2)
      end do
      ELSE
      DO  K1=1,3
      DRADC(K1,1)=DRDC(LA,K1,2)
      DRADC(K1,2)=DRDC(LA,K1,1)
      end do
      ENDIF
      IF (j(2).EQ.IB(LB,2)) THEN
      DO  K1=1,3
      DRBDC(K1,2)=DRDC(LB,K1,1)
      DRBDC(K1,3)=DRDC(LB,K1,2)
      end do
      ELSE
      DO K1=1,3
      DRBDC(K1,2)=DRDC(LB,K1,2)
      DRBDC(K1,3)=DRDC(LB,K1,1)
      end do
      ENDIF
      IF (j(3).EQ.IB(LC,2)) THEN
      DO K1=1,3
      DRCDC(K1,3)=DRDC(LC,K1,1)
      DRCDC(K1,4)=DRDC(LC,K1,2)
      end do
      ELSE
      DO K1=1,3
      DRCDC(K1,3)=DRDC(LC,K1,2)
      DRCDC(K1,4)=DRDC(LC,K1,1)
      end do
      ENDIF

      do k1=1,3
      dhddc(1,k1)=dh1rdc(1,k1)
      dhddc(2,k1)=dh1rdc(2,k1)
      dhddc(3,k1)=dh1rdc(3,k1)
      dhedc(1,k1+1)=dh2rdc(1,k1)
      dhedc(2,k1+1)=dh2rdc(2,k1)
      dhedc(3,k1+1)=dh2rdc(3,k1)
      end do

********************************************************************** 
*     write (64,*)j(1),j(2),j(3),j(4)
*     do k1=1,3
*     write (64,'(10f12.4)')(dh1rdc(k1,k2),k2=1,3),
*    $(dhdc(ld,k1,k2),k2=1,3),(dhddc(k1,k2),k2=1,4)
*     write (64,'(10f12.4)')(dh2rdc(k1,k2),k2=1,3),
*    $(dhdc(le,k1,k2),k2=1,3),(dhedc(k1,k2),k2=1,4)
*     end do
*     write (64,*)
********************************************************************** 
      HTRA=RLA+COSHD*(RLC*COSHE-RLB)
      HTRB=RLB-RLA*COSHD-RLC*COSHE
      HTRC=RLC+COSHE*(RLA*COSHD-RLB)
      HTHD=RLA*SINHD*(RLB-RLC*COSHE)
      HTHE=RLC*SINHE*(RLB-RLA*COSHD)
      HNRA=RLC*SINHD*SINHE
      HNRC=RLA*SINHD*SINHE
      HNHD=RLA*RLC*COSHD*SINHE
      HNHE=RLA*RLC*SINHD*COSHE

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      DO  K1=1,3
      DRDA(K1)=A(K1)/R4
      DO  K2=1,4
      DRVDC(K1,K2)=DRDA(K1)*DADC(K2)
      DTDC(K1,K2)=2.0*(DRADC(K1,K2)*HTRA+DRBDC(K1,K2)*HTRB+DRCDC(K1,K2
     $)*HTRC-DRVDC(K1,K2)*R4+DHDDC(K1,K2)*HTHD+DHEDC(K1,K2)*HTHE)
      DNDC(K1,K2)=2.0*(DRADC(K1,K2)*HNRA+DRCDC(K1,K2)*HNRC+DHDDC(K1,K2
     $)*HNHD+DHEDC(K1,K2)*HNHE)
      DARGTDC(i1,K1,K2)=(DTDC(K1,K2)-ARG*DNDC(K1,K2))/POEM
 
      ftmp = DARGTDC(i1,K1,K2)*detdar+
     $dargtdc(i1,k1,k2)*decodar+(detdhd+decodhd)*dhddc(k1,k2)+
     $(detdhe+decodhe)*dhedc(k1,k2)
      D(K1,J(K2))=D(K1,J(K2))+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(j(k2),k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/4
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do k2=1,4
                  ihu=j(k2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(la)
      ihu=idbo(la,i2)
      do k1=1,3
      ftmp = dbondc(la,k1,i2)*(detdboa+decodbola)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(la)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(la)
                  ihu=idbo(la,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(lb)
      ihu=idbo(lb,i2)
      do k1=1,3
         ftmp = dbondc(lb,k1,i2)*(detdbob+decodbolb)
     $   +dbopindc(lb,k1,i2)*detdbopib
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(lb)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(lb)
                  ihu=idbo(lb,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i2=1,idbo1(lc)
      ihu=idbo(lc,i2)
      do k1=1,3
      ftmp = dbondc(lc,k1,i2)*(detdboc+decodbolc)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(lc)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i2=1,idbo1(lc)
                  ihu=idbo(lc,i2)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      do i2=1,ia(j(2),2)
      iob=ia(j(2),2+i2)
      ncubo=nubon2(j(2),i2)
      if (bo(ncubo).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(ncubo)
      ihu=idbo(ncubo,i3)
      do k1=1,3
      ftmp = detdsbo1*dbondc(ncubo,k1,i3)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do
      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do

      do i2=1,ia(j(3),2)
      iob=ia(j(3),2+i2)
      ncubo=nubon2(j(3),i2)
      if (bo(ncubo).gt.0.0) then

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do i3=1,idbo1(ncubo)
      ihu=idbo(ncubo,i3)
      do k1=1,3
      ftmp = detdsbo1*dbondc(ncubo,k1,i3)
      d(k1,ihu)=d(k1,ihu)+ftmp

      if (Lvirial.eq.1) then
      do k1p=1,3
      virial_tmp(k1,k1p)=virial_tmp(k1,k1p)+ftmp*c(ihu,k1p)
      end do
      endif

      end do
      end do

      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/idbo1(ncubo)
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               do i3=1,idbo1(ncubo)
                  ihu=idbo(ncubo,i3)
                  atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               end do
            end do
         endif

      endif

      endif
      end do

   10 continue 

      return
      end 
********************************************************************** 
********************************************************************** 

      subroutine nonbon

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkch.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkdcell.blk"
#include "cbkenergies.blk"
#include "cbkff.blk"
#include "cbkia.blk"
#include "cbknonbon.blk"
#include "cbkpairs.blk"
#include "cbknvlown.blk"
#include "cellcoord.blk"
#include "control.blk"
#include "small.blk"

      dimension a(3),da(6)
      dimension virial_tmp(3,3),virialsym(6)
********************************************************************** 
*                                                                    *
*     Calculate vdWaals and Coulomb energies and derivatives         *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In nonbon'
c$$$      call timer(65)
c$$$      end if

      ew=0.0
      ep=0.0

      c1c=332.0638
      third=one/three
      fothird=4.0/3.0
      twothird=2.0/3.0
      h15=(vpar(29)-1.0)/vpar(29)

      nptmp=0
      nstmp=0
      do 10 ivl=1,nvpair-nvlself
c Use precomputed midpoint criterion to decide if interaction is owned.
      if (nvlown(ivl).eq.1) then

      i1=nvl1(ivl)
      i2=nvl2(ivl)

      call dista2(i1,i2,rr,a(1),a(2),a(3))
      if (rr.gt.swb.or.rr.lt.0.001) goto 10

      ity1=ia(i1,1)
      ity2=ia(i2,1)
      imol1=iag(i1,3+mbond)
      imol2=iag(i2,3+mbond)
      rr2=rr*rr

      sw=1.0
      sw1=0.0
      call taper(rr,rr2)
**********************************************************************
*                                                                    *
*     Calculate vdWaals energy                                       *
*                                                                    *
**********************************************************************
      p1=p1co(ity1,ity2)
      p2=p2co(ity1,ity2)
      p3=p3co(ity1,ity2)
      hulpw=(rr**vpar(29)+gamwco(ity1,ity2))
      rrw=hulpw**(1.0/vpar(29))
      h1=exp(p3*(1.0-rrw/p1))
      h2=exp(0.50*p3*(1.0-rrw/p1))

      ewh=p2*(h1-2.0*h2)
      rrhuw=rr**(vpar(29)-1.0)
      dewdr=(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**(-h15))

**********************************************************************
*                                                                    *
*     Calculate Coulomb energy                                       *
*                                                                    *
**********************************************************************
      q1q2=ch(i1)*ch(i2)
      hulp1=(rr2*rr+gamcco(ity1,ity2))
      eph=c1c*q1q2/(hulp1**third)
      depdr=-c1c*q1q2*rr2/(hulp1**fothird)
**********************************************************************
*                                                                    *
*     Taper correction                                               *
*                                                                    *
**********************************************************************
      ephtap=eph*sw
      depdrtap=depdr*sw+eph*sw1
      ewhtap=ewh*sw
      dewdrtap=dewdr*sw+ewh*sw1

*     write (64,*)i1,i2,p1,p2,p3,gamwco(ity1,ity2),vpar(29),ewh,ew
      ew=ew+ewhtap
      ep=ep+ephtap
      estrain(i1)=estrain(i1)+0.50*(ewhtap+ephtap) !1st atom energy
      estrain(i2)=estrain(i2)+0.50*(ewhtap+ephtap) !2nd atom energy

**********************************************************************
*                                                                    *
*     Calculate derivatives vdWaals energy to cartesian              *
*     coordinates                                                    *
*                                                                    *
**********************************************************************

      if (Lvirial.eq.1) then
         do k1=1,3
         do k2=1,3
            virial_tmp(k1,k2) = 0.0
         end do
         end do
      endif

      do k4=1,3
         ftmp = (dewdrtap+depdrtap)*(a(k4)/rr)
         d(k4,i1)=d(k4,i1)+ftmp
         d(k4,i2)=d(k4,i2)-ftmp
         if (Lvirial.eq.1) then
            do k1p=1,3
               virial_tmp(k4,k1p)=virial_tmp(k4,k1p)+
     $              ftmp*c(i1,k1p)-ftmp*c(i2,k1p)
            end do
         endif
      end do

      if (Lvirial.eq.1) then
         virialsym(1) = virial_tmp(1,1)
         virialsym(2) = virial_tmp(2,2)
         virialsym(3) = virial_tmp(3,3)
         virialsym(4) = virial_tmp(1,2)
         virialsym(5) = virial_tmp(1,3)
         virialsym(6) = virial_tmp(2,3)
         do k1 = 1,6
            virial(k1) = virial(k1) + virialsym(k1)
         end do

         if (Latomvirial.eq.1) then
            frac = 1.0d0/2
            do k1 = 1,6
               vtmp = virialsym(k1)*frac
               ihu=i1
               atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
               ihu=i2
               atomvirial(k1,ihu) = atomvirial(k1,ihu) + vtmp
            end do
         endif

      endif

      endif

   10 continue

      return
      end

********************************************************************** 
**********************************************************************
 
      subroutine efield
 
**********************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkch.blk"
#include "cbkcha.blk"
#include "cbkd.blk"
#include "cbkefield.blk"
#include "cbkenergies.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "small.blk"
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In efield'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
**********************************************************************
*                                                                    *
*     Electric field                                                 *
*                                                                    *
**********************************************************************
      efi=0.0
      efix=0.0
      efiy=0.0
      efiz=0.0
      c1c=332.0638       !Coulomb energy conversion
      c1=23.02           !conversion from kcal to eV
 
      if (ifieldx.eq.1) then
      do i1=1,na
      efih=vfieldx*c1*c1c*ch(i1)*c(i1,1)
      efix=efix+efih
      estrain(i1)=estrain(i1)+efih               !atom energy

      defidc=c1*c1c*vfieldx*ch(i1)
      d(1,i1)=d(1,i1)+defidc
      end do
      end if
 
      if (ifieldy.eq.1) then
      do i1=1,na
      efih=vfieldy*c1*c1c*ch(i1)*c(i1,2)
      efiy=efiy+efih
      estrain(i1)=estrain(i1)+efih               !atom energy

      defidc=c1*c1c*vfieldy*ch(i1)
      d(2,i1)=d(2,i1)+defidc
      end do
      end if
 
      if (ifieldz.eq.1) then
      do i1=1,na
      efih=vfieldz*c1*c1c*ch(i1)*c(i1,3)
      efiz=efiz+efih
      estrain(i1)=estrain(i1)+efih               !atom energy

      defidc=c1*c1c*vfieldz*ch(i1)
      d(3,i1)=d(3,i1)+defidc
      end do
      end if
 
      efi=efix+efiy+efiz
      return
      end
**********************************************************************
********************************************************************** 

      subroutine restraint

********************************************************************** 
#include "cbka.blk"
#include "cbkatomcoord.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkenergies.blk"
#include "cbkrestr.blk"
#include "cbktorang.blk"
#include "cbktorsion.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "small.blk"
#include "cbkinit.blk"
      dimension drda(3),j(4),dhrdc(3,3),dargdc(3,3)
********************************************************************** 
*                                                                    *
*     Calculate distance restraint energy                            *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In restraint'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      do i1=1,nrestra
      ih1=irstra(i1,1)
      ih2=irstra(i1,2)
      if (itend(i1).eq.0.or.(mdstep.gt.itstart(i1).and.mdstep.lt.
     $itend(i1))) then
      call dista2(ih1,ih2,rr,dx,dy,dz)
      diffr=rr-rrstra(i1)
*     diffr=rrstra(i1)
      exphu=exp(-vkrst2(i1)*(diffr*diffr))
      erh=vkrstr(i1)*(1.0-exphu)
      deresdr=2.0*vkrst2(i1)*diffr*vkrstr(i1)*exphu
*     deresdr=-2.0*vkrst2(i1)*diffr*vkrstr(i1)*exphu
      eres=eres+erh
      drda(1)=dx/rr
      drda(2)=dy/rr
      drda(3)=dz/rr
      do k1=1,3
      d(k1,ih1)=d(k1,ih1)+deresdr*drda(k1)
      d(k1,ih2)=d(k1,ih2)-deresdr*drda(k1)
      end do
      end if
      end do
      
**********************************************************************
*                                                                    *
*     Calculate angle restraint energy                               *
*                                                                    *
**********************************************************************
      do i1=1,nrestrav
      j(1)=irstrav(i1,1)
      j(2)=irstrav(i1,2)
      j(3)=irstrav(i1,3)
      ittr=0
*     do i2=1,nval
*     if (j(1).eq.iv(i2,2).and.j(2).eq.iv(i2,3).and.j(3).eq.iv(i2,4))
*    $ittr=i2
*     end do
*     if (ittr.eq.0) stop 'Wrong valence angle restraint'
      call calvalres(j(1),j(2),j(3),arg,hr,dhrdc,dargdc)
      vaval=hr*rdndgr
      diffv=-(vaval-vrstra(i1))*dgrrdn
      exphu=exp(-vkr2v(i1)*(diffv*diffv))
      erh=vkrv(i1)*(1.0-exphu)
      deresdv=-2.0*vkr2v(i1)*diffv*vkrv(i1)*exphu
      eres=eres+erh
      do k1=1,3
      do k2=1,3
      d(k1,j(k2))=d(k1,j(k2))+deresdv*dhrdc(k1,k2)
      end do
      end do

      end do

********************************************************************** 
*                                                                    *
*     Calculate torsion restraint energy                             *
*                                                                    *
********************************************************************** 
      do i1=1,nrestrat
      j(1)=irstrat(i1,1)
      j(2)=irstrat(i1,2)
      j(3)=irstrat(i1,3)
      j(4)=irstrat(i1,4)
      ittr=0
      do i2=1,ntor
      if (j(1).eq.it(i2,2).and.j(2).eq.it(i2,3).and.j(3).eq.it(i2,4)
     $.and.j(4).eq.it(i2,5)) ittr=i2
      if (j(4).eq.it(i2,2).and.j(3).eq.it(i2,3).and.j(2).eq.it(i2,4)
     $.and.j(1).eq.it(i2,5)) ittr=i2
      end do
      if (ittr.eq.0) then
      write (*,*)'Wrong torsion restraint'
      write (*,*)i1,j(1),j(2),j(3),j(4)
      stop 'Wrong torsion restraint'
      end if
      vtor=thg(ittr)
      difft=-(vtor-trstra(i1))*dgrrdn
      exphu=exp(-vkr2t(i1)*(difft*difft))
      erh=vkrt(i1)*(1.0-exphu)
      deresdt=2.0*vkr2t(i1)*difft*vkrt(i1)*exphu
      if (vtor.lt.zero) deresdt=-deresdt
      eres=eres+erh
      do k1=1,3
      do k2=1,4
      d(k1,j(k2))=d(k1,j(k2))+deresdt*dargtdc(ittr,k1,k2)
      end do
      end do

      end do
********************************************************************** 
*                                                                    *
*     Calculate mass centre restraint energy                         *
*                                                                    *
********************************************************************** 
      do i1=1,nrestram
      j1=irstram(i1,2)
      j2=irstram(i1,3)
      j3=irstram(i1,4)
      j4=irstram(i1,5)
      kdir=irstram(i1,1)
      cmx1=0.0 
      cmy1=0.0 
      cmz1=0.0 
      cmx2=0.0 
      cmy2=0.0 
      cmz2=0.0 
      summas1=0.0
      summas2=0.0
      do i2=j1,j2
      cmx1=cmx1+c(i2,1)*xmasat(i2)
      cmy1=cmy1+c(i2,2)*xmasat(i2)
      cmz1=cmz1+c(i2,3)*xmasat(i2)
      summas1=summas1+xmasat(i2)
      end do
      cmx1=cmx1/summas1
      cmy1=cmy1/summas1
      cmz1=cmz1/summas1
      if (mdstep.lt.2) then
      rmstrax(i1)=cmx1
      rmstray(i1)=cmy1
      rmstraz(i1)=cmz1
      end if
      if (kdir.le.3) then
      do i2=j3,j4
      cmx2=cmx2+c(i2,1)*xmasat(i2)
      cmy2=cmy2+c(i2,2)*xmasat(i2)
      cmz2=cmz2+c(i2,3)*xmasat(i2)
      summas2=summas2+xmasat(i2)
      end do
      cmx2=cmx2/summas2
      cmy2=cmy2/summas2
      cmz2=cmz2/summas2
      end if
      if (kdir.eq.1) dist=cmx1-cmx2
      if (kdir.eq.2) dist=cmy1-cmy2
      if (kdir.eq.3) dist=cmz1-cmz2
      if (kdir.eq.4) then
      distx=cmx1-rmstrax(i1)
      disty=cmy1-rmstray(i1)
      distz=cmz1-rmstraz(i1)
      dist=sqrt(distx*distx+disty*disty+distz*distz)
      end if
      dismacen(i1)=dist
      dist=dist-rmstra1(i1)
      erh=rmstra2(i1)*dist*dist
      deresdr=2.0*dist*rmstra2(i1)
*     exphu=exp(-rmstra3(i1)*(dist*dist))
*     erh=rmstra2(i1)*(1.0-exphu)
*     deresdr=2.0*rmstra3(i1)*dist*rmstra2(i1)*exphu
      eres=eres+erh
      if (kdir.le.3) then
      do i2=j1,j2
      d(kdir,i2)=d(kdir,i2)+deresdr*xmasat(i2)/summas1
      end do
      do i2=j3,j4
      d(kdir,i2)=d(kdir,i2)-deresdr*xmasat(i2)/summas2
      end do
      end if
      if (kdir.eq.4.and.mdstep.gt.5) then
      do i2=j1,j2
      d(1,i2)=d(1,i2)+deresdr*(distx/dist)*(xmasat(i2)/summas1)
      d(2,i2)=d(2,i2)+deresdr*(disty/dist)*(xmasat(i2)/summas1)
      d(3,i2)=d(3,i2)+deresdr*(distz/dist)*(xmasat(i2)/summas1)
      end do
      end if
      end do
********************************************************************** 
*                                                                    *
*     Calculate morphing energy                                      *
*                                                                    *
********************************************************************** 
      if (imorph.eq.1) then
      distot=zero
      do i1=1,na
      dmx=c(i1,1)-cmo(i1,1)
      dmy=c(i1,2)-cmo(i1,2)
      dmz=c(i1,3)-cmo(i1,3)
      dism=sqrt(dmx*dmx+dmy*dmy+dmz*dmz)
      distot=distot+dism
*     exphu=exp(-vmo2(i1)*(dism*dism))
*     erh=vmo1(i1)*(1.0-exphu)
      erh=vmo1(i1)*dism
      eres=eres+erh
*     deresddis=2.0*vmo2(i1)*dism*vmo1(i1)*exphu
      deresddis=vmo1(i1)
      drda1=dmx/dism
      drda2=dmy/dism
      drda3=dmz/dism
      d(1,i1)=d(1,i1)+deresddis*drda1
      d(2,i1)=d(2,i1)+deresddis*drda2
      d(3,i1)=d(3,i1)+deresddis*drda3
      end do
  
      end if
      
      
      return
      end
********************************************************************** 
******************************************************************** 

      subroutine calvalres (ja1,ja2,ja3,arg,hr,dhrdc,dargdc)

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
      dimension a(3),b(3),j(3),dradc(3,3),drbdc(3,3),dtdc(3,3),
     $dargdc(3,3),dndc(3,3),dadc(3),dbdc(3),dhrdc(3,3)
********************************************************************** 
*                                                                    *
*     Calculate valency angles and their derivatives to cartesian    *
*     coordinates  for restraint calculations                        *
*                                                                    *
********************************************************************** 
c$$$*     if (ndebug.eq.1) then
c$$$C*     open (65,file='fort.65',status='unknown',access='append')
c$$$*     write (65,*) 'In calvalres'
c$$$*     call timer(65)
c$$$*     close (65)
c$$$*     end if

      dadc(1)=-1.0
      dadc(2)=1.0
      dadc(3)=0.0
      dbdc(1)=0.0
      dbdc(2)=1.0
      dbdc(3)=-1.0
      do k1=1,3
      do k2=1,3
      dradc(k1,k2)=0.0
      drbdc(k1,k2)=0.0
      end do
      end do
**********************************************************************
*                                                                    *
*     Determine valency angle                                        *
*                                                                    *
**********************************************************************
      call dista2(ja1,ja2,rla,dx1,dy1,dz1)
      call dista2(ja2,ja3,rlb,dx2,dy2,dz2)
 
      a(1)=-dx1
      a(2)=-dy1
      a(3)=-dz1
      b(1)=dx2
      b(2)=dy2
      b(3)=dz2
      poem=rla*rlb
      tel=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      arg=tel/poem
      arg2=arg*arg
      s1ma22=1.0-arg2
      if (s1ma22.lt.1.0d-10) s1ma22=1.0d-10
      s1ma2=sqrt(s1ma22)
      if (arg.gt.1.0) arg=1.0
      if (arg.lt.-1.0) arg=-1.0
      hr=acos(arg)
**********************************************************************
*                                                                    *
*     Calculate derivative valency angle to cartesian coordinates    *
*                                                                    *
**********************************************************************
      do k1=1,3
      dradc(k1,1)=-a(k1)/rla
      dradc(k1,2)=a(k1)/rla
      end do

      do k1=1,3
      drbdc(k1,2)=b(k1)/rlb
      drbdc(k1,3)=-b(k1)/rlb
      end do

      do k1=1,3
      do k2=1,3
      dndc(k1,k2)=rla*drbdc(k1,k2)+rlb*dradc(k1,k2)
      dtdc(k1,k2)=a(k1)*dbdc(k2)+b(k1)*dadc(k2)
      dargdc(k1,k2)=(dtdc(k1,k2)-arg*dndc(k1,k2))/poem
      dhrdc(k1,k2)=-dargdc(k1,k2)/s1ma2
      end do
      end do
      
   10 continue

      return
      end
********************************************************************** 
******************************************************************** 

      subroutine calvalhb (ja1,ja2,ja3,ix,iy,iz,arg,hr,dhrdc,dargdc)

********************************************************************** 
#include "cbka.blk"
#include "cbkc.blk"
      dimension a(3),b(3),j(3),dradc(3,3),drbdc(3,3),dtdc(3,3),
     $dargdc(3,3),dndc(3,3),dadc(3),dbdc(3),dhrdc(3,3)
********************************************************************** 
*                                                                    *
*     Calculate valency angles and their derivatives to cartesian    *
*     coordinates  for hydrogen bond calculations                    *
*                                                                    *
********************************************************************** 
c$$$*     if (ndebug.eq.1) then
c$$$*     open (65,file='fort.65',status='unknown',access='append')
c$$$*     write (65,*) 'In calvalhb'
c$$$*     call timer(65)
c$$$*     close (65)
c$$$*     end if

      dadc(1)=-1.0
      dadc(2)=1.0
      dadc(3)=0.0
      dbdc(1)=0.0
      dbdc(2)=1.0
      dbdc(3)=-1.0
      do k1=1,3
      do k2=1,3
      dradc(k1,k2)=0.0
      drbdc(k1,k2)=0.0
      end do
      end do
**********************************************************************
*                                                                    *
*     Determine valency angle                                        *
*                                                                    *
**********************************************************************
      call dista2(ja1,ja2,rla,dx1,dy1,dz1)
      call dista2(ja2,ja3,rlb,dx2,dy2,dz2)

      a(1)=-dx1
      a(2)=-dy1
      a(3)=-dz1
      b(1)=dx2
      b(2)=dy2
      b(3)=dz2
      poem=rla*rlb
      tel=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      arg=tel/poem
      arg2=arg*arg
      s1ma22=1.0-arg2
      if (s1ma22.lt.1.0d-10) s1ma22=1.0d-10
      s1ma2=sqrt(s1ma22)
      if (arg.gt.1.0) arg=1.0
      if (arg.lt.-1.0) arg=-1.0
      hr=acos(arg)
**********************************************************************
*                                                                    *
*     Calculate derivative valency angle to cartesian coordinates    *
*                                                                    *
**********************************************************************
      do k1=1,3
      dradc(k1,1)=-a(k1)/rla
      dradc(k1,2)=a(k1)/rla
      end do

      do k1=1,3
      drbdc(k1,2)=b(k1)/rlb
      drbdc(k1,3)=-b(k1)/rlb
      end do

      do k1=1,3
      do k2=1,3
      dndc(k1,k2)=rla*drbdc(k1,k2)+rlb*dradc(k1,k2)
      dtdc(k1,k2)=a(k1)*dbdc(k2)+b(k1)*dadc(k2)
      dargdc(k1,k2)=(dtdc(k1,k2)-arg*dndc(k1,k2))/poem
      dhrdc(k1,k2)=-dargdc(k1,k2)/s1ma2
      end do
      end do
      
   10 continue

      return
      end
********************************************************************** 
********************************************************************** 

      subroutine caltor(ja1,ja2,ja3,ja4,ht)

********************************************************************** 
#include "cbka.blk"
#include "cbkenergies.blk"
#include "cbktregime.blk"
#include "control.blk"
#include "cbkinit.blk"
      DIMENSION  A(3),DRDA(3),DADC(4),DRADC(3,4),DRBDC(3,4),
     $DRCDC(3,4),DHDDC(3,4),DHEDC(3,4),DRVDC(3,4),DTDC(3,4),
     $DNDC(3,4)
      dimension j(4),dvdc1(3,3),dargdc1(3,3),dvdc2(3,3),dargdc2(3,3)
********************************************************************** 
*                                                                    *
*     Calculate torsion angle (for internal coordinates output)      *
*                                                                    *
********************************************************************** 
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In caltor'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      do k1=1,3
      do k2=1,4
      dhddc(k1,k2)=0.0
      dhedc(k1,k2)=0.0
      dradc(k1,k2)=0.0
      drbdc(k1,k2)=0.0
      drcdc(k1,k2)=0.0
      end do
      end do
      et=0.0
      eco=0.0
      dadc(1)=1.0
      dadc(2)=0.0
      dadc(3)=0.0
      dadc(4)=-1.0
      call dista2(ja1,ja2,rla,dx1,dy1,dz1)
      call dista2(ja2,ja3,rlb,dx2,dy2,dz2)
      call dista2(ja3,ja4,rlc,dx2,dy2,dz2)
      call dista2(ja1,ja4,r4,dx2,dy2,dz2)
      call calvalres(ja1,ja2,ja3,arg1,h1,dvdc1,dargdc1)
      call calvalres(ja2,ja3,ja4,arg2,h2,dvdc2,dargdc2)
********************************************************************** 
*                                                                    *
*     Determine torsion angle                                        *
*                                                                    *
********************************************************************** 
      d142=r4*r4
      coshd=cos(h1)
      coshe=cos(h2)
      sinhd=sin(h1)
      sinhe=sin(h2)
      poem=2.0*rla*rlc*sinhd*sinhe
      poem2=poem*poem
      tel=rla*rla+rlb*rlb+rlc*rlc-d142-2.0*(rla*rlb*coshd-rla*rlc*
     $coshd*coshe+rlb*rlc*coshe)
      arg=tel/poem
      if (arg.gt.1.0) arg=1.0
      if (arg.lt.-1.0) arg=-1.0
      arg2=arg*arg
      ht=acos(arg)*rdndgr
      k1=ja1
      k2=ja2
      k3=ja3
      k4=ja4
      call dista2(k3,k2,dis,x3,y3,z3)
      y32z32=y3*y3+z3*z3
      wort1=sqrt(y32z32)+1e-6
      wort2=sqrt(y32z32+x3*x3)+1e-6
      sinalf=y3/wort1
      cosalf=z3/wort1
      sinbet=x3/wort2
      cosbet=wort1/wort2
      call dista2(k1,k2,dis,x1,y1,z1)
      x1=x1*cosbet-y1*sinalf*sinbet-z1*cosalf*sinbet
      y1=y1*cosalf-z1*sinalf
      wort3=sqrt(x1*x1+y1*y1)+1e-6
      singam=y1/wort3
      cosgam=x1/wort3
      call dista2(k4,k2,dis,x4,y4,z4)
      x4=x4*cosbet-y4*sinalf*sinbet-z4*cosalf*sinbet
      y4=y4*cosalf-z4*sinalf
      y4=x4*singam-y4*cosgam
      if (y4.gt.0.0) ht=-ht
      if (ht.lt.-179.999999d0) ht=-179.999999d0
      if (ht.gt.179.999999d0) ht=179.999999d0

      return
      end 
********************************************************************** 

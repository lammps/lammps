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
************************************************************************

      subroutine taper(r,r2)

************************************************************************
#include "cbka.blk"
#include "cbkconst.blk"
#include "cbkenergies.blk"
#include "cbkinit.blk"
#include "cbknonbon.blk"
************************************************************************
*                                                                      *
*     Taper function for Coulomb interaction                           *
*                                                                      *
************************************************************************
      r3=r2*r
      SW=SWC7*R3*R3*R+SWC6*R3*R3+SWC5*R3*R2+SWC4*R2*R2+SWC3*R3+SWC2*R2+
     $SWC1*R+SWC0
      SW1=7.0D0*SWC7*R3*R3+6.0D0*SWC6*R3*R2+5.0D0*SWC5*R2*R2+
     $4.0D0*SWC4*R3+THREE*SWC3*R2+TWO*SWC2*R+SWC1
      return
      end
************************************************************************
************************************************************************

      subroutine tap7th

************************************************************************
#include "cbka.blk"
#include "cbkc.blk"
#include "cbkconst.blk"
#include "cbkd.blk"
#include "cbkff.blk"
#include "cbkh.blk"
#include "control.blk"

************************************************************************
*                                                                      *
*     7th order taper function setup                                   *
*                                                                      *
************************************************************************
c$$$      if (ndebug.eq.1) then
c$$$C      open (65,file='fort.65',status='unknown',access='append')
c$$$      write (65,*) 'In tap7th'
c$$$      call timer(65)
c$$$      close (65)
c$$$      end if
      D1=SWB-SWA
      D7=D1**7.0D0
      SWA2=SWA*SWA
      SWA3=SWA2*SWA
      SWB2=SWB*SWB
      SWB3=SWB2*SWB
 
************************************************************************
*     7th order taper function                                         *
************************************************************************
 
      SWC7=  20.0D0/D7
      SWC6= -70.0D0*(SWA+SWB)/D7
      SWC5=  84.0D0*(SWA2+3.0D0*SWA*SWB+SWB2)/D7
      SWC4= -35.0D0*(SWA3+9.0D0*SWA2*SWB+9.0D0*SWA*SWB2+SWB3)/D7
      SWC3= 140.0D0*(SWA3*SWB+3.0D0*SWA2*SWB2+SWA*SWB3)/D7
      SWC2=-210.0D0*(SWA3*SWB2+SWA2*SWB3)/D7
      SWC1= 140.0D0*SWA3*SWB3/D7
      SWC0=(-35.0D0*SWA3*SWB2*SWB2+21.0D0*SWA2*SWB3*SWB2+
     $7.0D0*SWA*SWB3*SWB3+SWB3*SWB3*SWB)/D7
 
      return
      END

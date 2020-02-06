! ------------ ----------------------------------------------------------
!   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
!   http://lammps.sandia.gov, Sandia National Laboratories
!   Steve Plimpton, sjplimp@sandia.gov
!
!   Copyright (2003) Sandia Corporation.  Under the terms of Contract
!   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
!   certain rights in this software.  This software is distributed under
!   the GNU General Public License.
!
!   See the README file in the top-level LAMMPS directory.
!   
!   Contributing author: Alexey N. Volkov, UA, avolkov1@ua.edu
!------------------------------------------------------------------------- 

module TPMM0 !**************************************************************************************
!
! TMD Library: Combined/Weighted potential of type 0
!
! Direct application of SST potential to calculation of segment-segment interaction
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!***************************************************************************************************

!use TMDCounters
use TubePotMono

implicit none

contains !******************************************************************************************
        
        integer*4 function TPMInteractionFSS ( Q, U, F1_1, F1_2, F2_1, F2_2, R1_1, R1_2, R2_1, R2_2, EType )
        real*8, intent(inout)                   :: Q, U
        real*8, dimension(0:2), intent(inout)   :: F1_1, F1_2, F2_1, F2_2
        real*8, dimension(0:2), intent(in)      :: R1_1, R1_2, R2_1, R2_2
        integer*4, intent(in)                   :: EType
        !-------------------------------------------------------------------------------------------
        real*8                                  :: Qa, Ua, Fd, L2
        real*8, dimension(0:2)                  :: F1_1a, F1_2a, F2_1a, F2_2a, R2_3, R2, Laxis2, F
	integer*4				:: IntSign
        !-------------------------------------------------------------------------------------------
!                C_TPM_4 = C_TPM_4 + 1
                R2 = 0.5d+00 * ( R2_1 + R2_2 )
                Laxis2 = R2_2 - R2_1
                L2 = S_V3norm3 ( Laxis2 )
                Laxis2 = Laxis2 / L2
                if ( EType < 2 ) then
                        TPMInteractionFSS = TPMInteractionF ( Q, U, F1_1, F1_2, F2_1, F2_2, Fd, R1_1, R1_2, R2_1, R2_2, 1 ) 
                        R2_3 = R2_2 + R2_2 - R2_1
                        IntSign = TPMInteractionF ( Qa, Ua, F1_1a, F1_2a, F2_1a, F2_2a, Fd, R1_1, R1_2, R2_2, R2_3, 1 ) 
			if ( IntSign > 0 ) then
				TPMInteractionFSS = 1
	                        call TPMSegmentForces ( F2_1a, F2_2a, F1_1a, F1_2a, R1_1, R1_2, R2, Laxis2, L2 )
        	                F = ( Fd - S_V3xV3 ( F2_2a, Laxis2 ) ) * Laxis2
                	        F2_2a = F2_2a + F
                        	F2_1a = F2_1a - F
			end if
                else
                        TPMInteractionFSS = TPMInteractionF ( Q, U, F1_1, F1_2, F2_1, F2_2, Fd, R1_1, R1_2, R2_1, R2_2, 2 ) 
                        R2_3 = R2_1 + R2_1 - R2_2
                        IntSign = TPMInteractionF ( Qa, Ua, F1_1a, F1_2a, F2_1a, F2_2a, Fd, R1_1, R1_2, R2_1, R2_3, 1 )
			if ( IntSign > 0 ) then
				TPMInteractionFSS = 1
	                        call TPMSegmentForces ( F2_1a, F2_2a, F1_1a, F1_2a, R1_1, R1_2, R2, Laxis2, L2 )
        	                F = ( - Fd - S_V3xV3 ( F2_1a, Laxis2 ) ) * Laxis2
                	        F2_1a = F2_1a + F
                        	F2_2a = F2_2a - F
			end if
                end if
		if ( IntSign > 0 ) then
	                Q = Q - Qa
                        if ( Q < 0.0d+00 ) Q = 0.0d+00
	                U = U - Ua
        	        F2_1 = F2_1 - F2_1a 
                	F2_2 = F2_2 - F2_2a 
                	F1_1 = F1_1 - F1_1a
	                F1_2 = F1_2 - F1_2a
		end if
        end function TPMInteractionFSS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMInteractionFW0 ( QQ, U, U1, U2, UU, F1, F2, F, G1, G2, R1, R2, N, NMAX, R )
        real*8, intent(inout)                           :: U, U1, U2
        integer*4, intent(in)                           :: N, NMAX
        real*8, dimension(0:NMAX-1), intent(out)        :: QQ, UU
        real*8, dimension(0:2), intent(out)             :: F1, F2
        real*8, dimension(0:2,0:NMAX-1), intent(out)    :: F, G1, G2
        real*8, dimension(0:2), intent(in)              :: R1, R2
        real*8, dimension(0:2,0:NMAX-1), intent(in)     :: R
        !-------------------------------------------------------------------------------------------
        integer*4                                       :: i, SType2, GeomID, EType
        real*8                                          :: Ua
        real*8, dimension(0:2)                          :: F1_1a, F1_2a, F2_1a, F2_2a
        real*8, dimension(0:2)                          :: R1a, R2a, Laxis1, Laxis2, L12, DR
        real*8                                          :: L1, L2, D1, D2, H, cosA, D, Dmina, Dminb
        !-------------------------------------------------------------------------------------------
                QQ = 0.0d+00
                U  = 0.0d+00
                U1 = 0.0d+00
                U2 = 0.0d+00
                UU = 0.0d+00
                F1 = 0.0d+00
                F2 = 0.0d+00
                F  = 0.0d+00
                G1 = 0.0d+00
                G2 = 0.0d+00
                TPMInteractionFW0 = 0
                do i = 0, N - 2
                        R1a = 0.5d+00 * ( R1 + R2 )
                        R2a = 0.5d+00 * ( R(0:2,i+1) + R(0:2,i) )
                        Laxis1 = R2 - R1
                        Laxis2 = R(0:2,i+1) - R(0:2,i)
                        L1 = S_V3norm3 ( Laxis1 )
                        L2 = S_V3norm3 ( Laxis2 )
                        Laxis1 = Laxis1 / L1
                        Laxis2 = Laxis2 / L2
                        L2 = 0.5d+00 * L2
                        L1 = 0.5d+00 * L1
                        GeomID = LineLine ( H, cosA, D1, D2, L12, R1a, Laxis1, R2a, Laxis2, TPGeomPrec )

                        DR = R1 - R(0:2,i)
                        call ApplyPeriodicBC ( DR )
                        Dmina = sqr ( DR(0) ) + sqr ( DR(1) ) + sqr ( DR(2) )
                        DR = R2 - R(0:2,i)
                        call ApplyPeriodicBC ( DR )
                        D = sqr ( DR(0) ) + sqr ( DR(1) ) + sqr ( DR(2) )
                        if ( D < Dmina ) Dmina = D
                        if ( GeomID == MD_LINES_NONPAR ) then
                                D = ( D2 - L2 ) * cosA
                                if ( D > D1 - L1 .and. D < D1 + L1 ) then
                                        D = sqr ( D2 - L2 ) * ( 1.0d+00 - sqr ( cosA ) ) + sqr ( H )
                                        if ( D < Dmina ) Dmina = D
                                end if
                        else
                                call LinePoint ( D, DR, R1, Laxis1, R(0:2,i) )
                                if ( D > 0.0d+00 .and. D < 2.0d+00 * L1 ) then
                                        DR = DR - R(0:2,i)
                                        call ApplyPeriodicBC ( DR )
                                        D = sqr ( DR(0) ) + sqr ( DR(1) ) + sqr ( DR(2) )
                                        if ( D < Dmina ) Dmina = D
                                end if
                        end if
                        
                        DR = R1 - R(0:2,i+1)
                        call ApplyPeriodicBC ( DR )
                        Dminb = sqr ( DR(0) ) + sqr ( DR(1) ) + sqr ( DR(2) )
                        DR = R2 - R(0:2,i+1)
                        call ApplyPeriodicBC ( DR )
                        D = sqr ( DR(0) ) + sqr ( DR(1) ) + sqr ( DR(2) )
                        if ( D < Dminb ) Dminb = D
                        if ( GeomID == MD_LINES_NONPAR ) then
                                D = ( D2 + L2 ) * cosA
                                if ( D > D1 - L1 .and. D < D1 + L1 ) then
                                        D = sqr ( D2 + L2 ) * ( 1.0d+00 - sqr ( cosA ) ) + sqr ( H )
                                        if ( D < Dminb ) Dminb = D
                                end if
                        else
                                call LinePoint ( D, DR, R1, Laxis1, R(0:2,i+1) )
                                if ( D > 0.0d+00 .and. D < 2.0d+00 * L1 ) then
                                        DR = DR - R(0:2,i+1)
                                        call ApplyPeriodicBC ( DR )
                                        D = sqr ( DR(0) ) + sqr ( DR(1) ) + sqr ( DR(2) )
                                        if ( D < Dminb ) Dminb = D
                                end if
                        end if
                        
                        if ( Dmina < Dminb ) then
                                EType = 1
                        else
                                EType = 2
                        end if
                        
                        if ( TPMInteractionFSS ( QQ(i), Ua, F1_1a, F1_2a, F2_1a, F2_2a, R1, R2, R(0:2,i), R(0:2,i+1), EType ) > 0 ) then
				TPMInteractionFW0 = 1
	                        U = U + Ua
        	                Ua = 0.25d+00 * Ua 
                	        U1 = U1 + Ua 
                        	U2 = U2 + Ua 
	                        UU(i) = UU(i) + Ua
        	                UU(i+1) = UU(i+1) + Ua
                        	F1 = F1 + F1_1a
                	        F2 = F2 + F1_2a
	                        F(0:2,i) = F(0:2,i) + F2_1a
        	                F(0:2,i+1) = F(0:2,i+1) + F2_2a
                	        G2(0:2,i) = F2_1a
                        	G1(0:2,i+1) = F2_2a
			end if
                end do
        end function TPMInteractionFW0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module TPMM0 !**********************************************************************************

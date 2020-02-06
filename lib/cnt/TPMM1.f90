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

module TPMM1 !**************************************************************************************
!
! TMD Library: Combined/Weighted potential of type 3
!
! Weighting functions are the same as in potential of type 2.
! Calculation of the combined potential is based on the 'extended' chain.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran.
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!***************************************************************************************************

!use TMDCounters
use TubePotMono

implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------

        ! Maximal length of a segment chain
        integer*4, parameter                            :: TPM_MAX_CHAIN = 100

!---------------------------------------------------------------------------------------------------
! Numerical parameters
!---------------------------------------------------------------------------------------------------

        ! Switching parameters
        real*8                                          :: TPMC123 = 1.0d+00    ! Non-dimensional
        real*8                                          :: TPMC3 = 10.0d+00     ! (A)

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------

        ! These global variables are used to speedup calculations
        real*8, dimension(0:2,0:TPM_MAX_CHAIN-1)        :: E1, E2, EE1, EE2
        real*8, dimension(0:2)                          :: Q1, Q2, Qe, Qe1, DR, Z1, Z2, S1, S2, Pe, Pe1
        real*8, dimension(0:TPM_MAX_CHAIN-1)            :: W, C
        real*8, dimension(0:2)                          :: RR, E10
        real*8                                          :: L10, D10
        
contains !******************************************************************************************
        
        subroutine PairWeight1 ( W, E1_1, E1_2, E2_1, E2_2, R2_1, R2_2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)                     :: W
        real*8, dimension(0:2), intent(out)     :: E1_1, E1_2, E2_1, E2_2
        real*8, dimension(0:2), intent(in)      :: R2_1, R2_2
        !-------------------------------------------------------------------------------------------
        real*8                                  :: D, L20, D20, t, dWdD
        real*8, dimension(0:2)                  :: E, E20
        !-------------------------------------------------------------------------------------------
                E = 0.5d+00 * ( R2_1 + R2_2 ) - RR
                call ApplyPeriodicBC ( E )
                D = E(0) * E(0) + E(1) * E(1) + E(2) * E(2)
                if ( D < D10 * D10 ) then
                        W = 1.0d+00
                        E1_1 = 0.0d+00
                        E1_2 = 0.0d+00
                        E2_1 = 0.0d+00
                        E2_2 = 0.0d+00
                        return
                end if
                E20 = 0.5d+00 * ( R2_2 - R2_1 ) 
                L20 = sqrt ( S_V3xx ( E20 ) + sqr ( TPMR2 ) )
                D20 = L10 + L20 + TPBRcutoff + RSkin
                if ( D > D20 * D20 ) then
                        W = 0.0d+00
                        E1_1 = 0.0d+00
                        E1_2 = 0.0d+00
                        E2_1 = 0.0d+00
                        E2_2 = 0.0d+00
                        return
                end if
                D = sqrt ( D )
                E = E / D
                E20 = E20 / L20
                D20 = D20 - D10
                t = ( D - D10 ) / D20
                W = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                dWdD = 3.0d+00 * t * ( t - 1.0d+00 ) / D20
                E1_1 = dWdD * ( t * E10 - E ) 
                E1_2 = dWdD * ( - t * E10 - E ) 
                E2_1 = dWdD * ( E + t * E20 ) 
                E2_2 = dWdD * ( E - t * E20 ) 
        end subroutine PairWeight1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function EndWeight1 ( W, E1_1, E1_2, E2_1, E2_2, R1_1, R1_2, R2_1, R2_2 ) !!!!!!!!
        real*8, intent(out)                     :: W
        real*8, dimension(0:2), intent(out)     :: E1_1, E1_2, E2_1, E2_2
        real*8, dimension(0:2), intent(in)      :: R1_1, R1_2, R2_1, R2_2
        !-------------------------------------------------------------------------------------------
        real*8                                  :: D, L20
        real*8                                  :: D1, D2, t, dWdD
        real*8, dimension(0:2)                  :: RR, E, E20
        !-------------------------------------------------------------------------------------------
                E = 0.5d+00 * ( R2_1 + R2_2 - ( R1_1 + R1_2 ) )
                call ApplyPeriodicBC ( E )
                D = S_V3norm3 ( E )
                E20 = 0.5d+00 * ( R2_2 - R2_1 ) 
                L20 = sqrt ( S_V3xx ( E20 ) + sqr ( TPMR2 ) )
                D1 = L10 + L20 + TPBRcutoff + RSkin
                if ( D < D1 ) then
                        EndWeight1 = 0
                        W = 1.0d+00
                        E1_1 = 0.0d+00
                        E1_2 = 0.0d+00
                        E2_1 = 0.0d+00
                        E2_2 = 0.0d+00
                        return
                end if
                D2 = D1 + TPMC3  
                if ( D > D2 ) then
                        EndWeight1 = 2
                        W = 0.0d+00
                        E1_1 = 0.0d+00
                        E1_2 = 0.0d+00
                        E2_1 = 0.0d+00
                        E2_2 = 0.0d+00
                        return
                end if
                EndWeight1 = 1
                E = E / D
                E20 = E20 / L20
                t = ( D - D1 ) / TPMC3
                W = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                dWdD = 3.0d+00 * t * ( t - 1.0d+00 ) / TPMC3
                E1_1 = dWdD * ( E10 - E ) 
                E1_2 = dWdD * ( - E10 - E ) 
                E2_1 = dWdD * ( E + E20 ) 
                E2_2 = dWdD * ( E - E20 ) 
        end function EndWeight1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMInteractionFC1 ( Q, U, F1, F2, P1, P2, Pe, Pe1, R1, R2, Q1, Q2, Qe, Qe1, EType ) 
        real*8, intent(out)                     :: Q, U
        real*8, dimension(0:2), intent(out)     :: F1, F2, P1, P2, Pe, Pe1
        real*8, dimension(0:2), intent(in)      :: R1, R2, Q1, Q2, Qe, Qe1
        integer*4, intent(in)                   :: EType
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2)                  :: M, QX, Me, F1a, F2a, P1a, P2a, F1b, F2b, P1b, P2b, ER1, ER2, EQe, EQe1
        real*8                                  :: W, W1, D, Qa, Qb, Ua, Ub, L, Pee, Peea, Peeb, DU
        integer*4                               :: IntSigna, IntSignb, CaseID
        !-------------------------------------------------------------------------------------------
                if ( EType == 0 ) then
!                        C_TPM_0 = C_TPM_0 + 1
                        TPMInteractionFC1 = TPMInteractionF ( Q, U, F1, F2, P1, P2, Pee, R1, R2, Q1, Q2, 0 )
                        Pe = 0.0d+00
                        Pe1 = 0.0d+00
                else if ( EType < 3 ) then
!                        C_TPM_1 = C_TPM_1 + 1
                        QX = 0.5d+00 * ( Q1 + Q2 )
                        M  = Q2 - Q1
                        L  = S_V3norm3 ( M )
                        M  = M / L
                        Me = Qe - QX
                        D = S_V3norm3 ( Me )
                        if ( EType == 1 ) then
                                TPMInteractionFC1 = TPMInteractionF ( Q, U, F1, F2, P1, P2, Pee, R1, R2, QX - D * M, QX, 1 )
                        else
                                TPMInteractionFC1 = TPMInteractionF ( Q, U, F1, F2, P1, P2, Pee, R1, R2, QX, QX + D * M, 2 )
                        end if
                        call TPMSegmentForces ( P1, P2, F1, F2, R1, R2, QX, M, L )
                        Pe  = ( Pee / D ) * Me
                        Pe1 = 0.0d+00
                        QX  = 0.5d+00 * Pe
                        P1  = P1 + QX
                        P2  = P2 + QX
                else
                        CaseID = EndWeight1 ( W, ER1, ER2, EQe, Eqe1, R1, R2, Qe, Qe1 )
                        if ( CaseID < 2 ) then
                                QX = 0.5d+00 * ( Q1 + Q2 )
                                M  = Q2 - Q1
                                L  = S_V3norm3 ( M )
                                M  =  M / L
                                Me = Qe - QX
                                D = S_V3norm3 ( Me )
                                if ( EType == 3 ) then
                                        IntSigna = TPMInteractionF ( Qa, Ua, F1a, F2a, P1a, P2a, Peea, R1, R2, QX - D * M, QX, 1 )
                                else
                                        IntSigna = TPMInteractionF ( Qa, Ua, F1a, F2a, P1a, P2a, Peea, R1, R2, QX, QX + D * M, 2 )
                                end if
                                call TPMSegmentForces ( P1a, P2a, F1a, F2a, R1, R2, QX, M, L )
                        end if
                        
                        if ( CaseID > 0 ) then
                                IntSignb = TPMInteractionF ( Qb, Ub, F1b, F2b, P1b, P2b, Peeb, R1, R2, Q1, Q2, 0 )
                        end if
                        
                        if ( CaseID == 0 ) then
!                                C_TPM_1 = C_TPM_1 + 1
                                TPMInteractionFC1 = IntSigna
                                Q   = Qa
                                U   = Ua
                                F1  = F1a
                                F2  = F2a
                                Pe  = ( Peea / D ) * Me
                                Pe1 = 0.0d+00
                                QX  = 0.5d+00 * Pe
                                P1  = P1a + QX
                                P2  = P2a + QX
                        else if ( CaseID == 2 ) then
!                                C_TPM_0 = C_TPM_0 + 1
                                TPMInteractionFC1 = IntSignb
                                Q   = Qb
                                U   = Ub
                                F1  = F1b
                                F2  = F2b
                                P1  = P1b
                                P2  = P2b
                                Pe  = 0.0d+00
                                Pe1 = 0.0d+00
                        else
!                                C_TPM_2 = C_TPM_2 + 1
                                TPMInteractionFC1 = 0
                                if ( IntSigna > 0 .or. IntSignb > 0 ) TPMInteractionFC1 = 1
                                W1  = 1.0d+00 - W
                                DU  = Ub - Ua 
                                Q   = W * Qa + W1 * Qb
                                U   = W * Ua + W1 * Ub
                                Pe  = ( W * Peea / D ) * Me
                                QX  = 0.5d+00 * Pe
                                F1  = W * F1a + W1 * F1b + DU * ER1
                                F2  = W * F2a + W1 * F2b + DU * ER2
                                P1  = W * P1a + W1 * P1b + QX
                                P2  = W * P2a + W1 * P2b + QX
                                Pe  = Pe - DU * EQe
                                Pe1 = - DU * EQe1
                        end if
                end if
        end function TPMInteractionFC1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMInteractionFW1 ( QQ, U, U1, U2, UU, F1, F2, F, Fe, G1, G2, R1, R2, N, NMAX, R, Re, EType ) 
        real*8, intent(out)                             :: U, U1, U2
        integer*4, intent(in)                           :: N, NMAX, EType
        real*8, dimension(0:NMAX-1), intent(out)        :: QQ, UU
        real*8, dimension(0:2), intent(out)             :: F1, F2, Fe
        real*8, dimension(0:2,0:NMAX-1), intent(out)    :: F, G1, G2
        real*8, dimension(0:2), intent(in)              :: R1, R2, Re
        real*8, dimension(0:2,0:NMAX-1), intent(in)     :: R
        !-------------------------------------------------------------------------------------------
        integer*4                                       :: i, j
        real*8                                          :: Q, WW, DD 
        !-------------------------------------------------------------------------------------------
                Q1 = 0.0d+00
                Q2 = 0.0d+00
                WW = 0.0d+00
                Z1  = 0.0d+00
                Z2  = 0.0d+00
                TPMInteractionFW1 = 0
                E10 = 0.5d+00 * ( R2 - R1 ) 
                L10 = sqrt ( S_V3xx ( E10 ) + sqr ( TPMR1 ) )
                D10 = TPMR1 + TPMR2 + TPMC123 * TPBRcutoff + RSkin
                E10 = E10 / L10
                RR  = 0.5d+00 * ( R1 + R2 )
                do i = 0, N - 2
                        call PairWeight1 ( W(i), E1(0:2,i), E2(0:2,i), EE1(0:2,i), EE2(0:2,i), R(0:2,i), R(0:2,i+1) )
                        Q1 = Q1 + W(i) * R(0:2,i)
                        Q2 = Q2 + W(i) * R(0:2,i+1)
                        WW = WW + W(i)
                        Z1  = Z1 + E1(0:2,i)
                        Z2  = Z2 + E2(0:2,i)
                end do
                if ( WW .le. TPGeomPrec ) return
                Q1 = Q1 / WW
                Q2 = Q2 / WW
                Z1  = Z1 / WW
                Z2  = Z2 / WW
                if ( EType == 1 ) then
                        Qe = R(0:2,0)
                        Qe1 = R(0:2,1)
                else if ( EType == 2 ) then
                        Qe = R(0:2,N-1)
                        Qe1 = R(0:2,N-2)
                else if ( EType == 3 ) then
                        Qe = Re
                        Qe1 = R(0:2,0)
                else if ( EType == 4 ) then
                        Qe = Re
                        Qe1 = R(0:2,N-1)
                else
                        Qe = 0.0d+00
                        Qe1 = 0.0d+00
                end if
                
                TPMInteractionFW1 = TPMInteractionFC1 ( Q, U, F1, F2, S1, S2, Pe, Pe1, R1, R2, Q1, Q2, Qe, Qe1, EType ) 
                if ( TPMInteractionFW1 == 0 ) return
                
                W(0:N-2)  = W(0:N-2) / WW
                E1(0:2,0:N-2) = E1(0:2,0:N-2) / WW
                E2(0:2,0:N-2) = E2(0:2,0:N-2) / WW
                EE1(0:2,0:N-2) = EE1(0:2,0:N-2) / WW
                EE2(0:2,0:N-2) = EE2(0:2,0:N-2) / WW
                G1(0:2,0:N-1) = 0.0d+00
                G2(0:2,0:N-1) = 0.0d+00
                U1 = 0.25d+00 * U
                U2 = U1
                UU = 0.0d+00
                do i = 0, N - 2
                        QQ(i)   = W(i) * Q
                        DD      = W(i) * U1
                        UU(i)   = UU(i) + DD
                        UU(i+1) = UU(i+1) + DD
                end do
                do i = 0, N - 2
                        C(i) = S_V3xV3 ( S1, R(0:2,i) ) + S_V3xV3 ( S2, R(0:2,i+1) )
                        F1 = F1 + C(i) * ( E1(0:2,i) - W(i) * Z1 )
                        F2 = F2 + C(i) * ( E2(0:2,i) - W(i) * Z2 )
                end do
                F(0:2,0) = W(0) * S1
                do j = 0, N - 2
                        if ( j == 0 ) then
                                DR = EE1(0:2,0) * ( 1.0d+00 - W(0) )
                        else 
                                DR = - W(j) * EE1(0:2,0)
                        end if
                        F(0:2,0) = F(0:2,0) + C(j) * DR
                end do
                do i = 1, N - 2
                        G1(0:2,i) = W(i-1) * S2
                        G2(0:2,i) = W(i) * S1
                        do j = 0, N - 2
                                if ( j == i ) then
                                        G1(0:2,i) = G1(0:2,i) - C(j) * W(j) * EE2(0:2,i-1)
                                        G2(0:2,i) = G2(0:2,i) + C(j) * ( EE1(0:2,j) - W(j) * EE1(0:2,i) )
                                else if ( j == i - 1 ) then
                                        G1(0:2,i) = G1(0:2,i) + C(j) * ( EE2(0:2,j) - W(j) * EE2(0:2,i-1) )
                                        G2(0:2,i) = G2(0:2,i) - C(j) * W(j) * EE1(0:2,i)
                                else 
                                        G1(0:2,i) = G1(0:2,i) - C(j) * W(j) * EE2(0:2,i-1) 
                                        G2(0:2,i) = G2(0:2,i) - C(j) * W(j) * EE1(0:2,i)
                                end if
                        end do
                        F(0:2,i) = G1(0:2,i) + G2(0:2,i)
                end do
                F(0:2,N-1) = W(N-2) * S2
                do j = 0, N - 2
                        if ( j == N - 2 ) then
                                DR = EE2(0:2,N-2) * ( 1.0d+00 - W(N-2) )
                        else 
                                DR = - W(j) * EE2(0:2,N-2)
                        end if
                        F(0:2,N-1) = F(0:2,N-1) + C(j) * DR
                end do
                Fe = 0.0d+00
                if ( EType == 1 ) then
                        F(0:2,0) = F(0:2,0) - Pe
                else if ( EType == 2 ) then
                        F(0:2,N-1) = F(0:2,N-1) - Pe
                else if ( EType == 3 ) then
                        F(0:2,0) = F(0:2,0) - Pe1
                        Fe = - Pe
                else if ( EType == 4 ) then
                        F(0:2,N-1) = F(0:2,N-1) - Pe1
                        Fe = - Pe
                end if
                G1(0:2,N-1) = F(0:2,N-1)
                G2(0:2,0) = F(0:2,0)
        end function TPMInteractionFW1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module TPMM1 !**********************************************************************************

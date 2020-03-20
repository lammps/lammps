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

module TPMForceField !************************************************************************************
!
! TMD Library: Calculation of the TMD force field
!
!---------------------------------------------------------------------------------------------------
!
! PGI Fortran, Intel Fortran
!
! Alexey N. Volkov, University of Alabama (avolkov1@ua.edu), Version 09.01.33, 2018
!
!***************************************************************************************************

use CNTPot
use TPMM0
use TPMM1
use iso_c_binding, only : c_int, c_double, c_char
implicit none

contains !******************************************************************************************
        
        subroutine TubeStretchingForceField ( U1, U2, F1, F2, S1, S2, X1, X2, R12, L12 ) !!!!!!!!!!!
        real(c_double), intent(inout)                           :: U1, U2       ! Interaction energies associated with nodes X1 and X2
        real(c_double), intent(inout), dimension(0:2)           :: F1, F2       ! Forces exerted on nodes X1 and X2
        real(c_double), intent(inout), dimension(0:2,0:2)       :: S1, S2       ! Contributions of nodes X1 and X2 to the virial stress tensor
        real(c_double), intent(in), dimension(0:2)              :: X1, X2       ! Coordinates of the segmnet nodes 
        real(c_double), intent(in)                              :: R12          ! Radius of nanotube the segment (X1,X2) belongs to
        real(c_double), intent(in)                              :: L12          ! Equilubrium length of segment (X1,X2)
        !-------------------------------------------------------------------------------------------
        integer(c_int)                                       :: ii, jj, Event
        real(c_double)                                          :: U, F, LL, S, Ubcl
        real(c_double), dimension(0:2)                          :: DX, FF
        !-------------------------------------------------------------------------------------------
                DX = X2 - X1 
                LL = S_V3norm3 ( DX )
                Event = CNTSTRCalc ( U, F, LL, R12, L12, 0, Ubcl )

                U = U / 2.0d+00
                FF = DX * F / LL
                
                F1 = F1 + FF
                U1 = U1 + U

                F2 = F2 - FF
                U2 = U2 + U

                ! Stress
                do ii = 0, 2
                        do jj = 0, 2
                                S = - 0.5d+00 * DX(ii) * FF(jj)
                                S1(ii,jj) = S1(ii,jj) + S 
                                S2(ii,jj) = S2(ii,jj) + S 
                        end do
                end do
        end subroutine TubeStretchingForceField !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TubeBendingForceField ( U1, U2, U3, F1, F2, F3, S1, S2, S3, X1, X2, X3, R123, L123, BBF2 ) 
        real(c_double), intent(inout)                           :: U1, U2, U3   ! Interaction energies associated with nodes X1, X2, and X3
        real(c_double), intent(inout), dimension(0:2)           :: F1, F2, F3   ! Forces exerted on nodes X1, X2, and X3
        real(c_double), intent(inout), dimension(0:2,0:2)       :: S1, S2, S3   ! Contributions of nodes X1, X2, and X3 to the virial stress tensor
        real(c_double), intent(in), dimension(0:2)              :: X1, X2, X3   ! Coordinates of nodes 
        real(c_double), intent(in)                              :: R123         ! Radius of nanotube the segment (X1,X2) belongs to
        real(c_double), intent(in)                              :: L123         ! Equilubrium length of segment (X1,X2) and (X2,X3) (It is assumed to be the same for both segments)
        integer(c_int), intent(inout)                        :: BBF2
        !-------------------------------------------------------------------------------------------
        integer(c_int)                                       :: ii, jj, Event
        real(c_double)                                          :: U, F, K, S, Ubcl
        real(c_double), dimension(0:2)                          :: G0, G1, G2
        !-------------------------------------------------------------------------------------------
                call BendingGradients ( K, G0, G1, G2, X1, X2, X3 )
                Event = CNTBNDCalc ( U, F, K, R123, L123, BBF2, Ubcl )

                if ( Event == CNTPOT_BBUCKLING ) then
                        BBF2 = 1
                else
                        BBF2 = 0
                end if

                U = U / 4.0d+00
                F = - F

                F1 = F1 + G0 * F 
                F2 = F2 + G1 * F 
                F3 = F3 + G2 * F 

                U1 = U1 + U
                U2 = U2 + 2.0d+00 * U
                U3 = U3 + U

                ! Stress
                do ii = 0, 2
                        do jj = 0, 2
                                S = 0.5d+00 * ( X1(ii) - X2(ii) ) * G0(jj)
                                S1(ii,jj) = S1(ii,jj) + S 
                                S2(ii,jj) = S2(ii,jj) + S 
                                S = 0.5d+00 * ( X3(ii) - X2(ii) ) * G2(jj)
                                S3(ii,jj) = S3(ii,jj) + S 
                                S2(ii,jj) = S2(ii,jj) + S 
                        end do
                end do
        end subroutine TubeBendingForceField !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! The purpose of subroutine SegmentTubeForceField is to calculate interaction forces 
        ! (as well potential nergies and componets of the virial stress tensor) between a segment
        ! (X1,X2) and a sequence of segments with node coordinates that belongs to a single CNT

        ! It is assumed that X contains ALL nodes of a single CNT that are included into the
        ! neighbor list of segment (X1,X2)

        ! The nodes in X are assumed to be ordered according to their physical appearence in the nanotube
        ! It means that (X(i),X(i+1)) are either correspond to a real segment or divided by a segments 
        ! that do not belong to a nanotube.

        ! Concept of the extendend chain:
        ! Let's consider a sequant of nodes (X1,X2,...,XN) forming continuous part of a nanotube.
        ! If node Xe preceeds X1 and Xe is the nanotube end, then the extended chain is (Xe,X1,...,XN) and Ee = 1.
        ! If node Xe follows XN and Xe is the nanotube end, then the extended chain is (X1,...,XN,Xe) and Ee = 2.
        ! In all other cases, extended chain coincides with (X1,...,XN) and Ee = 0
        ! If the extended chain contains additional node, then non-zero force is exterted on this node
        
        subroutine SegmentTubeForceField ( U1, U2, U, F1, F2, F, Fe, S1, S2, S, Se, X1, X2, R12, N, X, Xe, BBF, R, E1, E2, Ee, TPMType ) 
        integer(c_int), intent(in)                           :: N            ! Number of nodes in array X
        real(c_double), intent(inout)                           :: U1, U2       ! Interaction energies associated with nodes X1 and X2
        real(c_double), intent(inout), dimension(0:N-1)         :: U            ! Interaction energies associated with nodes X
        real(c_double), intent(inout), dimension(0:2)           :: F1, F2       ! Forces exerted on nodes X1 and X2
        real(c_double), intent(inout), dimension(0:2,0:N-1)     :: F            ! Forces exerted on nodes X
        real(c_double), intent(inout), dimension(0:2)           :: Fe           ! Force exerted on node Xe (can be updated only if Ee > 0)
        real(c_double), intent(inout), dimension(0:2,0:2)       :: S1, S2       ! Contributions of nodes X1 and X2 to the virial stress tensor
        real(c_double), intent(inout), dimension(0:2,0:2,0:N-1) :: S            ! Contributions of nodes X to the virial stress tensor
        real(c_double), intent(inout), dimension(0:2,0:2)       :: Se           ! Contributions of node Xe to the virial stress tensor (can be updated only if Ee > 0)
        real(c_double), intent(in), dimension(0:2)              :: X1, X2       ! Coordinates of the segmnet nodes 
        real(c_double), intent(in)                              :: R12          ! Radius of nanotube the segment (X1,X2) belongs to
        real(c_double), intent(in), dimension(0:2,0:N-1)        :: X            ! Coordinates of the nanotube nodes
        real(c_double), intent(in), dimension(0:2)              :: Xe           ! Additiona node of the extended chain if Ee > 0
        integer(c_int), intent(in), dimension(0:N-1)         :: BBF          ! Bending buckling flags (BBF(i) = 1 in a case of buckling in node i)
        real(c_double), intent(in)                              :: R            ! Radius of nanotube X
        integer(c_int), intent(in)                           :: E1, E2       ! E1 = 1 if the chnane node 0 is a CNT end; E1 = 2 if the chnane node N-1 is a CNT end;
        integer(c_int), intent(in)                           :: Ee           ! Parameter defining the type of the extended chain (0,1,2)
        integer(c_int), intent(in)                           :: TPMType      ! Type of the tubular potential (0 or 1)
        !-------------------------------------------------------------------------------------------
        integer(c_int)                                       :: k, ii, jj, IntSign
        integer(c_int)                                       :: BType, EType, LocalTPMType
        real(c_double), dimension(0:2,0:N-1)                    :: G1, G2
        real(c_double), dimension(0:N-1)                        :: QQ
        logical                                         :: EType1, EType2
        real(c_double), dimension(0:2)                          :: G, DG, DQ, XX
        real(c_double)                                          :: UT, DR, DS, DS1
        real(c_double)                                          :: xU1, xU2     ! Interaction energies associated with nodes X1 and X2
        real(c_double), dimension(0:N-1)                        :: xU           ! Interaction energies associated with nodes X
        real(c_double), dimension(0:2)                          :: xF1, xF2     ! Forces exerted on nodes X1 and X2
        real(c_double), dimension(0:2,0:N-1)                    :: xF           ! Forces exerted on nodes X
        real(c_double), dimension(0:2)                          :: xFe          ! Force exerted on node Xe (can be updated only if Ee > 0)
        !-------------------------------------------------------------------------------------------

                !U1 = 0.0d+00
                !U2 = 0.0d+00 
                !U = 0.0d+00 
                !F1 = 0.0d+00 
                !F2 = 0.0d+00 
                !F = 0.0d+00 
                !S1 = 0.0d+00 
                !S2 = 0.0d+00 
                !S = 0.0d+00

                ! Looking for a buckling point 
                BType = 0
                do k = 0, N - 1
                        if ( BBF(k) == 1 ) then
                                BType = 1
                                exit
                        end if
                end do

                ! Choosing the LocalTPMType and Etype.
                ! LocalTPMType is set to 0 if both ends of the chain are nanotube ends or the chain contains a buckling point. 
                ! Overwise, LocalTPMType = TPMType.
                if ( BType == 1 ) then
                        LocalTPMType = 0
                        EType = 0
                else
                        if ( E1 == 1 ) then ! First node in the chain is the tube end
                                EType1 = .true.
                        else
                                EType1 = .false.
                        end if
                        if ( E2 == 1 ) then ! Last node in the chain is the tube end
                                EType2 = .true.
                        else
                                EType2 = .false.
                        end if
                        if ( EType1 .and. EType2 ) then
                                LocalTPMType = 0
                        else
                                LocalTPMType = TPMType 
                                if ( EType1 ) then
                                        EType = 1
                                else if ( EType2 ) then
                                        EType = 2
                                else ! No tube ends in the chain
                                        EType = 0
                                end if
                        end if
                end if
                                        
                if ( LocalTPMType == 0 ) then
                        IntSign = TPMInteractionFW0 ( QQ, UT, xU1, xU2, xU, xF1, xF2, xF, G1, G2, X1, X2, N, N, X )
                else
                        if ( EType == 0 ) then
                                if ( Ee == 1 ) then ! First node in the extended chain is the tube end
                                        EType = 3
                                else if ( Ee == 2 ) then ! Last node in the extended chain is the tube end
                                        EType = 4
                                end if
                        end if
                        IntSign = TPMInteractionFW1 ( QQ, UT, xU1, xU2, xU, xF1, xF2, xF, xFe, G1, G2, X1, X2, N, N, X, Xe, EType )
                end if
                if ( IntSign == 0 ) return ! No interaction

                ! Final potential energies
                U1 = U1 + 0.5d+00 * xU1
                U2 = U2 + 0.5d+00 * xU2
                U(0:N-1) = U(0:N-1) + 0.5d+00 * xU(0:N-1)

                ! Contributions to the virial stresses tensor
                do ii = 0, 2
                        DR = 0.125d+00 * ( X2(ii) - X1(ii) ) 
                        do jj = 0, 2
                                DS = DR * ( xF2(jj) - xF1(jj) )
                                S1(ii,jj) = S1(ii,jj) + DS 
                                S2(ii,jj) = S2(ii,jj) + DS
                        end do
                end do
                XX = 0.5d+00 * ( X2 + X1 )
                if ( EType > 2 ) then
                        DQ =  Xe - XX
                        call ApplyPeriodicBC ( DQ )
                        DQ = DQ / 6.0d+00
                        do ii = 0, 2
                                do jj = 0, 2
                                        DS = DQ(ii) * xFe(jj)
                                        S1(ii,jj) = S1(ii,jj) + DS 
                                        S2(ii,jj) = S1(ii,jj) + DS
                                        Se(ii,jj) = Se(ii,jj) + DS
                                end do
                        end do
                end if
                do k = 0, N - 2
                        DQ = 0.5d+00 * ( X(0:2,k+1) + X(0:2,k) ) - XX 
                        call ApplyPeriodicBC ( DQ )
                        DQ = 0.125d+00 * DQ 
                        G  = G1(0:2,k+1) + G2(0:2,k) 
                        DG = G1(0:2,k+1) - G2(0:2,k)
                        do ii = 0, 2
                                DR = 0.125d+00 * ( X(ii,k+1) - X(ii,k) ) 
                                do jj = 0, 2
                                        DS  = DQ(ii) * G(jj)
                                        DS1 = DS + DR * DG(jj)
                                        S1(ii,jj) = S1(ii,jj) + DS 
                                        S2(ii,jj) = S2(ii,jj) + DS
                                        S(ii,jj,k) = S(ii,jj,k) + DS1
                                        S(ii,jj,k+1) = S(ii,jj,k+1) + DS1
                                end do
                        end do
                end do

                ! Final forces
                F1 = F1 + 0.5d+00 * xF1
                F2 = F2 + 0.5d+00 * xF2
                F(0:2,0:N-1) = F(0:2,0:N-1) + 0.5d+00 * xF(0:2,0:N-1)
                if ( EType > 2 ) then
                        Fe = Fe + 0.5d+00 * xFe
                end if

        end subroutine SegmentTubeForceField !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TPMForceField !**************************************************************************

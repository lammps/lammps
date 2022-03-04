! ------------ ----------------------------------------------------------
!   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
!   https://www.lammps.org/ Sandia National Laboratories
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

module TubePotTrue !********************************************************************************
!
! TMD Library: True tubular potential and transfer function
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!---------------------------------------------------------------------------------------------------
!
! This module implements calculation of the true potential and transfer functions for interaction
! between two cylinder segments of nanotubes by direct integration over the surfaces of both
! segments.
!
!***************************************************************************************************

use TPMGeom
use TubePotBase
use iso_c_binding, only : c_int, c_double, c_char
implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------

        integer(c_int), parameter    :: TPTNXMAX     = 257
        integer(c_int), parameter    :: TPTNEMAX     = 128

!---------------------------------------------------------------------------------------------------
! Types
!---------------------------------------------------------------------------------------------------

        type TPTSEG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double)                  :: X, Y, Z
        real(c_double)                  :: Psi, Theta, Phi              ! Euler's angles
        real(c_double)                  :: R                            ! Segment radius
        real(c_double)                  :: L                            ! Segment length
        integer(c_int)               :: NX, NE                          ! Number of nodes for numerical integration
        real(c_double)                  :: DX, DE                       ! Spacings
        real(c_double), dimension(0:2,0:2) :: M                         ! Transformation matrix
        real(c_double), dimension(0:TPTNXMAX-1,0:TPTNXMAX-1,0:2) :: Rtab! Node coordinates
        end type TPTSEG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------

        type(TPTSEG)    :: TPTSeg1, TPTSeg2     ! Two segments

contains !******************************************************************************************

        subroutine TPTSegAxisVector ( S, Laxis ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(in)                :: S
        real(c_double), dimension(0:2), intent(out)     :: Laxis
        !-------------------------------------------------------------------------------------------
                Laxis(0:2) = S%M(2,0:2)
        end subroutine TPTSegAxisVector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPTSegRadVector ( S, Lrad, Eps ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(in)                :: S
        real(c_double), dimension(0:2), intent(out)     :: Lrad
        real(c_double), intent(in)                      :: Eps
        !-------------------------------------------------------------------------------------------
        real(c_double)                                  :: Ce, Se
        !-------------------------------------------------------------------------------------------
                Ce = cos ( Eps )
                Se = sin ( Eps )
                Lrad(0) = Ce * S%M(0,0) + Se * S%M(1,0)
                Lrad(1) = Ce * S%M(0,1) + Se * S%M(1,1)
                Lrad(2) = Ce * S%M(0,2) + Se * S%M(1,2)
        end subroutine TPTSegRadVector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPTRadiusVector ( S, R, X, Eps ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(in)                :: S
        real(c_double), dimension(0:2), intent(out)     :: R
        real(c_double), intent(in)                      :: X, Eps
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2)                  :: Laxis, Lrad
        !-------------------------------------------------------------------------------------------
                call TPTSegAxisVector ( S, Laxis )
                call TPTSegRadVector ( S, Lrad, Eps )
                R(0) = S%X + X * Laxis(0) + S%R * Lrad(0)
                R(1) = S%Y + X * Laxis(1) + S%R * Lrad(1)
                R(2) = S%Z + X * Laxis(2) + S%R * Lrad(2)
        end subroutine TPTRadiusVector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPTCalcSegNodeTable ( S ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(inout)     :: S
        !-------------------------------------------------------------------------------------------
        real(c_double)                          :: X, Eps
        integer(c_int)                          :: i, j
        !-------------------------------------------------------------------------------------------
                X = - S%L / 2.0
                call RotationMatrix3 ( S%M, S%Psi, S%Theta, S%Phi )
                do i = 0, S%NX - 1
                        Eps = 0.0d+00
                        do j = 0, S%NE - 1
                                call TPTRadiusVector ( S, S%Rtab(i,j,0:2), X, Eps )
                                Eps = Eps + S%DE
                        end do
                        X = X + S%DX
                end do
        end subroutine TPTCalcSegNodeTable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPTSetSegPosition1 ( S, Rcenter, Laxis, L ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(inout)                     :: S
        real(c_double), dimension(0:2), intent(in)      :: Rcenter, Laxis
        real(c_double), intent(in)                      :: L
        !-------------------------------------------------------------------------------------------
                S%L  = L
                S%DX = L / ( S%NX - 1 )
                call EulerAngles ( S%Psi, S%Theta, Laxis )
                S%Phi= 0.0d+00
                S%X  = Rcenter(0)
                S%Y  = Rcenter(1)
                S%Z  = Rcenter(2)
                call TPTCalcSegNodeTable ( S )
        end subroutine TPTSetSegPosition1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPTSetSegPosition2 ( S, R1, R2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(inout)                     :: S
        real(c_double), dimension(0:2), intent(in)      :: R1, R2
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2)                  :: R, Laxis
        real(c_double)                                  :: L
        !-------------------------------------------------------------------------------------------
                R = 0.5 * ( R1 + R2 )
                Laxis = R2 - R1
                L = S_V3norm3 ( Laxis )
                Laxis = Laxis / L
                call TPTSetSegPosition1 ( S, R, Laxis, L )
        end subroutine TPTSetSegPosition2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function TPTCheckIntersection ( S1, S2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(in)        :: S1, S2
        !-------------------------------------------------------------------------------------------
        integer(c_int)                          :: i, j
        real(c_double)                          :: L1, L2, Displacement, D
        real(c_double), dimension(0:2)          :: Laxis, Q, R
        !-------------------------------------------------------------------------------------------
                L2 = S1%L / 2.0
                L1 = - L2
                call TPTSegAxisVector ( S1, Laxis )
                R(0) = S1%X
                R(1) = S1%Y
                R(2) = S1%Z
                do i = 0, S2%NX - 1
                        do j = 0, S2%NE - 1
                                call LinePoint ( Displacement, Q, R, Laxis, S2%Rtab(i,j,0:2) )
                                D = sqrt ( sqr ( Q(0) - S2%Rtab(i,j,0) ) + sqr ( Q(1) - S2%Rtab(i,j,1) ) &
                                        + sqr ( Q(2) - S2%Rtab(i,j,2) )  )
                                if ( Displacement > L1 .and. Displacement < L2 .and. D < S1%R ) then
                                        TPTCheckIntersection = 1
                                        return
                                end if
                        end do
                end do
                TPTCheckIntersection = 0
        end function TPTCheckIntersection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function TPTCalcPointRange ( S, Xmin, Xmax, Re ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(TPTSEG), intent(in)                        :: S
        real(c_double), intent(out)                     :: Xmin, Xmax
        real(c_double), dimension(0:2), intent(in)      :: Re
        !-------------------------------------------------------------------------------------------
        real(c_double)                                  :: Displacement, Distance
        real(c_double), dimension(0:2)                  :: Laxis, Q, R
        !-------------------------------------------------------------------------------------------
                call TPTSegAxisVector ( S, Laxis )
                R(0) = S%X
                R(1) = S%Y
                R(2) = S%Z
                call LinePoint ( Displacement, Q, R, Laxis, Re )
                Distance = sqrt ( sqr ( Q(0) - Re(0) ) + sqr ( Q(1) - Re(1) ) + sqr ( Q(2) - Re(2) ) ) - S%R
                if ( TPBRcutoff < Distance ) then
                        Xmin = 0.0d+00
                        Xmax = 0.0d+00
                        TPTCalcPointRange = 0
                        return
                end if
                Distance = sqrt ( TPBRcutoff * TPBRcutoff - Distance * Distance )
                Xmin = Displacement - Distance
                Xmax = Displacement + Distance
                TPTCalcPointRange = 1
        end function TPTCalcPointRange !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPTGetEnds ( R1_1, R1_2, R2_1, R2_2, X1_1, X1_2, X2_1, X2_2, H, A ) !!!!!!!!!!!!!
        real(c_double), dimension(0:2), intent(out)     :: R1_1, R1_2, R2_1, R2_2
        real(c_double), intent(in)                      :: X1_1, X1_2, X2_1, X2_2, H, A
        !-------------------------------------------------------------------------------------------
                R1_1(0) = 0.0d+00
                R1_1(1) = 0.0d+00
                R1_1(2) = X1_1
                R1_2(0) = 0.0d+00
                R1_2(1) = 0.0d+00
                R1_2(2) = X1_2
                R2_1(0) = H
                R2_1(1) = - X2_1 * sin ( A )
                R2_1(2) = X2_1 * cos ( A )
                R2_2(0) = H
                R2_2(1) = - X2_2 * sin ( A )
                R2_2(2) = X2_2 * cos ( A )
        end subroutine TPTGetEnds !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Tubular potential
!---------------------------------------------------------------------------------------------------

        integer(c_int) function TPTPointPotential ( Q, U, F, R, S ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the potential U and force F applied to an atom in position R and
        ! produced by the segment S.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)                     :: Q, U
        real(c_double), dimension(0:2), intent(out)     :: F
        real(c_double), dimension(0:2), intent(in)      :: R
        type(TPTSEG), intent(in)                        :: S
        !-------------------------------------------------------------------------------------------
        integer(c_int)                                  :: i, j
        real(c_double), dimension(0:2)                  :: RR, FF
        real(c_double)                                  :: QQ, UU, UUU, FFF, Rabs
        real(c_double)                                  :: Coeff, Xmin, Xmax, X
        !-------------------------------------------------------------------------------------------
                TPTPointPotential = 0
                Q = 0.0d+00
                U = 0.0d+00
                F = 0.0d+00
                if ( TPTCalcPointRange ( S, Xmin, Xmax, R ) == 0 ) return
                X = - S%L / 2.0
                do i = 0, S%NX - 1
                        if ( X > Xmin .and. X < Xmax ) then
                                QQ = 0.0d+00
                                UU = 0.0d+00
                                FF = 0.0d+00
                                do j = 0, S%NE - 1
                                        RR(0:2) = S%Rtab(i,j,0:2) - R(0:2)
                                        Rabs = S_V3norm3 ( RR )
                                        if ( Rabs < TPBRcutoff ) then
                                                QQ  = QQ + TPBQCalc0 ( Rabs )
                                                call TPBUCalc1 ( UUU, FFF, Rabs )
                                                UU  = UU + UUU
                                                FFF = FFF / Rabs
                                                FF  = FF + FFF * RR
                                                TPTPointPotential = 1
                                        end if
                                end do
                                if ( i == 0 .or. i == S%NX - 1 ) then
                                        Q = Q + 0.5d+00 * QQ
                                        U = U + 0.5d+00 * UU
                                        F = F + 0.5d+00 * FF
                                else
                                        Q = Q + QQ
                                        U = U + UU
                                        F = F + FF
                                end if
                        end if
                        X = X + S%DX
                end do
                Coeff = TPBD * S%DX * S%R * S%DE
                Q = Q * S%DX * S%R * S%DE
                U = U * Coeff
                F = F * Coeff
        end function TPTPointPotential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function TPTSectionPotential ( Q, U, F, M, S, i, Ssource ) !!!!!!!!!!!!!!!!!!
        ! This function returns the potential U, force F and torque M produced by the segment Ssource
        ! and applied to the i-th circular cross-section of the segment S.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)                     :: Q, U
        real(c_double), dimension(0:2), intent(out)     :: F, M
        type(TPTSEG), intent(in)                        :: S, Ssource
        integer(c_int), intent(in)                      :: i
        !-------------------------------------------------------------------------------------------
        integer(c_int)                                  :: j
        real(c_double), dimension(0:2)                  :: R, Fp, Mp, Lrad
        real(c_double)                                  :: Qp, Up, Eps
        real(c_double)                                  :: Coeff
        !-------------------------------------------------------------------------------------------
                TPTSectionPotential = 0
                Q = 0.0d+00
                U = 0.0d+00
                F = 0.0d+00
                M = 0.0d+00
                Eps = 0.0d+00
                do j = 0, S%NE - 1
                        call TPTSegRadVector ( S, Lrad, Eps )
                        if ( TPTPointPotential ( Qp, Up, Fp, S%Rtab(i,j,0:2), Ssource ) == 1 ) then
                                Q    = Q + Qp
                                U    = U + Up
                                F    = F + Fp
                                R(0) = S%Rtab(i,j,0) - S%X
                                R(1) = S%Rtab(i,j,1) - S%Y
                                R(2) = S%Rtab(i,j,2) - S%Z
                                call V3_V3xxV3 ( Mp, R, Fp )
                                M = M + Mp
                                TPTSectionPotential = 1
                        end if
                        Eps = Eps + S%DE
                end do
                Coeff = TPBD * S%R * S%DE
                Q     = Q * S%R * S%DE
                U     = U * Coeff
                F     = F * Coeff
                M     = M * Coeff
        end function TPTSectionPotential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function TPTSegmentPotential ( Q, U, F, M, S, Ssource ) !!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the potential U, force F and torque M produced by the segment
        ! Ssource and applied to the segment S.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)                     :: Q, U
        real(c_double), dimension(0:2), intent(out)     :: F, M
        type(TPTSEG), intent(in)                :: S, Ssource
        integer(c_int)                               :: i
        real(c_double), dimension(0:2)                  :: Fc, Mc
        real(c_double)                                  :: Qc, Uc
        !-------------------------------------------------------------------------------------------
                TPTSegmentPotential = 0
                Q = 0.0d+00
                U = 0.0d+00
                F = 0.0d+00
                M = 0.0d+00
                if ( TPTCheckIntersection ( S, Ssource ) == 1 ) then
                         TPTSegmentPotential = 2
                        return
                end if
                do i = 0, S%NX - 1
                        if ( TPTSectionPotential ( Qc, Uc, Fc, Mc, S, i, Ssource ) == 1 ) then
                                if ( i == 0 .or. i == S%NX - 1 ) then
                                        Q = Q + 0.5d+00 * Qc
                                        U = U + 0.5d+00 * Uc
                                        F = F + 0.5d+00 * Fc
                                        M = M + 0.5d+00 * Mc
                                else
                                        Q = Q + Qc
                                        U = U + Uc
                                        F = F + Fc
                                        M = M + Mc
                                end if
                                TPTSegmentPotential = 1
                        end if
                end do
                Q = Q * S%DX
                U = U * S%DX
                F = F * S%DX
                M = M * S%DX
        end function TPTSegmentPotential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Forces
!---------------------------------------------------------------------------------------------------

        subroutine TPTSegmentForces ( F1, F2, F, M, Laxis, L ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), dimension(0:2), intent(out)     :: F1, F2
        real(c_double), dimension(0:2), intent(in)      :: F, M, Laxis
        real(c_double), intent(in)                      :: L
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2)                  :: MM, FF, FFF
        !-------------------------------------------------------------------------------------------
                FF = 0.5d+00 * F
                MM = M / L
                call V3_V3xxV3 ( FFF, MM, Laxis )
                F1 = FF - FFF
                F2 = FF + FFF
        end subroutine TPTSegmentForces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function TPTInteractionF ( Q, U, F1_1, F1_2, F2_1, F2_2, R1_1, R1_2, R2_1, R2_2 )
        ! This function returns the potential and forces applied to the ends of segments.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)                     :: Q, U
        real(c_double), dimension(0:2), intent(out)     :: F1_1, F1_2, F2_1, F2_2
        real(c_double), dimension(0:2), intent(in)      :: R1_1, R1_2, R2_1, R2_2
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2)                  :: R1, R2, Laxis1, Laxis2, DR, F1, M1, F2, M2
        real(c_double)                                  :: L1, L2
        !-------------------------------------------------------------------------------------------
                R1 = 0.5 * ( R1_1 + R1_2 )
                R2 = 0.5 * ( R2_1 + R2_2 )
                Laxis1 = R1_2 - R1_1
                Laxis2 = R2_2 - R2_1
                L1 = S_V3norm3 ( Laxis1 )
                L2 = S_V3norm3 ( Laxis2 )
                Laxis1 = Laxis1 / L1
                Laxis2 = Laxis2 / L2
                DR = R2 - R1
                call TPTSetSegPosition1 ( TPTSeg1, R1, Laxis1, L1 )
                call TPTSetSegPosition1 ( TPTSeg2, R2, Laxis2, L2 )
                TPTInteractionF = TPTSegmentPotential ( Q, U, F1, M1, TPTSeg1, TPTSeg2 )
                if ( TPTInteractionF .ne. 1 ) return
                call V3_V3xxV3 ( M2, DR, F1 )
                F2 = - F1
                M2 = - M1 - M2
                call TPTSegmentForces ( F1_1, F1_2, F1, M1, Laxis1, L1 )
                call TPTSegmentForces ( F2_1, F2_2, F2, M2, Laxis2, L2 )
        end function TPTInteractionF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Initialization
!---------------------------------------------------------------------------------------------------

        subroutine TPTInit ( R1, R2, NX, NE ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: R1, R2
        integer(c_int), intent(in)   :: NX, NE
        !-------------------------------------------------------------------------------------------
                TPTSeg1%X     = 0.0d+00
                TPTSeg1%Y     = 0.0d+00
                TPTSeg1%Z     = 0.0d+00
                TPTSeg1%Psi   = 0.0d+00
                TPTSeg1%Theta = 0.0d+00
                TPTSeg1%Phi   = 0.0d+00
                TPTSeg1%R     = R1
                TPTSeg1%NX    = NX
                TPTSeg1%NE    = NE
                TPTSeg1%DE    = M_2PI / NE
                TPTSeg2%X     = 0.0d+00
                TPTSeg2%Y     = 0.0d+00
                TPTSeg2%Z     = 0.0d+00
                TPTSeg2%Psi   = 0.0d+00
                TPTSeg2%Theta = 0.0d+00
                TPTSeg2%Phi   = 0.0d+00
                TPTSeg2%R     = R2
                TPTSeg2%NX    = NX
                TPTSeg2%NE    = NE
                TPTSeg2%DE    = M_2PI / NE
        end subroutine TPTInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TubePotTrue !****************************************************************************

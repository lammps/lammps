module TPMGeom !************************************************************************************
!
! Geometry functions for TPM force field
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
!
!***************************************************************************************************

use TPMLib

implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------

        integer*4, parameter    :: MD_LINES_NONPAR      = 0
        integer*4, parameter    :: MD_LINES_PAR         = 1

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------
        
        ! Coordinates of the whole domain
        real*8                  :: DomXmin, DomXmax, DomYmin, DomYmax, DomZmin, DomZmax
        real*8                  :: DomLX, DomLY, DomLZ
        real*8                  :: DomLXhalf, DomLYhalf, DomLZhalf
        
        ! Boundary conditions 
        integer*4               :: BC_X                 = 0
        integer*4               :: BC_Y                 = 0
        integer*4               :: BC_Z                 = 0

        ! Skin parameter in NBL and related algorithms
        real*8                  :: Rskin                = 1.0d+00
        
contains !******************************************************************************************

        subroutine ApplyPeriodicBC ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This subroutine changes coortinates of the point accorning to periodic boundary conditions
        ! it order to makesure that the point is inside the computational cell
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2), intent(inout)   :: R
        !-------------------------------------------------------------------------------------------
                ! These commented lines implemment the more general, but less efficient algorithm
                !if ( BC_X == 1 ) R(0) = R(0) - DomLX * roundint ( R(0) / DomLX )
                !if ( BC_Y == 1 ) R(1) = R(1) - DomLY * roundint ( R(1) / DomLY )
                !if ( BC_Z == 1 ) R(2) = R(2) - DomLZ * roundint ( R(2) / DomLZ )
                if ( BC_X == 1 ) then
                        if ( R(0) .GT. DomLXHalf ) then
                                R(0) = R(0) - DomLX
                        else if ( R(0) .LT. - DomLXHalf ) then
                                R(0) = R(0) + DomLX
                        end if
                end if
                if ( BC_Y == 1 ) then
                        if ( R(1) .GT. DomLYHalf ) then
                                R(1) = R(1) - DomLY
                        else if ( R(1) .LT. - DomLYHalf ) then
                                R(1) = R(1) + DomLY
                        end if
                end if
                if ( BC_Z == 1 ) then
                        if ( R(2) .GT. DomLZHalf ) then
                                R(2) = R(2) - DomLZ
                        else if ( R(2) .LT. - DomLZHalf ) then
                                R(2) = R(2) + DomLZ
                        end if
                end if
        end subroutine ApplyPeriodicBC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine LinePoint ( Displacement, Q, R1, L1, R0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function calculates the point Q of projection of point R0 on line (R1,L1)
        ! Q = R1 + Disaplacement * L1
        !-------------------------------------------------------------------------------------------
        real*8, intent(inout)                   :: Displacement
        real*8, dimension(0:2), intent(inout)   :: Q
        real*8, dimension(0:2), intent(in)      :: R1, L1, R0
        !--------------------------------------------------------------------------------------------
                Q = R0 - R1
                ! Here we take into account periodic boundaries
                call ApplyPeriodicBC ( Q )
                Displacement = S_V3xV3 ( Q, L1 )
                Q = R1 + Displacement * L1
        end subroutine LinePoint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function LineLine ( H, cosA, D1, D2, L12, R1, L1, R2, L2, Prec ) !!!!!!!!!!!!!!!!!
        ! This function determines the neares distance H between two lines (R1,L1) and (R2,L2)
        !-------------------------------------------------------------------------------------------
        ! Input values:
        !      R1, L1, point and direction of line 1
        !      R2, L2, point and direction of line 2
        !      Prec, precision for the case L1 * L2 = 0 (parallel lines)
        ! Return values:
        !      H, minimal distance between lines
        !      cosA, cosine of angle between lines
        !      D1, D2, displacemets
        !      L12, unit vector directed along the closes distance
        !-------------------------------------------------------------------------------------------      
        real*8, intent(inout)                   :: H, cosA, D1, D2
        real*8, dimension(0:2), intent(out)     :: L12
        real*8, dimension(0:2), intent(in)      :: R1, L1, R2, L2
        !-------------------------------------------------------------------------------------------
        real*8, intent(in)                      :: Prec 
        real*8, dimension(0:2)                  :: Q1, Q2, R
        real*8                                  :: C, DD1, DD2, C1, C2
        !-------------------------------------------------------------------------------------------
                cosA = S_V3xV3 ( L1, L2 )
                C     = 1.0 - sqr ( cosA )
                if ( C < Prec ) then ! Lines are parallel to each other
                        LineLine = MD_LINES_PAR
                        return
                end if
                LineLine = MD_LINES_NONPAR
                R = R2 - R1
                ! Here we take into account periodic boundaries
                call ApplyPeriodicBC ( R )
                DD1 = S_V3xV3 ( R, L1 )
                DD2 = S_V3xV3 ( R, L2 )
                D1 = ( cosA * DD2 - DD1 ) / C
                D2 = ( DD2 - cosA * DD1 ) / C
                Q1 = R1 - D1 * L1
                Q2 = R2 - D2 * L2
                L12 = Q2 - Q1
                ! Here we take into account periodic boundaries
                call ApplyPeriodicBC ( L12 )
                H = S_V3norm3 ( L12 )
                if ( H < Prec ) then ! Lines intersect each other
                        C1 = signum ( D1 )
                        C2 = signum ( D1 ) * signum ( cosA )
                        Q1 = C1 * L1
                        Q2 = C2 * L2
                        call V3_V3xxV3 ( L12, Q1, Q2 )
                        call V3_ort ( L12 )
                else ! No intersection
                        L12 = L12 / H
                end if
        end function LineLine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module TPMGeom !********************************************************************************

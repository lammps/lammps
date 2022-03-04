module TMDGen3D !***********************************************************************************
!
! Generator of 3D CNT samples for TPM force field
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 13.00, 2020
!
!---------------------------------------------------------------------------------------------------

use TMDGenData

implicit none

contains !******************************************************************************************
        
        real*8 function MinimalDistance3D ( S1, S2, H, cosA, P, Q ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the minimum distance between two line segments in 3D
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)                     :: S1, S2
        real*8, intent(in)                      :: H, cosA
        real*8, dimension(0:1), intent(in)      :: P, Q
        !-------------------------------------------------------------------------------------------
        real*8                                  :: H2, cosA2, D
        real*8, dimension(0:1)                  :: P1, Q1
        integer*4, dimension(0:1,0:1)           :: KA
        integer*4                               :: i, j, K
        !-------------------------------------------------------------------------------------------
                if ( ( P(0) * P(1) .le. 0.0d+00 ) .and. ( Q(0) * Q(1) .le. 0.0d+00 ) ) then
                        MinimalDistance3D = H
                        S1 = 0.5d+00 * ( P(0) + P(1) )
                        S2 = 0.5d+00 * ( Q(0) + Q(1) )
                        return
                end  if
                do i = 0, 1 
                        P1(i) = P(i) * cosA
                        Q1(i) = Q(i) * cosA
                end do
                KA = 1
                K  = 0
                do i = 0, 1
                        if ( ( Q1(i) .ge. P(0) ) .and. ( Q1(i) .le. P(1) ) ) then
                                D = sqr ( Q(i) )
                                if ( K == 0 ) then
                                        MinimalDistance3D = D
                                        S1 = Q1(i)
                                        S2 = Q(i)
                                        K = 1
                                else if ( D < MinimalDistance3D ) then
                                        MinimalDistance3D = D
                                        S1 = Q1(i)
                                        S2 = Q(i)
                                end if
                                KA(0,i) = 0
                                KA(1,i) = 0
                        end if
                        if ( ( P1(i) .ge. Q(0) ) .and. ( P1(i) .le. Q(1) ) ) then
                                D = sqr ( P(i) )
                                if ( K == 0 ) then
                                        MinimalDistance3D = D
                                        S1 = P(i)
                                        S2 = P1(i)
                                        K = 1
                                else if ( D < MinimalDistance3D ) then
                                        MinimalDistance3D = D
                                        S1 = P(i)
                                        S2 = P1(i)
                                end if
                                KA(i,0) = 0
                                KA(i,1) = 0
                        end if
                end do
                H2 = sqr ( H )
                cosA2 = 2.0d+00 * cosA
                if ( K == 1 ) MinimalDistance3D = H2 + MinimalDistance3D * ( 1.0d+00 - sqr ( cosA ) )
                do i = 0, 1
                        do j = 0, 1
                                if ( KA(i,j) == 1 ) then
                                        D = H2 + sqr ( P(i) ) + sqr ( Q(j) ) - P(i) * Q(j) * cosA2 
                                        if ( K == 0 ) then
                                                MinimalDistance3D = D
                                                S1 = P(i)
                                                S2 = Q(j)
                                                K = 1
                                        else if ( D < MinimalDistance3D ) then
                                                MinimalDistance3D = D
                                                S1 = P(i)
                                                S2 = Q(j)
                                        end if
                                end if
                        end do
                end do
                MinimalDistance3D = dsqrt ( MinimalDistance3D )
        end function MinimalDistance3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
        subroutine RandTube3D ( X, L ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This subroutine generates a random tube in an isotropic 3D sample
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2), intent(out)     :: X, L
        !-------------------------------------------------------------------------------------------
        real*8                                  :: CT, ST, E
        !-------------------------------------------------------------------------------------------
                if ( BC_X0 == 0 ) then
                        X(0)= LS0 * randnumber () 
                else
                        X(0)= LS0 * ( 0.5d+00 - 1.0d+00 * randnumber () )
                end if
                if ( BC_Y0 == 0 ) then
                        X(1)= LS0 * randnumber ()
                else
                        X(1)= LS0 * ( 0.5d+00 - 1.0d+00 * randnumber () )
                end if
                if ( BC_Z0 == 0 ) then
                        X(2)= HS0 *randnumber () 
                else
                        X(2)= HS0 * ( 0.5d+00 - 1.0d+00 * randnumber () )
                end if
                CT  = 1.0d+00 - 2.0d+00 * randnumber () 
                ST  = sqrt ( 1.0d+00 - sqr ( CT ) )
                E   = M_2PI * randnumber ()
                L(0)= CT
                L(1)= ST * cos ( E )
                L(2)= ST * sin ( E )
        end subroutine RandTube3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        logical function AddTubeToSample3D ( MS ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function adds the last generated tube to the existing sample, if possible.
        ! In a case of periodic boundaries, this version is valid only f the tube length is smaller 
        ! than the half of the sample.
        !-------------------------------------------------------------------------------------------
        real*8, intent(inout)   :: MS
        !-------------------------------------------------------------------------------------------
        integer*4               :: i, m
        real*8                  :: Dmin, LT2, H, cosA, D1, D2, S1, S2
        real*8, dimension(0:2)  :: X, L12
        real*8, dimension(0:1)  :: P, Q
        !-------------------------------------------------------------------------------------------
                AddTubeToSample3D = .false.
                if ( .not. IsTubeInside ( NT ) ) return
                
                LT2 = 0.5d+00 * LT(NT)
                do m = 0, NT - 1
                        X = CT(NT,0:2)
                        if ( LineLine ( H, cosA, D1, D2, L12, X, DT(NT,0:2), CT(m,0:2), DT(m,0:2), GeomPrec ) == MD_LINES_NONPAR ) then
                                P(0) = D1 - LT2
                                P(1) = D1 + LT2
                                Q(0) = D2 - 0.5d+00 * LT(m)
                                Q(1) = D2 + 0.5d+00 * LT(m)
                                Dmin = MinimalDistance3D ( S1, S2, H, cosA, P, Q ) 
                        else
                                call LinePoint ( H, L12, CT(m,0:2), DT(m,0:2), X )
                                L12 = L12 - X
                                call ApplyPeriodicBC ( L12 )
                                Dmin = S_V3norm3 ( L12 )
                        end if
                        if ( Dmin < RT(NT) + RT(m) + DeltaT ) return
                end do

                MS = MS + TubeMass ( NT )
                NT = NT + 1
                AddTubeToSample3D = .true.
        end function AddTubeToSample3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine Generator3D () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This subroutine implements the whole fgenerator of 3D samples
        !-------------------------------------------------------------------------------------------
        integer*4       :: NA, NT0
        real*8          :: MS
        real*8          :: X1, X2, Y1, Y2, Z1, Z2
        !-------------------------------------------------------------------------------------------
                NT = 0
                MS = 0.0d+00
                NT0 = int ( MS0 / ( M_2PI * RT0 * LT0 * TPBM * TPBD ) )
                do 
                        if ( NT == MAX_TUBE ) then
                                print *, 'Error in [Generator3D]: MAX_TUBE is too small' 
                                stop
                        end if
                        if ( MS .ge. MS0 ) exit
                        NA = 0
                        ! Trying to add the tube to the sample
                        ! The maximal number of attempts is equal to NAmax
                        RT(NT) = RT0
                        LT(NT) = LT0 
                        do 
                                if ( NA == NAmax ) exit
                                call RandTube3D ( CT(NT,0:2), DT(NT,0:2) )
                                if ( AddTubeToSample3D ( MS ) ) then
                                        print '(a,i10,a,i10,a,i10)', 'Tube ', NT, '(Appr.', NT0, ' total): Attempt ', NA
                                        if ( BC_X0 == 0 ) then
                                                X1 = CT(NT,0) - 0.5d+00 * LT(NT) * DT(NT,0)
                                                X2 = CT(NT,0) + 0.5d+00 * LT(NT) * DT(NT,0)
                                                if ( DomXmin > X1 ) DomXmin = X1
                                                if ( DomXmin > X2 ) DomXmin = X2
                                                if ( DomXmax < X1 ) DomXmax = X1
                                                if ( DomXmax < X2 ) DomXmax = X2
                                        end if
                                        if ( BC_Y0 == 0 ) then
                                                Y1 = CT(NT,1) - 0.5d+00 * LT(NT) * DT(NT,1)
                                                Y2 = CT(NT,1) + 0.5d+00 * LT(NT) * DT(NT,1)
                                                if ( DomYmin > Y1 ) DomYmin = Y1
                                                if ( DomYmin > Y2 ) DomYmin = Y2
                                                if ( DomYmax < Y1 ) DomYmax = Y1
                                                if ( DomYmax < Y2 ) DomYmax = Y2
                                        end if
                                        if ( BC_Z0 == 0 ) then
                                                Z1 = CT(NT,2) - 0.5d+00 * LT(NT) * DT(NT,2)
                                                Z2 = CT(NT,2) + 0.5d+00 * LT(NT) * DT(NT,2)
                                                if ( DomZmin > Z1 ) DomZmin = Z1
                                                if ( DomZmin > Z2 ) DomZmin = Z2
                                                if ( DomZmax < Z1 ) DomZmax = Z1
                                                if ( DomZmax < Z2 ) DomZmax = Z2
                                        end if
                                        exit
                                end if
                                NA = NA + 1
                        end do
                end do
                MS0 = MS

                if ( BC_X0 == 0 ) DomLX = DomXmax - DomXmin
                if ( BC_Y0 == 0 ) DomLY = DomYmax - DomYmin
                if ( BC_Z0 == 0 ) DomLZ = DomZmax - DomZmin

                VS0 = ( DomXmax - DomXmin ) * ( DomYmax - DomYmin ) * ( DomZmax - DomZmin )
                DS0 = MS0 / VS0 * ( K_MDDU / 1.0d+03 )
        end subroutine Generator3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
end module TMDGen3D !*******************************************************************************

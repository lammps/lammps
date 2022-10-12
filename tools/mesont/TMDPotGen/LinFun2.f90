module LinFun2 !************************************************************************************
!
! Bi-linear functions and their derivatives.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!***************************************************************************************************

implicit none
        
contains !******************************************************************************************

        real*8 function CalcLinFun1_0 ( i, X, N, P, F ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer*4, intent(in)                           :: i, N
        real*8, intent(in)                              :: X
        real*8, dimension(0:N-1), intent(in)            :: P
        real*8, dimension(0:N-1), intent(inout)         :: F
        integer*4                                       :: i1
        real*8                                          :: A, A0
        !-------------------------------------------------------------------------------------------
                i1 = i - 1
                A0 = ( P(i) - X ) / ( P(i) - P(i1) )
                A  = 1.0d+00 - A0
                CalcLinFun1_0 = A0 * F(i1) + A * F(i)
        end function CalcLinFun1_0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine CalcLinFun1_1 ( S, Sx1, i, X, N, P, F, Fx ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)                             :: S, Sx1
        integer*4, intent(in)                           :: i, N
        real*8, intent(in)                              :: X
        real*8, dimension(0:N-1), intent(in)            :: P
        real*8, dimension(0:N-1), intent(inout)         :: F, Fx
        integer*4                                       :: i1
        real*8                                          :: A, A0
        !-------------------------------------------------------------------------------------------
                i1 = i - 1
                A0 = ( P(i) - X ) / ( P(i) - P(i1) )
                A  = 1.0d+00 - A0
                S = A0 * F(i1) + A * F(i)
                Sx1 = A0 * Fx(i1) + A * Fx(i)
        end subroutine CalcLinFun1_1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real*8 function CalcLinFun2_0 ( i, j, X, Y, N1, N2, P1, P2, F ) !!
        integer*4, intent(in)                           :: i, j, N1, N2
        real*8, intent(in)                              :: X, Y
        real*8, dimension(0:N1-1), intent(in)           :: P1
        real*8, dimension(0:N2-1), intent(in)           :: P2
        real*8, dimension(0:N1-1,0:N2-1), intent(inout) :: F
        integer*4                                       :: i1, j1
        real*8                                          :: A, A0, B, B0, G, G0
        !-------------------------------------------------------------------------------------------
                i1 = i - 1
                j1 = j - 1
                A0 = ( P1(i) - X ) / ( P1(i) - P1(i1) )
                A  = 1.0d+00 - A0
                B0 = ( P2(j) - Y ) / ( P2(j) - P2(j1) )
                B  = 1.0d+00 - B0
                G  = B0 * F(i,j1) + B * F(i,j)
                G0 = B0 * F(i1,j1) + B * F(i1,j)
                CalcLinFun2_0 = A0 * G0 + A * G
        end function CalcLinFun2_0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine CalcLinFun2_1 ( S, Sx1, Sy1, i, j, X, Y, N1, N2, P1, P2, F, Fx, Fy ) !!!!!!!!!!!!
        real*8, intent(out)                             :: S, Sx1, Sy1
        integer*4, intent(in)                           :: i, j, N1, N2
        real*8, intent(in)                              :: X, Y
        real*8, dimension(0:N1-1), intent(in)           :: P1
        real*8, dimension(0:N2-1), intent(in)           :: P2
        real*8, dimension(0:N1-1,0:N2-1), intent(inout) :: F, Fx, Fy
        integer*4                                       :: i1, j1
        real*8                                          :: A, A0, B, B0, G, G0
        !-------------------------------------------------------------------------------------------
                i1 = i - 1
                j1 = j - 1
                A0 = ( P1(i) - X ) / ( P1(i) - P1(i1) )
                A  = 1.0d+00 - A0
                B0 = ( P2(j) - Y ) / ( P2(j) - P2(j1) )
                B  = 1.0d+00 - B0
                
                G  = B0 * F(i,j1) + B * F(i,j)
                G0 = B0 * F(i1,j1) + B * F(i1,j)
                S = A0 * G0 + A * G
                
                G  = B0 * Fx(i,j1) + B * Fx(i,j)
                G0 = B0 * Fx(i1,j1) + B * Fx(i1,j)
                Sx1 = A0 * G0 + A * G
                
                G  = B0 * Fy(i,j1) + B * Fy(i,j)
                G0 = B0 * Fy(i1,j1) + B * Fy(i1,j)
                Sy1 = A0 * G0 + A * G
                
        end subroutine CalcLinFun2_1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module LinFun2 !********************************************************************************

module Spline2 !************************************************************************************
!
! TMD Library: Two-dimensional cubic spline function
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!***************************************************************************************************

use Spline1

implicit none

contains !******************************************************************************************

        subroutine CreateSpline2 ( CL, CD, CR, CU, N1, N2, N, P1, P2, F, Fxx, Fyy, Fxxyy, FF, MM, DD, K0, K1, K2 )
        integer*4, intent(in)                           :: CL, CD, CR, CU, N1, N2, N
        real*8, dimension(0:N1-1), intent(in)           :: P1
        real*8, dimension(0:N2-1), intent(in)           :: P2
        real*8, dimension(0:N1-1,0:N2-1), intent(inout) :: F, Fxx, Fyy, Fxxyy
        real*8, dimension(0:N-1), intent(inout)         :: FF, MM, DD, K0, K1, K2
        integer*4                                       :: II
        !-------------------------------------------------------------------------------------------
                do II = 0, N2 - 1
                        FF(0:N1-1) = F(0:N1-1,II)
                        MM(0)    = Fxx(0,II)
                        MM(N1-1) = Fxx(N1-1,II)
                        call CreateSpline1 ( CL, CR, N1, P1, FF, MM, DD, K0, K1, K2 )
                        Fxx(0:N1-1,II) = MM(0:N1-1)
                end do
                do II = 0, N1 - 1
                        MM(0) = Fyy(II,0)
                        MM(N-1) = Fyy(II,N2-1)
                        FF(0:N2-1) = F(II,0:N2-1)
                        call CreateSpline1 ( CD, CU, N2, P2, FF, MM, DD, K0, K1, K2 )
                        Fyy(II,0:N2-1) = MM(0:N2-1)
                end do
                FF(0:N1-1) = Fyy(0:N1-1,0 )
                call CreateSpline1 ( 3, 3, N1, P1, FF, MM, DD, K0, K1, K2 )
                Fxxyy(0:N1-1,0) = MM(0:N1-1)
                FF(0:N1-1) = Fyy(0:N1-1,N2-1 )
                call CreateSpline1 ( 3, 3, N1, P1, FF, MM, DD, K0, K1, k2 )
                Fxxyy(0:N1-1,N2-1) = MM(0:N1-1)
                do II = 1, N1 - 2
                        MM(0) = Fxxyy(II,0)
                        MM(N-1) = Fxxyy(II,N2-1)
                        FF(0:N2-1) = Fxx(II,0:N2-1)
                        call CreateSpline1 ( 2 , 2, N2, P2, FF, MM, DD, K0, K1, K2 )
                        Fxxyy(II,0:N2-1) = MM(0:N2-1)
                end do
        end subroutine CreateSpline2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CreateSpline2Ext ( CL, CD, CR, CU, N1, N1A, N2, N2A, N, P1, P2, F, Fxx, Fyy, Fxxyy, FF, MM, DD, K0, K1, K2 )
        integer*4, intent(in)                           :: CL, CD, CR, CU, N1, N1A, N2, N2A, N
        real*8, dimension(0:N1-1), intent(in)           :: P1
        real*8, dimension(0:N2-1), intent(in)           :: P2
        real*8, dimension(0:N1-1,0:N2-1), intent(inout) :: F, Fxx, Fyy, Fxxyy
        real*8, dimension(0:N-1), intent(inout)         :: FF, MM, DD, K0, K1, K2
        integer*4                                       :: II
        !-------------------------------------------------------------------------------------------
                Fxx = 0.0d+00
                Fyy = 0.0d+00
                Fxxyy = 0.0d+00
                
                do II = 0, N2A
                        FF(0:N1-1) = F(0:N1-1,II)
                        MM(0)    = Fxx(0,II)
                        MM(N1-1) = Fxx(N1-1,II)
                        call CreateSpline1 ( CL, CR, N1, P1, FF, MM, DD, K0, K1, K2 )
                        Fxx(0:N1-1,II) = MM(0:N1-1)
                end do
                
                do II = N2A + 1, N2 - 1
                        FF(0:N1-N1A-1) = F(N1A:N1-1,II)
                        MM(0) = Fxx(N1A,II)
                        MM(N1-N1A-1) = Fxx(N1-1,II)
                        call CreateSpline1 ( CL, CR, N1 - N1A, P1, FF, MM, DD, K0, K1, K2 )
                        Fxx(N1A:N1-1,II) = MM(0:N1-N1A-1)
                end do

                do II = 0, N1A - 1
                        MM(0) = Fyy(II,0)
                        MM(N2A) = Fyy(II,N2A)
                        FF(0:N2A) = F(II,0:N2A)
                        call CreateSpline1 ( CD, CU, N2A + 1, P2, FF, MM, DD, K0, K1, K2 )
                        Fyy(II,0:N2A) = MM(0:N2A)
                end do
                
                do II = N1A, N1 - 1
                        MM(0) = Fyy(II,0)
                        MM(N-1) = Fyy(II,N2-1)
                        FF(0:N2-1) = F(II,0:N2-1)
                        call CreateSpline1 ( CD, CU, N2, P2, FF, MM, DD, K0, K1, K2 )
                        Fyy(II,0:N2-1) = MM(0:N2-1)
                end do
                
                FF(0:N1-1) = Fyy(0:N1-1,0)
                call CreateSpline1 ( 3, 3, N1, P1, FF, MM, DD, K0, K1, K2 )
                Fxxyy(0:N1-1,0) = MM(0:N1-1)
                
                FF(0:N1A) = Fyy(0:N1A,N2A)
                call CreateSpline1 ( 3, 3, N1A + 1, P1, FF, MM, DD, K0, K1, K2 )
                Fxxyy(0:N1A,N2A) = MM(0:N1A)
                
                FF(0:N1-N1A-1) = Fyy(N1A:N1-1,N2-1 )
                call CreateSpline1 ( 3, 3, N1-N1A, P1, FF, MM, DD, K0, K1, K2 )
                Fxxyy(N1A:N1-1,N2-1) = MM(0:N1-N1A-1)
                
                do II = 1, N1A
                        MM(0) = Fxxyy(II,0)
                        MM(N2A) = Fxxyy(II,N2A)
                        FF(0:N2A) = Fxx(II,0:N2A)
                        call CreateSpline1 ( 2 , 2, N2A + 1, P2, FF, MM, DD, K0, K1, K2 )
                        Fxxyy(II,0:N2A) = MM(0:N2A)
                end do
                
                do II = N1A + 1, N1 - 2
                        MM(0) = Fxxyy(II,0)
                        MM(N-1) = Fxxyy(II,N2-1)
                        FF(0:N2-1) = Fxx(II,0:N2-1)
                        call CreateSpline1 ( 2 , 2, N2, P2, FF, MM, DD, K0, K1, K2 )
                        Fxxyy(II,0:N2-1) = MM(0:N2-1)
                end do
                
        end subroutine CreateSpline2Ext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        real*8 function CalcSpline2_0 ( i, j, X, Y, N1, N2, P1, P2, F, Fxx, Fyy, Fxxyy ) !!!!!!!!!!!
        integer*4, intent(in)                           :: i, j, N1, N2
        real*8, intent(in)                              :: X, Y
        real*8, dimension(0:N1-1), intent(in)           :: P1
        real*8, dimension(0:N2-1), intent(in)           :: P2
        real*8, dimension(0:N1-1,0:N2-1), intent(inout) :: F, Fxx, Fyy, Fxxyy
        integer*4                                       :: i1, j1
        real*8                                          :: T, Gy_0, Gy_1, Gxxy_0, Gxxy_1
        !-------------------------------------------------------------------------------------------
                i1 = i - 1
                j1 = j - 1
                T  = P2(j) - P2(j1)
                Gy_0    = ValueSpline1_0 ( Y, P2(j), P2(j1), F(i,j),    F(i,j1),    Fyy(i,j),    Fyy(i,j1),    T )
                Gy_1    = ValueSpline1_0 ( Y, P2(j), P2(j1), F(i1,j),   F(i1,j1),   Fyy(i1,j),   Fyy(i1,j1),   T )
                Gxxy_0  = ValueSpline1_0 ( Y, P2(j), P2(j1), Fxx(i,j),  Fxx(i,j1),  Fxxyy(i,j),  Fxxyy(i,j1),  T )
                Gxxy_1  = ValueSpline1_0 ( Y, P2(j), P2(j1), Fxx(i1,j), Fxx(i1,j1), Fxxyy(i1,j), Fxxyy(i1,j1), T )
                CalcSpline2_0 = ValueSpline1_0 ( X, P1(i), P1(i1), Gy_0, Gy_1,Gxxy_0, Gxxy_1, P1(i) - P1(i1) )
        end function CalcSpline2_0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CalcSpline2_1 ( S, Sx1, Sy1, i, j, X, Y, N1, N2, P1, P2, F, Fxx, Fyy, Fxxyy ) !!!
        real*8, intent(out)                             :: S, Sx1, Sy1
        integer*4, intent(in)                           :: i, j, N1, N2
        real*8, intent(in)                              :: X, Y
        real*8, dimension(0:N1-1), intent(in)           :: P1
        real*8, dimension(0:N2-1), intent(in)           :: P2
        real*8, dimension(0:N1-1,0:N2-1), intent(inout) :: F, Fxx, Fyy, Fxxyy
        integer*4                                       :: i1, j1
        real*8                                          :: T, Gy_0, Gy_1, Gxxy_0, Gxxy_1
        real*8                                          :: Gyy_0, Gyy_1, Gxxyy_0, Gxxyy_1
        !-------------------------------------------------------------------------------------------
                i1 = i - 1
                j1 = j - 1
                T  = P2(j) - P2(j1)
                call ValueSpline1_1 ( Gy_0, Gyy_0, Y, P2(j), P2(j1), F(i,j),    F(i,j1),    Fyy(i,j),    Fyy(i,j1),    T )
                call ValueSpline1_1 ( Gy_1, Gyy_1, Y, P2(j), P2(j1), F(i1,j),   F(i1,j1),   Fyy(i1,j),   Fyy(i1,j1),   T )
                call ValueSpline1_1 ( Gxxy_0, Gxxyy_0, Y, P2(j), P2(j1), Fxx(i,j),  Fxx(i,j1),  Fxxyy(i,j),  Fxxyy(i,j1),  T )
                call ValueSpline1_1 ( Gxxy_1, Gxxyy_1, Y, P2(j), P2(j1), Fxx(i1,j), Fxx(i1,j1), Fxxyy(i1,j), Fxxyy(i1,j1), T )
                call ValueSpline1_1 ( S, Sx1, X, P1(i), P1(i1), Gy_0, Gy_1,Gxxy_0, Gxxy_1, P1(i) - P1(i1) )
                Sy1 = ValueSpline1_0 ( X, P1(i), P1(i1), Gyy_0, Gyy_1,Gxxyy_0, Gxxyy_1, P1(i) - P1(i1) )
        end subroutine CalcSpline2_1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module Spline2 !********************************************************************************

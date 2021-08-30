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

module Spline1 !************************************************************************************
!
! One-dimensional cubic spline function.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!***************************************************************************************************
use iso_c_binding, only : c_int, c_double, c_char
implicit none

contains !******************************************************************************************

        real(c_double) function ValueSpline1_0 ( X, Xi, Xi_1, Yi, Yi_1, Mi, Mi_1, Hi_1 ) !!!!!!!!!!!!
        real(c_double), intent(in)      :: X, Xi, Xi_1, Yi, Yi_1, Mi, Mi_1, Hi_1
        real(c_double)                  :: H26, HL, HR
        !-------------------------------------------------------------------------------------------
                H26     = Hi_1 * Hi_1 / 6.0
                Hl      = X - Xi_1
                Hr      = Xi - X
                ValueSpline1_0 = ( ( Mi_1 * Hr * Hr * Hr + Mi * Hl * Hl * Hl ) / 6.0 + ( Yi_1 - Mi_1 * H26 ) * Hr &
                        + ( Yi - Mi * H26 ) * Hl ) / Hi_1
        end function ValueSpline1_0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine ValueSpline1_1 ( S, S1, X, Xi, Xi_1, Yi, Yi_1, Mi, Mi_1, Hi_1 ) !!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: S, S1
        real(c_double), intent(in)      :: X, Xi, Xi_1, Yi, Yi_1, Mi, Mi_1, Hi_1
        real(c_double)                  :: H6, H26, HL, HR, HL2, HR2
        !-------------------------------------------------------------------------------------------
                H6  = Hi_1 / 6.0d+00
                H26 = Hi_1 * H6
                HL  = X - Xi_1
                HR  = Xi - X
                HL2 = HL * HL
                HR2 = HR * HR
                S   = ( ( Mi_1 * HR2 * Hr + Mi * HL2 * Hl ) / 6.0 + ( Yi_1 - Mi_1 * H26 ) * HR + ( Yi - Mi * H26 ) * HL ) / Hi_1
                S1  = ( ( Mi * HL2 - Mi_1 * HR2 ) / 2.0d+00 + Yi - Yi_1 ) / Hi_1 - H6 * ( Mi - Mi_1 )
        end subroutine ValueSpline1_1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine sprogonka3 ( N, K0, K1, K2, F, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(c_int), intent(in)                   :: N
        real(c_double), dimension(0:N-1), intent(in)    :: K0, K1, K2
        real(c_double), dimension(0:N-1), intent(inout) :: F, X
        real(c_double)                                  :: D
        integer(c_int)                               :: i
        !-------------------------------------------------------------------------------------------
                X(0) = F(0) / K1(0)
                F(0) = - K2(0) / K1(0)
                do i = 1, N - 1
                        D    = - ( K1(i) + F(i-1) * K0(i) )
                        X(i) = ( K0(i) * X(i-1) - F(i) ) / D
                        F(i) = K2(i) / D
                end do
                do i = N - 2, 0, -1
                        X(i) = X(i) + F(i) * X(i+1)
                end do
        end subroutine sprogonka3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CreateSpline1 ( CL, CR, N, P, F, M, D, K0, K1, K2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(c_int), intent(in)                   :: CL, CR, N
        real(c_double), dimension (0:N-1), intent(in)   :: P, F
        real(c_double), dimension (0:N-1), intent(inout):: M, D, K0, K1, K2
        integer(c_int)                               :: i
        real(c_double)                                  :: Z
        !-------------------------------------------------------------------------------------------
                do i = 1, N - 1
                        K0(i) = P(i) - P(i-1)
                        K1(i) = ( F(i) - F(i-1) ) / K0(i)
                end do
                select case ( CL )
                        case (1)
                                K1(0) = 2.0d+00 / 3.0d+00
                                K2(0) = 1.0d+00 / 3.0d+00
                                D (0) = 2 * ( K1(1) - M(0) ) / K0(1)
                        case (2)
                                K1(0) = 1.0d+00
                                K2(0) = 0.0d+00
                                D(0)  = M(0)
                        case (3)
                                K1(0) = 1.0d+00
                                K2(0) = 0.0d+00
                                D(0)  = 0.0d+00
                end select
                Z = K1(N-1)
                do i = 1, N - 2
                        D(i)  = 6.0d+00 * ( K1(i+1) - K1(i) )
                        K2(i) = K0(i+1)
                        K1(i) = 2.0d+00 * ( K2(i) + K0(i) )
                end do
                select case ( CR )
                        case (1)
                                D(N-1)  = 2.0d+00 * ( M(N-1) - Z ) / K0(N-1)
                                K1(N-1) = 2.0d+00 / 3.0d+00
                                K0(N-1) = 1.0d+00 / 3.0d+00
                        case (2)
                                K1(N-1) = 1.0d+00
                                K0(N-1) = 0.0d+00
                                D(N-1)  = M(N-1)
                        case (3)
                                K1(N-1) = 1.0d+00
                                K0(N-1) = 0.0d+00
                                D(N-1)  = 0.0d+00
                end select
                call sprogonka3 ( N, K0, K1, K2, D, M )
        end subroutine CreateSpline1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(c_double) function CalcSpline1_0 ( i, X, N, P, F, M ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(c_int), intent(in)                   :: i, N
        real(c_double), intent(in)                      :: X
        real(c_double), dimension(0:N-1), intent(in)    :: P, F, M
        integer(c_int)                               :: j
        real(c_double)                                  :: HL, HR, H, H6, H26, HR2, HL2, HRH, HLH
        !-------------------------------------------------------------------------------------------
                j   = i - 1
                HL  = X - P(j)
                HR  = P(i) - X
                H   = P(i) - P(j)
                H6  = H / 6.0d+00
                H26 = H * H6
                HL2 = HL * HL
                HR2 = HR * HR
                HLH = HL / H
                HRH = HR / H
                CalcSpline1_0 = ( M(j) * HR2 * HRH + M(i) * HL2 * HLH ) / 6.0d+00 + ( F(j) - M(j) * H26 ) * HRH &
                        + ( F(i) - M(i) * H26 ) * HLH
        end function CalcSpline1_0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CalcSpline1_1 ( S, S1, i, X, N, P, F, M ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)                     :: S, S1
        integer(c_int), intent(in)                   :: i, N
        real(c_double), intent(in)                      :: X
        real(c_double), dimension(0:N-1), intent(in)    :: P, F, M
        integer(c_int)                               :: j
        real(c_double)                                  :: HL, HR, H, H6, H26, HR2, HL2, HRH, HLH
        !-------------------------------------------------------------------------------------------
                j   = i - 1
                HL  = X - P(j)
                HR  = P(i) - X
                H   = P(i) - P(j)
                H6  = H / 6.0d+00
                H26 = H * H6
                HL2 = HL * HL
                HR2 = HR * HR
                HLH = HL / H
                HRH = HR / H
                S  = ( M(j) * HR2 * HRH + M(i) * HL2 * HLH ) / 6.0d+00 + ( F(j) - M(j) * H26 ) * HRH + ( F(i) - M(i) * H26 ) * HLH
                S1 = ( ( M(i) * HL2 - M(j) * HR2 ) / 2.0d+00 + F(i) - F(j) ) / H - H6 * ( M(i) - M(j) )
        end subroutine CalcSpline1_1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CalcSpline1_2 ( S, S1, S2, i, X, N, P, F, M ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)                     :: S, S1, S2
        integer(c_int), intent(in)                   :: i, N
        real(c_double), intent(in)                      :: X
        real(c_double), dimension(0:N-1), intent(in)    :: P, F, M
        integer(c_int)                               :: j
        real(c_double)                                  :: HL, HR, H, H6, H26, HR2, HL2, HRH, HLH
        !-------------------------------------------------------------------------------------------
                j   = i - 1
                HL  = X - P(j)
                HR  = P(i) - X
                H   = P(i) - P(j)
                H6  = H / 6.0d+00
                H26 = H * H6
                HL2 = HL * HL
                HR2 = HR * HR
                HLH = HL / H
                HRH = HR / H
                S  = ( M(j) * HR2 * HRH + M(i) * HL2 * HLH ) / 6.0d+00 + ( F(j) - M(j) * H26 ) * HRH + ( F(i) - M(i) * H26 ) * HLH
                S1 = ( ( M(i) * HL2 - M(j) * HR2 ) / 2.0d+00 + F(i) - F(j) ) / H - H6 * ( M(i) - M(j) )
                S2 = M(j) * HRH + M(i) * HLH
        end subroutine CalcSpline1_2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Spline1 !********************************************************************************

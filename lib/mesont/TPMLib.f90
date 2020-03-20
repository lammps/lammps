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

module TPMLib !*************************************************************************************
!
! TMD Library: Basic constants, types, and mathematical functions
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

!---------------------------------------------------------------------------------------------------
! Mathematical constants
!---------------------------------------------------------------------------------------------------

        real(c_double), parameter  :: M_PI_2            = 1.57079632679489661923
        real(c_double), parameter  :: M_PI              = 3.14159265358979323846
        real(c_double), parameter  :: M_3PI_2           = 4.71238898038468985769
        real(c_double), parameter  :: M_2PI             = 6.28318530717958647692
        real(c_double), parameter  :: M_PI_180          = 0.017453292519943295769

!---------------------------------------------------------------------------------------------------
! Physical unit constants
!---------------------------------------------------------------------------------------------------

        real(c_double), parameter  :: K_AMU             = 1.66056E-27                           ! a.m.u. (atomic mass unit, Dalton)
        real(c_double), parameter  :: K_EV              = 1.60217646e-19                        ! eV (electron-volt)

        real(c_double), parameter  :: K_MDLU            = 1.0E-10                               ! MD length unit (m)
        real(c_double), parameter  :: K_MDEU            = K_EV                                  ! MD energy unit (J)
        real(c_double), parameter  :: K_MDMU            = K_AMU                                 ! MD mass unit (kg)
        real(c_double), parameter  :: K_MDFU            = K_MDEU / K_MDLU                       ! MD force unit (N)
        real(c_double), parameter  :: K_MDCU            = K_MDEU / K_MDMU                       ! MD specific heat unit (J/(kg*K))

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------

        integer(c_int)          :: StdUID            = 31

contains !******************************************************************************************

!---------------------------------------------------------------------------------------------------
! Simple mathematical functions
!---------------------------------------------------------------------------------------------------
        
        real(c_double) function rad ( X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: X
        !-------------------------------------------------------------------------------------------
                rad = X * M_PI_180
        end function rad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(c_double) function sqr ( X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: X
        !-------------------------------------------------------------------------------------------
                sqr = X * X
        end function sqr !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function signum ( X )  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: X
        !-------------------------------------------------------------------------------------------   
                if ( X > 0 ) then
                        signum = 1
                else if ( X < 0 ) then
                        signum = -1
                else
                        signum = 0
                end if
        end function signum !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Vector & matrix functions
!---------------------------------------------------------------------------------------------------

        real(c_double) function S_V3xx ( V ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), dimension(0:2), intent(in)      :: V
        !-------------------------------------------------------------------------------------------
                S_V3xx = V(0) * V(0) + V(1) * V(1) + V(2) * V(2)
        end function S_V3xx !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(c_double) function S_V3xV3 ( V1, V2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), dimension(0:2), intent(in)      :: V1, V2
        !-------------------------------------------------------------------------------------------
                S_V3xV3 = V1(0) * V2(0) + V1(1) * V2(1) + V1(2) * V2(2)
        end function S_V3xV3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(c_double) function S_V3norm3 ( V ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), dimension(0:2), intent(in)      :: V
        !-------------------------------------------------------------------------------------------
                S_V3norm3 = dsqrt ( V(0) * V(0) + V(1) * V(1) + V(2) * V(2) )
        end function S_V3norm3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine V3_ort ( V ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Vector production
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2), intent(inout)   :: V
        !-------------------------------------------------------------------------------------------
        real(c_double)                                  :: Vabs 
        !-------------------------------------------------------------------------------------------
                Vabs = S_V3norm3 ( V )
                V(0) = V(0) / Vabs
                V(1) = V(1) / Vabs
                V(2) = V(2) / Vabs
        end subroutine V3_ort !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine V3_V3xxV3 ( V, V1, V2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Vector production
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2), intent(out)     :: V
        real(c_double), dimension(0:2), intent(in)      :: V1, V2
        !-------------------------------------------------------------------------------------------
                V(0) = V1(1) * V2(2) - V1(2) * V2(1)
                V(1) = V1(2) * V2(0) - V1(0) * V2(2)
                V(2) = V1(0) * V2(1) - V1(1) * V2(0)
        end subroutine V3_V3xxV3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Handling of spherical and Euler angles
!---------------------------------------------------------------------------------------------------
        
        subroutine RotationMatrix3  ( M, Psi, Tet, Phi ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Ksi, Tet and Phi are Euler angles
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2,0:2), intent(out) :: M
        real(c_double), intent(in)                      :: Psi, Tet, Phi
        !-------------------------------------------------------------------------------------------
        real(c_double)                                  :: cK, sK, cT, sT, cP, sP
        !-------------------------------------------------------------------------------------------
                cK = dcos ( Psi )
                sK = dsin ( Psi )
                cT = dcos ( Tet )
                sT = dsin ( Tet )
                cP = dcos ( Phi )
                sP = dsin ( Phi )
                M(0,0) = cP * cK - sK * sP * cT
                M(0,1) = cP * sK + sP * cT * cK
                M(0,2) = sP * sT
                M(1,0) = - sP * cK - cP * cT * sK
                M(1,1) = - sP * sK + cP * cT * cK
                M(1,2) = cP * sT
                M(2,0) = sT * sK
                M(2,1) = - sT * cK
                M(2,2) = cT
        end subroutine RotationMatrix3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine EulerAngles ( Psi, Tet, L ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)                     :: Tet, Psi
        real(c_double), dimension(0:2), intent(in)      :: L
        !-------------------------------------------------------------------------------------------
                Tet = acos ( L(2) )
                Psi = atan2 ( L(1), L(0) )
                if ( Psi > M_3PI_2 ) then
                        Psi = Psi - M_3PI_2
                else 
                        Psi = Psi + M_PI_2
                end if
        end subroutine EulerAngles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! File inout and output
!---------------------------------------------------------------------------------------------------

        integer(c_int) function OpenFile ( Name, Params, Path ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        character*(*), intent(in)       :: Name, Params, Path
        !-------------------------------------------------------------------------------------------
        integer(c_int)                       :: Fuid
        character*512                   :: FullName, Msg, Name1, Action1, Status1, Form1, Position1
        !-------------------------------------------------------------------------------------------
                OpenFile = StdUID + 1
                if ( Params(1:1) == 'r' ) then
                        Action1   = 'read'
                        Status1   = 'old'
                        Position1 = 'rewind'
                else if ( Params(1:1) == 'w' ) then
                        Action1 = 'write'
                        Status1 = 'replace'
                        Position1 = 'rewind'
                else if ( Params(1:1) == 'a' ) then
                        Action1 = 'write'
                        Status1 = 'old'
                        Position1 = 'append'
                endif
                if ( Params(2:2) == 'b' ) then
                        Form1 = 'binary'
                else
                        Form1 = 'formatted'
                endif
                open ( unit = OpenFile, file = Name, form = Form1, action = Action1, status = Status1, position = Position1 )
                StdUID = StdUID + 1
                return
        end function OpenFile !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CloseFile ( Fuid ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(c_int), intent(inout)        :: Fuid
        !-------------------------------------------------------------------------------------------
                if ( Fuid < 0 ) return
                close ( unit = Fuid )
                Fuid = -1
        end subroutine CloseFile !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module TPMLib !*********************************************************************************

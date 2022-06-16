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

module TubePotBase !********************************************************************************
!
! Non-bonded pair interaction potential and transfer functions for atoms composing nanotubes.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 09.01, 2017
!
!---------------------------------------------------------------------------------------------------
!
! This module contains basic parameters for all modules involved into calculations of tubular
! potentials.
!
! It includes definitions of
!   -- TPBU, Lennard-Jones (12-6) potential
!   -- TPBQ, Transfer function
!
! All default values are adjusted for non-bonded carbon-carbon interaction in carbon nanotubes.
!
!***************************************************************************************************

use TPMLib
use iso_c_binding, only : c_int, c_double, c_char
implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------

        ! Types of the potential with respect to the breathing mode
        integer(c_int), parameter            :: TP_POT_MONO_R           = 0
        integer(c_int), parameter            :: TP_POT_POLY_R           = 1

        ! Maximal number of elements in corresponding tables
        integer(c_int), parameter            :: TPBNMAX                 = 2001

        ! Numerical constants
        real(c_double), parameter               :: TPbConstD            = 5.196152422706632d+00 ! = 3.0**1.5

        ! Mass of C atom
        real(c_double), parameter               :: TPBMc                = 12.0107d+00            ! (Da)

        ! Parameters of the Van der Waals interaction between carbon atoms in graphene sheets, see
        ! Stuart S.J., Tutein A.B., Harrison J.A., J. Chem. Phys. 112(14), 2000
        real(c_double), parameter               :: TPBEcc               = 0.00284d+00           ! (eV)
        real(c_double), parameter               :: TPBScc               = 3.4d+00               ! (A)

        ! Lattice parameter and surface number density of atoms for a graphene sheet, see
        ! Dresselhaus et al, Carbon 33(7), 1995
        real(c_double), parameter               :: TPBAcc               = 1.421d+00             ! (A)
        real(c_double), parameter               :: TPBDcc               = 4.0d+00 / ( TPBConstD * TPBAcc * TPBAcc ) ! (1/A^2)

        ! Specific heat of carbon nanotubes
        real(c_double), parameter               :: TPBSHcc              = 600.0d+00 / K_MDCU    ! (eV/(Da*K))

        ! Cutoff distances for the interactomic potential and transfer function.
        ! Changes in these parameters can result in necessity to change some numerical parameters too.
        real(c_double), parameter               :: TPBRmincc            = 0.001d+00 * TPBScc    ! (A)
        real(c_double), parameter               :: TPBRcutoffcc         = 3.0d+00 * TPBScc      ! (A)
        real(c_double), parameter               :: TPBRcutoff1cc        = 2.16d+00 * TPBScc     ! (A)

        ! Parameters of the transfer function for non-bonded interaction between carbon atoms
        real(c_double), parameter               :: TPBQScc              = 7.0d+00               ! (A)
        real(c_double), parameter               :: TPBQRcutoff1cc       = 8.0d+00               ! (A)

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------

        ! Set to .true. to generate diagnostic and warning messages
        logical         :: TPErrCheck           = .true.
        character*512   :: TPErrMsg             = ''

        real(c_double)  :: TPGeomPrec           = 1.0d-06       ! Geometric precision, see TPInt
        integer(c_int)  :: TPPotType            = TP_POT_MONO_R ! Type of the potential with respect to the breathing mode

        ! Parameters of the interatomic potential and atoms distribution at the surface
        ! of the tube

        real(c_double)  :: TPBM                 = TPBMc          ! Mass of an atom (Da)
        real(c_double)  :: TPBE                 = TPBEcc         ! Depth of the energy well in (12-6) LJ interatomic potential (eV)
        real(c_double)  :: TPBS                 = TPBScc         ! Sigma parameter of (12-6) LJ interatomic potential (A)
        real(c_double)  :: TPBD                 = TPBDcc         ! Numerical density of atoms at the tube surface (1/A^2)
        real(c_double)  :: TPBSH                = TPBSHcc        ! Specific heat (eV/(Da*K))

        real(c_double)  :: TPBRmin              = TPBRmincc      ! (A)
        real(c_double)  :: TPBRcutoff           = TPBRcutoffcc   ! (A)
        real(c_double)  :: TPBRcutoff1          = TPBRcutoff1cc  ! (A)

        ! Parameters of the transfer function

        real(c_double)  :: TPBQS                = TPBQScc       ! Sigma parameter of the transfer function (A)
        real(c_double)  :: TPBQRcutoff1         = TPBQRcutoff1cc! (A)

        ! Auxiliary variables

        real(c_double)  :: TPBE4, TPBE24, TPBDRcutoff, TPBQDRcutoff
        real(c_double)  :: TPBQR0                               ! Constant-value distance for the transfer function (A)

        ! Table of inter-particle potential, force, and transfer function

        integer(c_int)                          :: TPBN                 = TPBNMAX
        real(c_double)                          :: TPBDR
        real(c_double), dimension(0:TPBNMAX-1)  :: TPBQ
        real(c_double), dimension(0:TPBNMAX-1)  :: TPBU, TPBdUdR

contains !******************************************************************************************

        integer(c_int) function TPBsizeof () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                TPBsizeof = 8 * ( size ( TPBQ ) + size ( TPBU ) + size ( TPBdUdR ) )
        end function TPBsizeof !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Interpolation
!---------------------------------------------------------------------------------------------------

        real(c_double) function TPBQInt0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real(c_double)                  :: Z, RR
        integer(c_int)                  :: i
        !-------------------------------------------------------------------------------------------
                if ( R < TPBRmin ) then
                        !call PrintStdLogMsg ( TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10)' ) ': R < Rmin: R=', R, ', Rmin=', TPBRmin
                        !call Error ( 'TPBQInt0', TPErrMsg )
                elseif ( R > TPBRcutoff ) then
                        TPBQInt0 = 0.0d+00
                        return
                endif
                RR   = ( R - TPBRmin ) / TPBDR
                i    = int ( RR )
                RR   = RR - i
                Z    = 1.0d+00 - RR
                TPBQInt0 = TPBQ(i) * Z + TPBQ(i+1) * RR
        end function TPBQInt0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(c_double) function TPBUInt0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real(c_double)                  :: Z, RR
        integer(c_int)                  :: i
        !-------------------------------------------------------------------------------------------
                if ( R < TPBRmin ) then
                        !call PrintStdLogMsg ( TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10)' ) ': R < Rmin: R=', R, ', Rmin=', TPBRmin
                        !call Error ( 'TPBUInt0', TPErrMsg )
                elseif ( R > TPBRcutoff ) then
                        TPBUInt0 = 0.0d+00
                        return
                endif
                RR   = ( R - TPBRmin ) / TPBDR
                i    = int ( RR )
                RR   = RR - i
                Z    = 1.0d+00 - RR
                TPBUInt0 = TPBU(i) * Z + TPBU(i+1) * RR
        end function TPBUInt0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPBUInt1 ( U, dUdR, R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: U, dUdR
        real(c_double), intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real(c_double)                  :: Z, RR
        integer(c_int)                  :: i
        !-------------------------------------------------------------------------------------------
                if ( R < TPBRmin ) then
                        !call PrintStdLogMsg ( TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10)' ) ': R < Rmin: R=', R, ', Rmin=', TPBRmin
                        !call Error ( 'TPBUInt1', TPErrMsg )
                elseif ( R > TPBRcutoff ) then
                        TPBU = 0.0d+00
                        TPBdUdR = 0.0d+00
                        return
                endif
                RR   = ( R - TPBRmin ) / TPBDR
                i    = int ( RR )
                RR   = RR - i
                Z    = 1.0d+00 - RR
                U    = TPBU(i) * Z + TPBU(i+1) * RR
                dUdR = TPBdUdR(i) * Z + TPBdUdR(i+1) * RR
        end subroutine TPBUInt1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Calculation
!---------------------------------------------------------------------------------------------------

        real(c_double) function TPBQCalc0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real(c_double)                  :: Z, t, S
        !-------------------------------------------------------------------------------------------
                if ( R > TPBRcutoff ) then
                        TPBQCalc0 = 0.0d+00
                else if ( R < TPBQR0 ) then
                        TPBQCalc0 = 1.0d+00
                else
                        Z = TPBQS / R
                        Z = Z * Z * Z
                        Z = Z * Z
                        TPBQCalc0 = 4.0d+00 * ( 1.0d+00 - Z ) * Z
                        if ( R > TPBQRcutoff1 ) then
                                t = ( R - TPBQRcutoff1 ) / TPBQDRcutoff
                                S = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                TPBQCalc0 = TPBQCalc0 * S
                        endif
                endif
        end function TPBQCalc0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(c_double) function TPBUCalc0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real(c_double)                  :: Z, t, S
        !-------------------------------------------------------------------------------------------
                if ( R > TPBRcutoff ) then
                        TPBUCalc0 = 0.0d+00
                else
                        Z = TPBS / R
                        Z = Z * Z * Z
                        Z = Z * Z
                        TPBUCalc0 = TPBE4 * ( Z - 1.0d+00 ) * Z
                        if ( R > TPBRcutoff1 ) then
                                t = ( R - TPBRcutoff1 ) / TPBDRcutoff
                                S = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                TPBUCalc0 = TPBUCalc0 * S
                        endif
                endif
        end function TPBUCalc0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPBUCalc1 ( U, dUdR, R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: U, dUdR
        real(c_double), intent(in)      :: R
        real(c_double)                  :: Z, t, S, dSdR
        !-------------------------------------------------------------------------------------------
                if ( R > TPBRcutoff ) then
                        U    = 0.0d+00
                        dUdR = 0.0d+00
                else
                        Z    = TPBS / R
                        Z    = Z * Z * Z
                        Z    = Z * Z
                        U    = TPBE4 * ( Z - 1.0d+00 ) * Z
                        dUdR = TPBE24 * ( 2.0d+00 * Z - 1.0d+00 ) * Z / R
                        if ( R > TPBRcutoff1 ) then
                                t    = ( R - TPBRcutoff1 ) / TPBDRcutoff
                                S    = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                dSdR = 6.0d+00 * t * ( t - 1.0d+00 ) / TPBDRcutoff
                                dUdR = dUdR * S + U * dSdR
                                U     = U * S
                        endif
                endif
        end subroutine TPBUCalc1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPBSegmentForces ( F1, F2, F, M, Laxis, L ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), dimension(0:2), intent(out)     :: F1, F2
        real(c_double), dimension(0:2), intent(in)      :: F, M, Laxis
        real(c_double), intent(in)                      :: L
        !-------------------------------------------------------------------------------------------
        real(c_double), dimension(0:2)                  :: FF, MM, FFF
        !-------------------------------------------------------------------------------------------
                FF = 0.5d+00 * F
                MM = M / L
                call V3_V3xxV3 ( FFF, MM, Laxis )
                F1 = FF - FFF
                F2 = FF + FFF
        end subroutine TPBSegmentForces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Initialization
!---------------------------------------------------------------------------------------------------

        subroutine TPBInit () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double)          :: R
        integer(c_int)       :: i
        !-------------------------------------------------------------------------------------------
                TPBE4 = 4.0d+00 * TPBE
                TPBE24 = - 24.0d+00 * TPBE
                TPBDRcutoff = TPBRcutoff - TPBRcutoff1
                TPBQDRcutoff = TPBRcutoff - TPBQRcutoff1
                TPBQR0 = TPBQS * 2.0d+00 ** ( 1.0d+00 / 6.0d+00 )
                TPBDR = ( TPBRcutoff - TPBRmin ) / ( TPBN - 1 )
                R = TPBRmin
                do i = 0, TPBN - 1
                        TPBQ(i) = TPBQCalc0 ( R )
                        call TPBUCalc1 ( TPBU(i), TPBdUdR(i), R )
                        R = R + TPBDR
                enddo
        end subroutine TPBInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TubePotBase !****************************************************************************

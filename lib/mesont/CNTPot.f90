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

module CNTPot !*************************************************************************************
!
! Mesoscopic potential for internal modes in CNTs.
!
!---------------------------------------------------------------------------------------------------
!
! Carbon nanotubes internal potentials:
!       CNTSTRH0, harmonic stretching potential of type 0 with constant Young's modulus
!       CNTSTRH1, harmonic stretching potential of type 1 with variable Young's modulus
!       CNTSTRNH0, non-harmonic stretching with fracture potential of type 0
!       CNTSTRNH1, non-harmonic stretching with fracture potential of type 1
!       CNTBNDH, harmonic bending potential
!       CNTBNDHB, harmonic bending-buckling potential
!       CNTBNDHBF, harmonic bending-buckling potential with fracture
!       CNTTRS, torsion potential
!       CNTBRT, breathing potential
!
! The functional form and force constants of harmonic stretching, bending and
! torsion potentials are taken from:
! L.V. Zhigilei, Ch. Wei, D. Srivastava, Phys. Rev. B 71, 165417 (2005)
!
! The model of stress-strain curve for the non-harmonic potential with fracture
! is developed and parameterized using the following constants
! -- Young's modulus
! -- maximum linear strain (only for the NH potential of type 1)
! -- tensile strength (or fracture strain)
! -- strain at failure (or fracture strain)
! -- maximum strain.
! All these parameters are assumed to be independent of CNT radius or chriality type.
! In this model, the true strain at failure CNTSTREft and true tensile strength
! CNTSTRSft are slightly different from the imposed values CNTSTREf and CNTSTRSf.
!
! The non-harmonic stretching potentials of types 0 and 1 are different from
! each other by the functional form of the stress-strain curve
!
! Different parameterizations of CNTSTRH0, CNTSTRNH0 and CNTSTRNH1 potentials can be chosen,
! see subroutine CNTSTRSetParameterization
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 13.00, 2020
!
!***************************************************************************************************

use TPMLib
use iso_c_binding, only : c_int, c_double, c_char
implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------
        integer(c_int), parameter            :: CNTPOT_STRETCHING    = 0
        integer(c_int), parameter            :: CNTPOT_SBUCKLING     = 1
        integer(c_int), parameter            :: CNTPOT_SFRACTURE     = 2

        integer(c_int), parameter            :: CNTPOT_BENDING       = 3
        integer(c_int), parameter            :: CNTPOT_BBUCKLING     = 4
        integer(c_int), parameter            :: CNTPOT_BFRACTURE     = 5

        ! Harmonic stretching model (constant Young's modulus)
        integer(c_int), parameter            :: CNTSTRMODEL_H0       = 0
        ! Harmonic stretching model (Young's modulus depends on radius)
        integer(c_int), parameter            :: CNTSTRMODEL_H1       = 1
        ! Non-harmonic stretching with fracture, potential of type 0
        integer(c_int), parameter            :: CNTSTRMODEL_NH0F     = 2
        ! Non-harmonic stretching without fracture, potential of type 1
        integer(c_int), parameter            :: CNTSTRMODEL_NH1      = 3
        ! Non-harmonic stretching with fracture, potential of type 1
        integer(c_int), parameter            :: CNTSTRMODEL_NH1F     = 4
        ! Harmonic stretching model + axial buckling
        integer(c_int), parameter            :: CNTSTRMODEL_H1B      = 5
        ! Harmonic stretching model + axial buckling + hysteresis
        integer(c_int), parameter            :: CNTSTRMODEL_H1BH     = 6

        integer(c_int), parameter            :: CNTBNDMODEL_H        = 0     ! Harmonic bending model
        integer(c_int), parameter            :: CNTBNDMODEL_HB       = 1     ! Harmonic bending - buckling model
        integer(c_int), parameter            :: CNTBNDMODEL_HBF      = 2     ! Harmonic bending - buckling - fracture model
        integer(c_int), parameter            :: CNTBNDMODEL_HBH      = 3     ! Harmonic bending - buckling + Hysteresis

        integer(c_int), parameter            :: CNTPOTNMAX           = 4000  ! Maximum number of points in the interpolation tables

!---------------------------------------------------------------------------------------------------
! Parameters of potentials
!---------------------------------------------------------------------------------------------------

        ! Stretching potential

        ! Type of the bending model
        integer(c_int)                       :: CNTSTRModel  = CNTSTRMODEL_H1
        ! Type of parameterization
        integer(c_int)                       :: CNTSTRParams = 0
        ! Type of dependence of the Young's modulus on tube radius
        integer(c_int)                       :: CNTSTRYMT = 0

        ! Parameters of non-harmonic potential and fracture model
        real(c_double)                          :: CNTSTRR0     = 6.8d+00       ! Reference radius of nanotubes (A)
                                                                                ! (this parameter is not used for the model
                                                                                ! parametrization, but only for calculation of the
                                                                                ! force constant in eV/A)
        real(c_double)                          :: CNTSTRD0     = 3.4d+00       ! CNT wall thickness (A)
        real(c_double)                          :: CNTSTREmin   = -0.4d+00      ! Minimum strain in tabulated potential
        real(c_double)                          :: CNTSTREmax   = 0.13d+00      ! Maximum strain in tabulated potential.
                                                                                ! Simultaneously, U=0 if E> CNTSTREmax
        real(c_double)                          :: CNTSTREl     = 5.0d-02       ! Maximum linear strain
        real(c_double)                          :: CNTSTREf     = 12.0d-02      ! Strain at failure
        real(c_double)                          :: CNTSTRS0     = 0.850e+12     ! Young's modulus (Pa)
        real(c_double)                          :: CNTSTRSl                     ! Maximum linear stress (Pa)
        real(c_double)                          :: CNTSTRSf     = 75.0d+09      ! Tensile strength (Pa)
        real(c_double)                          :: CNTSTRF0                     ! Elastic force constant (eV/A**2)
        real(c_double)                          :: CNTSTRFl                     ! Maximal linear force, (eV/A**2)
        real(c_double)                          :: CNTSTRFf                     ! Tensile force at failure (eV/A**2)
        real(c_double)                          :: CNTSTRSi                     ! Maximum stress (not used in the model) (Pa)
        real(c_double)                          :: CNTSTRDf                     ! dF/dE at failure

        real(c_double)                          :: CNTSTRAA, CNTSTRBB           !
        real(c_double)                          :: CNTSTRAAA, CNTSTRBBB         ! Auxiliary constants
        real(c_double)                          :: CNTSTRUl, CNTSTRUf           !

        ! Axial buckling - hysteresis approach
        real(c_double)                          :: CNTSTREc     = -0.0142d+00   ! The minimum buckling strain
        real(c_double)                          :: CNTSTREc1    = -0.04d+00     ! Critical axial buckling strain
        real(c_double)                          :: CNTSTREc2    = -0.45d+00     ! Maximum buckling strain

        ! Bending potential

        integer(c_int)                       :: CNTBNDModel  = CNTBNDMODEL_H    ! Type of the bending model
        ! Buckling model parameters
        real(c_double)                          :: CNTBNDN      = 1.0d+00       ! Buckling exponent
        real(c_double)                          :: CNTBNDB      = 0.68d+00      ! Buckling number
        real(c_double)                          :: CNTBNDR      = 275.0d+00     ! Critical radius of curvature (A)
                                                                                ! This is the mean value for (10,10) SWCNT
        real(c_double)                          :: CNTBNDTF     = M_PI * 120.0d+00 / 180.0d+00 ! Fracture buckling angle (rad)
        real(c_double)                          :: CNTBNDN1
        real(c_double)                          :: CNTBNDC2

contains !******************************************************************************************

!---------------------------------------------------------------------------------------------------
! Stretching potential
!---------------------------------------------------------------------------------------------------

        subroutine CNTSTRSetParameterization ( PType ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This subroutine setups parameters for further parameterization of stretching models
        ! References:
        ! [1] Yu M.-F. et al., Phys. Rev. Lett. 84(24), 5552 (2000)
        ! [2] Liew K.M. et al., Acta Materialia 52, 2521 (2004)
        ! [3] Mielke S.L. et al., Chem. Phys. Lett. 390, 413 (2004)
        ! [4] Zhigilei L.V. et al., Phys. Rev. B 71, 165417 (2005)
        ! [5] Kelly B.T., Physics of graphite, 1981
        !-------------------------------------------------------------------------------------------
        integer(c_int), intent(in) :: PType
        !-------------------------------------------------------------------------------------------
                select case ( PType )
                        case ( 0 ) ! This parametrization is based on averaged exp. data of Ref. [1]
                                CNTSTRR0     = 6.8d+00       ! Ref. [1]
                                CNTSTRD0     = 3.4d+00       ! Ref. [1]
                                CNTSTREmin   = -0.4d+00      ! Chosen arbitrary
                                CNTSTREmax   = 3.64d-02      ! = CNTSTREf + 0.005
                                CNTSTREl     = 2.0d-02       ! Chosen arbitrary
                                CNTSTREf     = 3.14d-02      ! Ref. [1]
                                CNTSTRS0     = 1.002e+12     ! Ref. [1]
                                CNTSTRSf     = 30.0d+09      ! Ref. [1]
                        case ( 1 ) ! This parameterization is taken from Ref. [2] for (10,10) CNTs.
                                   ! These values are obtained in MD simulations with REBO potential.
                                   ! Values of Young's modulus, tensile strength and stress here
                                   ! are close to those obtained in Ref. [3] for pristine (defectless)
                                   ! (5,5) CNT in semi-empirical QM calculations based on PM3 model
                                CNTSTRR0     = 6.785d+00     ! Calculated with the usual formula for (10,10) CNT
                                CNTSTRD0     = 3.35d+00      ! Ref. [2]
                                CNTSTREmin   = -0.4d+00      ! Chosen arbitrary
                                CNTSTREmax   = 28.4d-02      ! = CNTSTREf + 0.005
                                CNTSTREl     = 5.94d-02      ! Ref. [2]
                                CNTSTREf     = 27.9d-02      ! Corresponds to maximum strain in Ref. [2]
                                CNTSTRS0     = 1.031e+12     ! Ref. [2]
                                CNTSTRSf     = 148.5d+09     ! Corresponds to tensile strength in Ref. [2]
                        case ( 2 ) ! This parametrization is taken from Ref. [3] for (5,5) CNTs
                                   ! with one atom vacancy defect obtained with the semi-empirical QM PM3 model
                                CNTSTRR0     = 3.43d+00      ! Ref. [3]
                                CNTSTRD0     = 3.4d+00       ! Ref. [3]
                                CNTSTREmin   = -0.4d+00      ! Chosen arbitrary
                                CNTSTREmax   = 15.8d-02      ! = CNTSTREf + 0.005
                                CNTSTREl     = 6.00d-02      ! Chosen similar to Ref. [2]
                                CNTSTREf     = 15.3d-02      ! Ref. [3]
                                CNTSTRS0     = 1.100e+12     ! Ref. [3]
                                CNTSTRSf     = 100.0d+09     ! Ref. [3]
                        case ( 3 ) ! This special parameterization changes only the value of Young's modulus
                                   ! in accordance with the stretching constant in Ref. [4]
                                CNTSTRS0     = ( 86.64d+00 + 100.56d+00 * CNTSTRR0 ) * K_MDFU &
                                        / ( M_2PI * CNTSTRR0 * CNTSTRD0 * 1.0d-20 ) ! Ref. [4]
                        case ( 4 ) ! This special parameterization changes only the value of Young's modulus
                                   ! making it equal to the in-plane Young's modulus of graphite
                                CNTSTRR0     = 6.785d+00     ! Calculated with the usual formula for (10,10) CNT
                                CNTSTRD0     = 3.4d+00       ! Ref. [1]
                                CNTSTRS0     = 1.06e+12      ! Ref. [5]
                end select
        end subroutine CNTSTRSetParameterization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Stretching without fracture, harmonic potential
        !

        integer(c_int) function CNTSTRH0Calc ( U, dUdL, L, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Young's modulus is independent of R.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)     :: U, dUdL
        real(c_double), intent(in)      :: L, R0, L0
        real(c_double)                  :: E
        !-------------------------------------------------------------------------------------------
                E    = ( L - L0 ) / L0
                dUdL = R0 * CNTSTRF0 * E
                U    = 0.5d+00 * L0 * E * dUdL
                CNTSTRH0Calc = CNTPOT_STRETCHING
        end function CNTSTRH0Calc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function CNTSTRH1Calc ( U, dUdL, L, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Young's modulus depends on R, see [4].
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)     :: U, dUdL
        real(c_double), intent(in)      :: L, R0, L0
        real(c_double)                  :: E, K
        !-------------------------------------------------------------------------------------------
                E    = ( L - L0 ) / L0
                K    = 86.64d+00 + 100.56d+00 * R0
                dUdL = K * E
                U    = 0.5d+00 * L0 * E * dUdL
                CNTSTRH1Calc = CNTPOT_STRETCHING
        end function CNTSTRH1Calc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Stretching without fracture, harmonic potential, with axial buckling without hysteresis
        !

        integer(c_int) function CNTSTRH1BCalc ( U, dUdL, L, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Young's modulus depends on R, see [4].
        ! Axial buckling without hysteresis.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)     :: U, dUdL
        real(c_double), intent(in)      :: L, R0, L0
        real(c_double)                  :: E, K, Kbcl, dUbcl, d, ud
        !-------------------------------------------------------------------------------------------
                E    = ( L - L0 ) / L0
                K    = 86.64d+00 + 100.56d+00 * R0
                Kbcl = -10.98d+00 * L0
                if ( E .gt. CNTSTREc ) then ! Harmonic stretching
                        dUdL = K * E
                        U    = 0.5d+00 * L0 * E * dUdL
                        CNTSTRH1BCalc = CNTPOT_STRETCHING
                else if ( E .gt. CNTSTREc2 ) then ! Axial buckling
                        dUbcl = 0.5d+00 * L0 * K * CNTSTREc * CNTSTREc - Kbcl * CNTSTREc
                        U = Kbcl * E + dUbcl
                        dUdL = Kbcl / L0
                        CNTSTRH1BCalc = CNTPOT_STRETCHING
                else    ! Return to harmonic potential
                        d = -0.0142794
                        dUdL = K * ( d + E - CNTSTREc2 )
                        dUbcl = 0.5d+00 * L0 * K * CNTSTREc * CNTSTREc - Kbcl * CNTSTREc + Kbcl * CNTSTREc2
                        Ud = 0.5d+00 * L0 * K * d * d
                        U = 0.5d+00 * L0 *  (d+E-CNTSTREc2) * dUdL + dUbcl - Ud
                        CNTSTRH1BCalc = CNTPOT_STRETCHING
                end if
        end function CNTSTRH1BCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Stretching without fracture, harmonic potential, with axial buckling with hysteresis
        !

        integer(c_int) function CNTSTRH1BHCalc ( U, dUdL, L, R0, L0, ABF, Ebuc ) !!!!!!!!!!!!!!!!!!!
        ! Young's modulus depends on R, see [4]
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)     :: U, dUdL, Ebuc
        real(c_double), intent(in)      :: L, R0, L0
        integer(c_int), intent(in)   :: ABF
        !-------------------------------------------------------------------------------------------
        real(c_double)                  :: E, K, dUbcl, Ebcl, Kbcl, Edu
        real(c_double)                  :: C, DE, t
        !-------------------------------------------------------------------------------------------
                E = ( L - L0 ) / L0
                K = 86.64d+00 + 100.56d+00 * R0
                Kbcl = -10.98d+00 * L0
                if ( E .gt. CNTSTREc ) then ! Harmonic potential - no buckling
                        dUdL = K * E
                        U = 0.5d+00 * L0 * E * dUdL
                        CNTSTRH1BHCalc = CNTPOT_STRETCHING
                        Ebuc = 0.0d+00
                else if ( E .gt. CNTSTREc1 ) then  ! Above minimal buckling strain, but not at critical strain
                        if ( ABF .eq. 0 ) then ! Not buckled. Continue harmonic potential
                                dUdL = K * E
                                U = 0.5d+00 * L0 * E * dUdL
                                CNTSTRH1BHCalc = CNTPOT_STRETCHING
                                Ebuc = 0.0d+00
                        else   ! Relaxing from buckled state. Use buckling potential
                                dUbcl = 0.5d+00 * L0 * K * CNTSTREc * CNTSTREc - Kbcl * CNTSTREc
                                U = Kbcl * E + dUbcl
                                dUdL = Kbcl / L0
                                CNTSTRH1BHCalc = CNTPOT_SBUCKLING
                                Ebuc = 0.0d+00
                        end if
                else if( E .gt. CNTSTREc2 ) then ! Axial buckling strain region
                        if ( ABF .eq. 0 ) then ! Newly buckled
                                dUbcl = 0.5d+00 * L0 * K * CNTSTREc * CNTSTREc - Kbcl * CNTSTREc
                                U = Kbcl * E + dUbcl
                                dUdL = Kbcl / L0
                                CNTSTRH1BHCalc = CNTPOT_SBUCKLING
                                Ebuc = 0.5d+00 * L0 * K * CNTSTREc1 * CNTSTREc1 - Kbcl * CNTSTREc1 - dUbcl
                        else ! Already buckled
                                dUbcl = 0.5d+00 * L0 * K * CNTSTREc * CNTSTREc - Kbcl * CNTSTREc
                                U = Kbcl * E + dUbcl
                                dUdL = Kbcl / L0
                                CNTSTRH1BHCalc = CNTPOT_SBUCKLING
                                Ebuc = 0.0d+00
                        end if
                else  ! Maximum strain and return to harmonic potential
                        dUdL = K * E
                        U = 0.5d+00 * L0 * E * dUdL
                        CNTSTRH1BHCalc = CNTPOT_STRETCHING
                        Ebuc = 0.0d+00
                end if
        end function CNTSTRH1BHCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Stretching with fracture, non-harmonic potential of type 0
        !

        integer(c_int) function CNTSTRNH0FCalc ( U, dUdL, L, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: U, dUdL
        real(c_double), intent(in)      :: L, R0, L0
        real(c_double)                  :: E, DE, t
        !-------------------------------------------------------------------------------------------
                E = ( L - L0 ) / L0
                if ( E < CNTSTREf ) then
                        dUdL = ( CNTSTRAA - CNTSTRBB * E ) * E
                        U    = ( CNTSTRAAA - CNTSTRBBB * E ) * E * E
                        CNTSTRNH0FCalc = CNTPOT_STRETCHING
                else
                        dUdL = 0.0d+00
                        U    = 0.0d+00
                        CNTSTRNH0FCalc = CNTPOT_SFRACTURE
                end if
                U    = L0 * R0 * U
                dUdL = R0 * dUdL
        end function CNTSTRNH0FCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CNTSTRNH0Init () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double)          :: S
        !-------------------------------------------------------------------------------------------
                S        = M_2PI * CNTSTRD0 * 1.0e-20 / K_MDFU
                CNTSTRSl = CNTSTRS0 * CNTSTREl
                CNTSTRF0 = CNTSTRS0 * S
                CNTSTRFl = CNTSTRSl * S
                CNTSTRFf = CNTSTRSf * S
                CNTSTRAA = CNTSTRF0
                CNTSTRBB = ( CNTSTRF0 * CNTSTREf - CNTSTRFf ) / ( CNTSTREf * CNTSTREf )
                CNTSTRAAA= CNTSTRAA / 2.0d+00
                CNTSTRBBB= CNTSTRAA / 3.0d+00
                CNTSTRUl = 0.0d+00
                CNTSTRUf = ( CNTSTRAAA - CNTSTRBBB * CNTSTREf ) * CNTSTREf * CNTSTREf
                ! These two values are not defined yet
                CNTSTRSi = 0.0d+00
                CNTSTRDf = 0.0d+00
        end subroutine CNTSTRNH0Init !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Stretching without fracture, non-harmonic potential of type 1
        !

        integer(c_int) function CNTSTRNH1Calc ( U, dUdL, L, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: U, dUdL
        real(c_double), intent(in)      :: L, R0, L0
        real(c_double)                  :: E, C, DE, t
        !-------------------------------------------------------------------------------------------
                E = ( L - L0 ) / L0
                if ( E < CNTSTREl ) then
                        dUdL = CNTSTRF0 * E
                        U    = 0.5d+00 * E * dUdL
                        CNTSTRNH1Calc = CNTPOT_STRETCHING
                else
                        DE   = E - CNTSTREl
                        C    = 1.0 + CNTSTRBB * DE
                        dUdL = CNTSTRFl + CNTSTRAA * ( 1.0d+00 - 1.0d+00 / C )
                        U    = CNTSTRUl + CNTSTRAAA * DE - CNTSTRBBB * dlog ( C )
                end if
                CNTSTRNH1Calc = CNTPOT_STRETCHING
                U    = L0 * R0 * U
                dUdL = R0 * dUdL
        end function CNTSTRNH1Calc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Stretching with fracture, non-harmonic potential of type 1
        !

        integer(c_int) function CNTSTRNH1FCalc ( U, dUdL, L, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: U, dUdL
        real(c_double), intent(in)      :: L, R0, L0
        real(c_double)                  :: E, C, DE, t
        !-------------------------------------------------------------------------------------------
                E = ( L - L0 ) / L0
                if ( E < CNTSTREl ) then
                        dUdL = CNTSTRF0 * E
                        U    = 0.5d+00 * E * dUdL
                        CNTSTRNH1FCalc = CNTPOT_STRETCHING
                else if ( E < CNTSTREf ) then
                        DE   = E - CNTSTREl
                        C    = 1.0 + CNTSTRBB * DE
                        dUdL = CNTSTRFl + CNTSTRAA * ( 1.0d+00 - 1.0d+00 / C )
                        U    = CNTSTRUl + CNTSTRAAA * DE - CNTSTRBBB * dlog ( C )
                        CNTSTRNH1FCalc = CNTPOT_STRETCHING
                else
                        dUdL = 0.0d+00
                        U    = 0.0d+00
                        CNTSTRNH1FCalc = CNTPOT_SFRACTURE
                end if
                U    = L0 * R0 * U
                dUdL = R0 * dUdL
        end function CNTSTRNH1FCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CNTSTRNH1Init () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double)          :: S, C, E, t
        integer(c_int)       :: i, CaseID
        !-------------------------------------------------------------------------------------------
                S        = M_2PI * CNTSTRD0 * 1.0e-20 / K_MDFU
                CNTSTRSl = CNTSTRS0 * CNTSTREl
                CNTSTRF0 = CNTSTRS0 * S
                CNTSTRFl = CNTSTRSl * S
                CNTSTRFf = CNTSTRSf * S
                CNTSTRAA = ( CNTSTRFf - CNTSTRFl ) * ( CNTSTREf * CNTSTRF0 - CNTSTRFl ) / ( CNTSTREf * CNTSTRF0 - CNTSTRFf )
                CNTSTRBB = CNTSTRF0 / CNTSTRAA
                CNTSTRAAA= CNTSTRFl + CNTSTRAA
                CNTSTRBBB= CNTSTRAA / CNTSTRBB
                CNTSTRSi = CNTSTRSl + CNTSTRAA / S
                C        =  1.0 + CNTSTRBB * ( CNTSTREf - CNTSTREl )
                CNTSTRDf = CNTSTRF0 / C / C
                CNTSTRUl = 0.5d+00 * CNTSTRFl * CNTSTREl
                CNTSTRUf = CNTSTRUl + ( CNTSTRFl + CNTSTRAA ) * ( CNTSTREf - CNTSTREl ) - CNTSTRAA * dlog ( C ) / CNTSTRBB
        end subroutine CNTSTRNH1Init !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! General
        !

        integer(c_int) function CNTSTRCalc ( U, dUdL, L, R0, L0 , ABF, Ebuc ) !!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)     :: U, dUdL, Ebuc
        real(c_double), intent(in)      :: L, R0, L0
        integer(c_int), intent(in)      :: ABF
        !-------------------------------------------------------------------------------------------
                Ebuc = 0.0d+00
                select case ( CNTSTRModel )
                        case ( CNTSTRMODEL_H0 )
                                CNTSTRCalc = CNTSTRH0Calc ( U, dUdL, L, R0, L0 )
                        case ( CNTSTRMODEL_H1 )
                                CNTSTRCalc = CNTSTRH1Calc ( U, dUdL, L, R0, L0 )
                        case ( CNTSTRMODEL_NH0F )
                                CNTSTRCalc = CNTSTRNH0FCalc ( U, dUdL, L, R0, L0 )
                        case ( CNTSTRMODEL_NH1 )
                                CNTSTRCalc = CNTSTRNH1Calc ( U, dUdL, L, R0, L0 )
                        case ( CNTSTRMODEL_NH1F )
                                CNTSTRCalc = CNTSTRNH1FCalc ( U, dUdL, L, R0, L0 )
                        case ( CNTSTRMODEL_H1B )
                                CNTSTRCalc = CNTSTRH1BCalc ( U, dUdL, L, R0, L0 )
                        case ( CNTSTRMODEL_H1BH )
                                CNTSTRCalc = CNTSTRH1BHCalc ( U, dUdL, L, R0, L0, ABF, Ebuc )
                end select
        end function CNTSTRCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CNTSTRInit ( STRModel, STRParams, YMType, Rref ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(c_int), intent(in)   :: STRModel, STRParams, YMType
        real(c_double), intent(in)      :: Rref
        !-------------------------------------------------------------------------------------------
                CNTSTRModel = STRModel
                CNTSTRParams = STRParams
                CNTSTRYMT = YMType
                if ( STRModel .ne. CNTSTRMODEL_H1 ) then
                        call CNTSTRSetParameterization ( STRParams )
                        if ( YMType == 2 ) then
                                call CNTSTRSetParameterization ( 4 )
                        else if ( YMType == 1 ) then
                                CNTSTRR0 = Rref
                                call CNTSTRSetParameterization ( 3 )
                        end if
                        if ( STRModel == CNTSTRMODEL_NH0F ) then
                                call CNTSTRNH0Init ()
                        else
                                call CNTSTRNH1Init ()
                        end if
                end if
        end subroutine CNTSTRInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Bending potentials
!---------------------------------------------------------------------------------------------------

        subroutine BendingGradients ( K, G0, G1, G2, R0, R1, R2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(inout)                   :: K
        real(c_double), dimension(0:2), intent(inout)   :: G0, G1, G2
        real(c_double), dimension(0:2), intent(in)      :: R0, R1, R2
        real(c_double), dimension(0:2)                  :: DR0, DR2
        real(c_double)                                  :: L0, L2
        !-------------------------------------------------------------------------------------------
                DR0 = R0 - R1
                DR2 = R2 - R1
                L0 = S_V3norm3 ( DR0 )
                L2 = S_V3norm3 ( DR2 )
                DR0 = DR0 / L0
                DR2 = DR2 / L2
                K = S_V3xV3 ( DR0, DR2 )
                G0 = DR2 - K * DR0
                G2 = DR0 - K * DR2
                G0 = G0 / L0
                G2 = G2 / L2
                G1 = - ( G0 + G2 )
        end subroutine BendingGradients !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function CNTBNDHCalc ( U, dUdC, C, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Bending model of type 0:Harmonic bending potential.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)             :: U, dUdC
        real(c_double), intent(in)              :: C, R0, L0
        real(c_double)                          :: E, K
        !-------------------------------------------------------------------------------------------
                E = 1.0d+00 - C
                K  = 2.0d+00 * ( 63.8d+00 * R0**2.93d+00 ) / L0
                U    = K * ( 1.0d+00 + C ) / E
                dUdC = 2.0d+00 * K / ( E * E )
                CNTBNDHCalc = CNTPOT_BENDING
        end function CNTBNDHCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function CNTBNDHBCalc ( U, dUdC, C, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Bending model of type 1: Harmonic bending potential with buckling.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)             :: U, dUdC
        real(c_double), intent(in)              :: C, R0, L0
        real(c_double)                          :: E1, E2, C2, Kbnd, Kbcl, Theta, DUbcl
        !-------------------------------------------------------------------------------------------
                E1 = 1.0d+00 - C
                E2 = 1.0d+00 + C
                ! Calculate the square of curvature
                C2 = 4.0d+00 * E2 / ( L0 * L0 * E1 )
                ! Check the condition for buckling
                if ( C2 .ge. CNTBNDC2 ) then ! Buckling takes place
                        Theta= M_PI - acos ( C )
                        Kbnd = 63.8d+00 * R0**2.93d+00
                        Kbcl = CNTBNDB * Kbnd / CNTBNDR
                        DUbcl= Kbnd * ( CNTBNDB * ( M_PI - 2.0d+00 * atan ( 2.0 * CNTBNDR / L0 ) ) - 0.5d+00 * L0 / CNTBNDR ) &
                                / CNTBNDR
                        U    = Kbcl * abs( Theta )**CNTBNDN - DUbcl
                        dUdC = Kbcl * CNTBNDN * abs( Theta )**CNTBNDN1 / sqrt ( 1.0d+00 - C * C )
                        CNTBNDHBCalc = CNTPOT_BBUCKLING
                else ! Harmonic bending
                        Kbnd = 2.0d+00 * ( 63.8d+00 * R0**2.93d+00 ) / L0
                        U    = Kbnd * E2 / E1
                        dUdC = 2.0d+00 * Kbnd / ( E1 * E1 )
                        CNTBNDHBCalc = CNTPOT_BENDING
                end if
        end function CNTBNDHBCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function CNTBNDHBFCalc ( U, dUdC, C, R0, L0 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)             :: U, dUdC
        real(c_double), intent(in)              :: C, R0, L0
        real(c_double)                          :: E1, E2, C2, Kbnd, Kbcl, Theta, DUbcl
        !-------------------------------------------------------------------------------------------
                E1 = 1.0d+00 - C
                E2 = 1.0d+00 + C
                ! Calculate the square of curvature
                C2 = 4.0d+00 * E2 / ( L0 * L0 * E1 )
                ! Check the condition for buckling
                if ( C2 .ge. CNTBNDC2 ) then ! Buckling takes place
                        Theta= M_PI - acos ( C )
                        if ( Theta > CNTBNDTF ) then ! Fracture takes place
                                U    = 0.0d+00
                                dUdC = 0.0d+00
                                CNTBNDHBFCalc = CNTPOT_BFRACTURE
                        else
                                Kbnd = 63.8d+00 * R0**2.93d+00
                                Kbcl = CNTBNDB * Kbnd / CNTBNDR
                                DUbcl= Kbnd * ( CNTBNDB * ( M_PI - 2.0d+00 * atan ( 2.0 * CNTBNDR / L0 ) ) - &
                                  0.5d+00 * L0 / CNTBNDR ) / CNTBNDR
                                U    = Kbcl * abs ( Theta )**CNTBNDN - DUbcl
                                dUdC = Kbcl * CNTBNDN * abs ( Theta )**CNTBNDN1 / sqrt ( 1.0d+00 - C * C )
                                CNTBNDHBFCalc = CNTPOT_BBUCKLING
                        end if
                else ! Harmonic bending
                        Kbnd = 2.0d+00 * ( 63.8d+00 * R0**2.93d+00 ) / L0
                        U    = Kbnd * E2 / E1
                        dUdC = 2.0d+00 * Kbnd / ( E1 * E1 )
                        CNTBNDHBFCalc = CNTPOT_BENDING
                end if
        end function CNTBNDHBFCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer(c_int) function CNTBNDHBHCalc ( U, dUdC, C, R0, L0, BBF, Ebuc ) !!!!!!!!!!!!!!!!!!!!
        ! Bending model of type 1: Harmonic bending potential with buckling with hysteresis approach.
        !-------------------------------------------------------------------------------------------
        real(c_double), intent(out)             :: U, dUdC, Ebuc
        real(c_double), intent(in)              :: C , R0, L0
        integer(c_int), intent(in)           :: BBF
        real(c_double)                          :: E1, E2, C2, Kbnd, Kbcl,Theta,DUbcl, Ubcl, Cmin,Rmax
        !-------------------------------------------------------------------------------------------
                Rmax = 340.0d+00
                Cmin = 1.0/(Rmax*Rmax)
                E1 = 1.0d+00 - C
                E2 = 1.0d+00 + C
                ! Calculate the square of curvature
                C2 = 4.0d+00 * E2 / ( L0 * L0 * E1 )
                Theta = M_PI - acos ( C )
                if ( C2 .lt. Cmin ) then  ! Harmonic bending
                        Kbnd = 2.0d+00 * ( 63.8d+00 * R0**2.93d+00 ) / L0
                        U    = Kbnd * E2 / E1
                        dUdC = 2.0d+00 * Kbnd / ( E1 * E1 )
                        CNTBNDHBHCalc = CNTPOT_BENDING
                        Ebuc = 0.0
                else if ( C2 .ge. Cmin .and. C2 .lt. CNTBNDC2 ) then  ! Potential depends on buckling flag of a node
                        if ( BBF .eq. 0 ) then ! Not buckled yet. Continue harmonic bending
                                Kbnd = 2.0d+00 * ( 63.8d+00 * R0**2.93d+00 ) / L0
                                U    = Kbnd * E2 / E1
                                dUdC = 2.0d+00 * Kbnd / ( E1 * E1 )
                                CNTBNDHBHCalc = CNTPOT_BENDING
                                Ebuc = 0.0d+00
                        else ! Already has been buckled or is buckled. Use buckling potential until Cmin.
                                Theta= M_PI - acos ( C )
                                Kbnd = 63.8d+00 * R0**2.93d+00
                                Kbcl = CNTBNDB * Kbnd / CNTBNDR
                                DUbcl= 2.0d+00*Kbnd * &
                                        (1.0d+00+cos(l0/Rmax+M_PI))/(1.0d+00-cos(l0/Rmax+M_PI))/L0-Kbcl*abs(l0/Rmax)**CNTBNDN
                                U    = Kbcl * abs( Theta )**CNTBNDN + DUbcl
                                dUdC = Kbcl * CNTBNDN * abs( Theta )**CNTBNDN1 / sqrt ( 1.0d+00 - C * C )
                                Ebuc = 0.0d+00
                                CNTBNDHBHCalc = CNTPOT_BBUCKLING
                        end if
                   else ! Greater than buckling critical point
                        if ( BBF .eq. 1 ) then ! Already buckled
                                Theta= M_PI - acos ( C )
                                Kbnd = 63.8d+00 * R0**2.93d+00
                                Kbcl = CNTBNDB * Kbnd / CNTBNDR
                                DUbcl= 2.0d+00*Kbnd * &
                                        (1.0d+00+cos(l0/Rmax+M_PI))/(1.0d+00-cos(l0/Rmax+M_PI))/L0-Kbcl*abs(l0/Rmax)**CNTBNDN
                                U    = Kbcl * abs( Theta )**CNTBNDN + DUbcl
                                dUdC = Kbcl * CNTBNDN * abs( Theta )**CNTBNDN1 / sqrt ( 1.0d+00 - C * C )
                                Ebuc = 0.0d00
                                CNTBNDHBHCalc = CNTPOT_BBUCKLING
                        else ! Newly buckled
                                Theta= M_PI - acos ( C )
                                Kbnd = 63.8d+00 * R0**2.93d+00
                                Kbcl = CNTBNDB * Kbnd / CNTBNDR
                                DUbcl= 2.0d+00*Kbnd * &
                                        (1.0d+00+cos(l0/Rmax+M_PI))/(1.0d+00-cos(l0/Rmax+M_PI))/L0-Kbcl*abs(l0/Rmax)**CNTBNDN
                                U    = Kbcl * abs( Theta )**CNTBNDN + DUbcl
                                dUdC = Kbcl * CNTBNDN * abs( Theta )**CNTBNDN1 / sqrt ( 1.0d+00 - C * C )
                                Ebuc = 2.0d+00*Kbnd * (1.0d+00+cos(l0/CNTBNDR+M_PI)) / (1.0d+00-cos(l0/CNTBNDR+M_PI))/L0 &
                                        - Kbcl * abs ( l0 / CNTBNDR ) ** CNTBNDN - dUbcl
                                CNTBNDHBHCalc = CNTPOT_BBUCKLING
                        end if
                end if
        end function CNTBNDHBHCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! General
        !

        integer(c_int) function CNTBNDCalc ( U, dUdC, C, R0, L0, BBF, Ebuc ) !!!!!!!!!!!!!!!!!!!!!!!
        real(c_double), intent(out)             :: U, dUdC, Ebuc
        real(c_double), intent(in)              :: C, R0, L0
        integer(c_int), intent(in)           :: BBF
        !-------------------------------------------------------------------------------------------
                Ebuc = 0.0d+00
                select case ( CNTBNDModel )
                        case ( CNTBNDMODEL_H )
                                CNTBNDCalc = CNTBNDHCalc ( U, dUdC, C, R0, L0 )
                        case ( CNTBNDMODEL_HB )
                                CNTBNDCalc = CNTBNDHBCalc ( U, dUdC, C, R0, L0 )
                        case ( CNTBNDMODEL_HBF )
                                CNTBNDCalc = CNTBNDHBFCalc ( U, dUdC, C, R0, L0 )
                        case ( CNTBNDMODEL_HBH )
                                CNTBNDCalc = CNTBNDHBHCalc ( U, dUdC, C, R0, L0, BBF, Ebuc )
                end select
        end function CNTBNDCalc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine CNTBNDInit ( BNDModel ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(c_int), intent(in)   :: BNDModel
        real(c_double)                  :: A, E
        integer(c_int)               :: i
        !-------------------------------------------------------------------------------------------
                CNTBNDModel= BNDModel
                CNTBNDN1   = CNTBNDN - 1.0d+00
                CNTBNDC2   = 1.0d+00 / ( CNTBNDR * CNTBNDR )
        end subroutine CNTBNDInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Module initialization
!---------------------------------------------------------------------------------------------------

        subroutine InitCNTPotModule ( STRModel, STRParams, YMType, BNDModel, Rref ) !!!!!!!!!!!!!!!!
        integer(c_int), intent(in)   :: STRModel, STRParams, YMType, BNDModel
        real(c_double), intent(in)      :: Rref
        !-------------------------------------------------------------------------------------------
                call CNTSTRInit ( STRModel, STRParams, YMType, Rref )
                call CNTBNDInit ( BNDModel )
        end subroutine InitCNTPotModule !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module CNTPot !*********************************************************************************

module TubePotBase !********************************************************************************
!
! TMD Library: Non-Bonded pair interaction potential and transfer functions for atoms composing 
! nanotubes.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
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
! All default values are adjusted for non-bonded carbob-carbon interaction in carbon nanotubes.
!
!***************************************************************************************************

use TPMLib

implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------
        
        ! Types of the potential with respect to the breathing mode
        integer*4, parameter            :: TP_POT_MONO_R        = 0
        integer*4, parameter            :: TP_POT_POLY_R        = 1

        ! Maximal number of elements in corresponding tables
        integer*4, parameter            :: TPBNMAX              = 2001

        ! Numerical constants
        real*8, parameter               :: TPbConstD            = 5.196152422706632d+00 ! = 3.0**1.5

        ! Mass of C atom
        real*8, parameter               :: TPBMc                = 12.0107d+00            ! (Da)
        
        ! Parameters of the Van der Waals inteaction between carbon atoms in graphene sheets, see
        ! Stuart S.J., Tutein A.B., Harrison J.A., J. Chem. Phys. 112(14), 2000
        real*8, parameter               :: TPBEcc               = 0.00284d+00           ! (eV)
        real*8, parameter               :: TPBScc               = 3.4d+00               ! (A)

        ! Lattice parameter and numerical density of atoms for a graphene sheet, see
        ! Dresselhaus et al, Carbon 33(7), 1995
        real*8, parameter               :: TPBAcc               = 1.421d+00             ! (A)
        real*8, parameter               :: TPBDcc               = 4.0d+00 / ( TPBConstD * TPBAcc * TPBAcc ) ! (1/A^2)
        
        ! Specific heat of carbon nanotubes
        real*8, parameter               :: TPBSHcc              = 600.0d+00 / K_MDCU    ! (eV/(Da*K))
        
        ! Cutoff distances for interactomic potential and transfer function
        ! Changes in these parameters can result in necessity to change some numerical parameters too.
        real*8, parameter               :: TPBRmincc            = 0.001d+00 * TPBScc    ! (A)
        real*8, parameter               :: TPBRcutoffcc         = 3.0d+00 * TPBScc      ! (A)
        real*8, parameter               :: TPBRcutoff1cc        = 2.16d+00 * TPBScc     ! (A)
        
        ! Parameters of the transfer function for non-bonded interaction between carbon atoms
        !real*8, parameter               :: TPBQScc             = TPBScc                ! (A)
        !real*8, parameter               :: TPBQRcutoff1cc      = 2.16d+00 * TPBScc     ! (A)
        real*8, parameter               :: TPBQScc              = 7.0d+00               ! (A)
        real*8, parameter               :: TPBQRcutoff1cc       = 8.0d+00               ! (A)
        
!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------
        
        logical                         :: TPErrCheck           = .true.                ! Set to .true. to generate diagnostic and warning messages
        character*512                   :: TPErrMsg             = ''                    ! Typically, this variable is set up in F_tt ()
        
        real*8                          :: TPGeomPrec           = 1.0d-06               ! Geometric precision, see TPInt
        integer*4                       :: TPPotType            = TP_POT_MONO_R         ! Type of the potential with respect to the breathing mode
        
        ! Physical parameters of the interatomic potential and atoms distribution at the surface
        ! of the tube
        
        real*8                          :: TPBM                 = TPBMc                 ! Mass of an atom, Da
        real*8                          :: TPBE                 = TPBEcc                ! Depth of the energy well in LJ (12-6) interatomic potential (eV)
        real*8                          :: TPBS                 = TPBScc                ! Sigma parameter of LJ (12-6) interatomic potential (A)
        real*8                          :: TPBD                 = TPBDcc                ! Numerical density of atoms at the tube surface (1/A^2)
        real*8                          :: TPBSH                = TPBSHcc               ! Specific heat (eV/(Da*K))
        
        real*8                          :: TPBRmin              = TPBRmincc             ! (A)
        real*8                          :: TPBRcutoff           = TPBRcutoffcc          ! (A)
        real*8                          :: TPBRcutoff1          = TPBRcutoff1cc         ! (A)

        ! Physical parameters of the transfer function
        
        real*8                          :: TPBQS                = TPBQScc               ! Sigma parameter of the transfer function (A)
        real*8                          :: TPBQRcutoff1         = TPBQRcutoff1cc        ! (A)
        
        ! Auxilary variables
                
        real*8                          :: TPBE4, TPBE24, TPBDRcutoff, TPBQDRcutoff
        real*8                          :: TPBQR0                                       ! Constant-value distance for the transfer function (A)
        
        ! Table of inter-particle potential, force, and transfer function
        
        integer*4                       :: TPBN                 = TPBNMAX
        real*8                          :: TPBDR
        real*8, dimension(0:TPBNMAX-1)  :: TPBQ
        real*8, dimension(0:TPBNMAX-1)  :: TPBU, TPBdUdR
        
contains !******************************************************************************************

        integer*4 function TPBsizeof () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !TPBsizeof = sizeof ( TPBU ) + sizeof ( TPBdUdR )
                TPBsizeof = 8 * ( size ( TPBQ ) + size ( TPBU ) + size ( TPBdUdR ) )
        end function TPBsizeof !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!---------------------------------------------------------------------------------------------------
! Interpolation
!---------------------------------------------------------------------------------------------------

        real*8 function TPBQInt0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: R 
        !-------------------------------------------------------------------------------------------
        real*8                  :: Z, RR
        integer*4               :: i
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

        real*8 function TPBUInt0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: R 
        !-------------------------------------------------------------------------------------------
        real*8                  :: Z, RR
        integer*4               :: i
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
        real*8, intent(out)     :: U, dUdR
        real*8, intent(in)      :: R 
        !-------------------------------------------------------------------------------------------
        real*8                  :: Z, RR
        integer*4               :: i
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
        
        real*8 function TPBQCalc0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real*8                  :: Z, t, S
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
        
        real*8 function TPBUCalc0 ( R ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: R
        !-------------------------------------------------------------------------------------------
        real*8                  :: Z, t, S
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
        real*8, intent(out)     :: U, dUdR
        real*8, intent(in)      :: R
        real*8                  :: Z, t, S, dSdR
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
        real*8, dimension(0:2), intent(out)     :: F1, F2
        real*8, dimension(0:2), intent(in)      :: F, M, Laxis
        real*8, intent(in)                      :: L
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2)                  :: FF, MM, FFF
        !-------------------------------------------------------------------------------------------
                FF = 0.5d+00 * F
                MM = M / L
                call V3_V3xxV3 ( FFF, MM, Laxis )
                F1 = FF - FFF
                F2 = FF + FFF
        end subroutine TPBSegmentForces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Printing
!---------------------------------------------------------------------------------------------------

!        subroutine TPBPrint ( FileName ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        character*(*), intent(in)       :: FileName
!        !-------------------------------------------------------------------------------------------
!       integer*4                       :: Fuid
!        integer*4                       :: i
!        real*8                          :: R
!        !-------------------------------------------------------------------------------------------
!                Fuid = OpenFile ( FileName, "wt", outputpath )
!                write ( Fuid, '(a)' ) 'TITLE="TPB Potentials"'
!                write ( Fuid, '(a)' ) 'VARIABLES="R" "Q" "U" "dUdR"'
!                write ( Fuid, '(a)' ) 'ZONE'
!                R = TPBRmin
!                do i = 0, TPBN - 1
!                        write ( Fuid, '(4e22.12)' ) R, TPBQ(i), TPBU(i), TPBDUDR(i)
!                        R = R + TPBDR
!                end do
!                call CloseFile ( Fuid )
!        end subroutine TPBPrint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Initialization
!---------------------------------------------------------------------------------------------------
        
        subroutine TPBInit () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8          :: R
        integer*4       :: i
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

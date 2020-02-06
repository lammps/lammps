module TubePotMono !********************************************************************************
!
! TMD Library: Approximate tubular potentials and transfer functions for mono-radius tubes
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
!
!---------------------------------------------------------------------------------------------------
!
! Four potentials and transfer functions are calculated in this module:
!
! 1. SSTP (segment - semi-infinite tube parallel). It gives a linear density of the potential along 
!    the segment axis which produced by a parallel semi-infinite tube. 2D tables for this potential 
!    are generated at initialization or can be loaded from a file
!
! 2. STP (segment - tube parallel). It gives a linear density of the potential along the segment axis 
!    which produced by a parallel infinite tubes. This is only a particular case of the SSTP potential, 
!    but it is considered separately for computational effiency. 1D tables of this potential are taken 
!    from 2D tables of SSTP potential.
!
! 3. SST (segment - semi-infinite tube). It gives a potential for a segment produced by a arbitrary-
!    oriented semi-infinite tube. Data of this potential can not be kept in 2D tabels, therefore all 
!    data are calculated 'on fly' with the help of SSTP potential and numerical integration along the 
!    segment axis
!
! 4. ST (segment - tube). It gives a potential for a segment produced by a arbitrary-oriented 
!    infinitely long tube. 2D tables for this potential are generated at initialization or can be 
!    loaded from a file
!
!***************************************************************************************************

use TPMLib
use TPMGeom
use TubePotBase
use TubePotTrue
use LinFun2
use Spline2

implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------

        integer*4, parameter                    :: TPMNZMAX     = 129
        integer*4, parameter                    :: TPMNEMAX     = 128
        
        integer*4, parameter                    :: TPMNHMAX     = 1001  
        integer*4, parameter                    :: TPMNXMAX     = 1001 
        integer*4, parameter                    :: TPMNMAX      = 1001 
        
!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------
        
        integer*4                               :: TPMStartMode = 1
        character*512                           :: TPMSSTPFile  = 'TPMSSTP.xrs'
        character*512                           :: TPMAFile     = 'TPMA.xrs'
        
        integer*4                               :: TPMNZ        = TPMNZMAX
        integer*4                               :: TPMNZ1       = TPMNZMAX - 1
        integer*4                               :: TPMNE        = TPMNEMAX
        integer*4                               :: TPMNE1       = TPMNEMAX - 1
        
        integer*4                               :: TPMNH        = TPMNHMAX
        integer*4                               :: TPMNH1       = TPMNHMAX - 1
        integer*4                               :: TPMNX        = TPMNXMAX
        integer*4                               :: TPMNX1       = TPMNXMAX - 1
        
        integer                                 :: TPMChiIndM   ! Chirality index M 
        integer                                 :: TPMChiIndN   ! Chirality index N
        real*8                                  :: TPMR1
        real*8                                  :: TPMR2
        
        real*8                                  :: TPMHmax
        real*8                                  :: TPMDH
        
        ! Parameters of empirical correction functions
        
        integer*4                               :: TPMAN        = 20
        real*8                                  :: TPMAHmin
        real*8                                  :: TPMAHmax
        real*8                                  :: TPMADH
        real*8, dimension(0:TPMNHMAX-1)         :: TPMAH, TPMAF, TPMAFxx
        
        ! Fitting parameters that depend on the SWCNT chirality

        real*8                                  :: TPMCaA       = 0.22d+00 ! 0.22 for (10,10) CNTs
        real*8                                  :: TPMCeA       = 0.35d+00 ! 0.35 for (10,10) CNTs
        real*8                                  :: TPMAHmin0    = 10.0d+00 ! 10.0 A for (10,10) CNTs             
        
        ! Parameters of SSTP integrator
        
        real*8                                  :: TPMDE
        real*8, dimension(0:TPMNEMAX-1)         :: TPMCE, TPMSE
        
        ! Additional parameters for SSTP potential

        real*8                                  :: TPMSSTPDelta = 0.25d+00
        integer*4                               :: TPMSSTPNH
        integer*4                               :: TPMSSTPNX
        
        real*8                                  :: TPMSSTPX1
        real*8                                  :: TPMSSTPXmax
        real*8                                  :: TPMSSTPDX
        
        real*8, dimension(0:TPMNHMAX-1,0:TPMNXMAX-1) :: TPMSSTPG
        real*8, dimension(0:TPMNHMAX-1,0:TPMNXMAX-1) :: TPMSSTPF, TPMSSTPFxx, TPMSSTPFyy, TPMSSTPFxxyy
        real*8, dimension(0:TPMNHMAX-1)         :: TPMSSTPH
        real*8, dimension(0:TPMNXMAX-1)         :: TPMSSTPX
        
        ! Additional parameters for STP potential
        
        ! In calcuation of this potential also some parameters of SSTP potential are used
        ! In particular, STP potential has no its own integrator. All data comes from SSTP integrator.
        ! It does not result in any computational inefficiency unless the STP potential is used without SSTP one.
        
        integer*4                               :: TPMNN        = 10
        real*8, dimension(0:TPMNHMAX-1)         :: TPMSTPG
        real*8, dimension(0:TPMNHMAX-1)         :: TPMSTPF, TPMSTPFxx

        ! Parameters for ST potential

        real*8                                  :: TPMSTDelta   = 1.0d+00       ! Minimal gap dh for ST-potential
        integer*4                               :: TPMSTNXS     = 10            ! Number of subdivisions for every grid step in ST-integrator
        real*8                                  :: TPMSTXmax
        real*8                                  :: TPMSTH1
        real*8                                  :: TPMSTH2
        real*8                                  :: TPMSTDH12
        
        real*8, dimension(0:TPMNHMAX-1,0:TPMNXMAX-1) :: TPMSTG
        real*8, dimension(0:TPMNHMAX-1,0:TPMNXMAX-1) :: TPMSTF, TPMSTFxx, TPMSTFyy, TPMSTFxxyy
        real*8, dimension(0:TPMNHMAX-1) :: TPMSTH
        real*8, dimension(0:TPMNXMAX-1) :: TPMSTX
        
        ! Switching parameters
        
        ! Height switch (at H=0 in SST-potential)
        integer*4                               :: TPMHSwitch   = 0             ! 1, use h-switch; 0, do not use the switch
        real*8                                  :: TPMHS        = 3.0d+00       ! Switch height, Angstrom
        
        ! Angle switch
        integer*4                               :: TPMASwitch   = 0             ! 1, use a-switch; 0, do not use the switch
        real*8                                  :: TPMAS        = 3.0d+00       ! Switch angle, degree
        real*8                                  :: TPMASMin
        real*8                                  :: TPMASMax
        real*8                                  :: TPMASDelta

        ! These variables are used to print error message if intertube force filed fails
        integer*4                               :: Err_CNT1 = 0, Err_CNT1_Node = 0, Err_CNT2 = 0, Err_CNT2_Node1 = 0, Err_CNT2_Node2 = 0, Err_EType = 0

contains !******************************************************************************************

        integer*4 function TPMsizeof () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                TPMsizeof = 8 * ( size ( TPMAH ) + size ( TPMAF ) + size ( TPMAFxx ) &
                        + size ( TPMCE ) + size ( TPMSE ) + size ( TPMSSTPG ) + size ( TPMSSTPF ) &
                        + size ( TPMSSTPFxx ) + size ( TPMSSTPFyy ) + size ( TPMSSTPFxxyy ) &
                        + size ( TPMSSTPH ) + size ( TPMSSTPX ) + size ( TPMSTPG ) + size ( TPMSTPF ) &
                        + size ( TPMSTPFxx ) + size ( TPMSTG ) + size ( TPMSTF ) + size ( TPMSTFxx ) &
                        + size ( TPMSTFyy ) + size ( TPMSTFxxyy ) + size ( TPMSTH ) + size ( TPMSTX ) )
        end function TPMsizeof !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Printing error message if intertube force field fails
!---------------------------------------------------------------------------------------------------

        subroutine PrintTPErrMsg () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !write ( TPErrMsg, fmt = '(a,i8,a,i8,a,i8,a,i8,a,i8,a,i1)' ) 'CNT ', Err_CNT1, ' [', Err_CNT1_Node,'] with CNT ', Err_CNT2, ' [', Err_CNT2_Node1, ', ', Err_CNT2_Node2, '] E=', Err_EType
                !call PrintStdLogMsg ( TPErrMsg )
        end subroutine PrintTPErrMsg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! SSTP: Linear potential density for the tube interacting with parallel semi-infinte tube
!---------------------------------------------------------------------------------------------------

        subroutine TPMSSTPIntegrator ( Q, U, H, D ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function calculates the transfer function Q and potential U between an infinitely long
        ! tube and a cross-section of another parallel tube for given height H and displacemnet D.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U
        real*8, intent(in)      :: H, D
        !-------------------------------------------------------------------------------------------
        integer*4               :: i, j, k
        real*8                  :: C, Zmin, Zmax, DZ, R1X, R1Y, R2X, R2Y, R2Z, R, Rcutoff2
        !-------------------------------------------------------------------------------------------
                Q = 0.0d+00
                U = 0.0d+00
                Zmin = D - TPBRcutoff
                Zmax = D + TPBRcutoff
                Rcutoff2 = TPBRcutoff * TPBRcutoff
                if ( Zmin < 0.0d+00 ) Zmin = 0.0d+00
                DZ = ( Zmax - Zmin ) / TPMNZ1
                do i = 0, TPMNE1 ! Integration over the section of the first tube
                        R1X = TPMR1 * TPMCE(i)
                        R1Y = TPMR1 * TPMSE(i)
                        do j = 0, TPMNE1 ! !Integration over the section of the second tube
                                R2X = TPMR1 * TPMCE(j) + H
                                R2Y = TPMR1 * TPMSE(j)
                                R2Z = Zmin - D
                                do k = 0, TPMNZ1 ! Integration over the length of the second tube
                                        R = sqr ( R2X - R1X ) + sqr ( R2Y - R1Y ) + sqr ( R2Z ) 
                                        if ( R < Rcutoff2 ) then
                                                R = dsqrt ( R )
                                                if ( k == 0 .or. k == TPMNZ1 ) then
                                                        Q = Q + 0.5d+00 * TPBQCalc0 ( R )
                                                        U = U + 0.5d+00 * TPBUCalc0 ( R )
                                                else
                                                        Q = Q + TPBQCalc0 ( R )
                                                        U = U + TPBUCalc0 ( R )
                                                end if
                                        end if
                                        R2Z = R2Z + DZ
                                end do
                        end do
                end do
                C = sqr ( TPMDE ) * TPMR1 * TPMR2 * DZ
                Q = Q * C
                U = U * sqr ( TPBD ) * C
        end subroutine TPMSSTPIntegrator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMSSTPInt0 ( Q, U, H, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the transfer function Q and potential U for the SSTP potential
        ! calculated with interpolation in the table without switch
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U
        real*8, intent(in)      :: H, X
        !-------------------------------------------------------------------------------------------
        integer*4               :: i, j
        real*8                  :: XX
        !-------------------------------------------------------------------------------------------
                i = 1 + int ( H / TPMDH )
                j = 1 + int ( ( X + TPMSSTPXMax ) / TPMSSTPDX )
                if ( ( i < TPMSSTPNH ) .and. ( j > TPMSSTPNX ) ) then
                        !call PrintTPErrMsg ()        
                        !!call PrintStdLogMsg (TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10)' ) 'Tubes intersect each other: H=', H, ', X=', X
                        !call Error ( 'TPMSSTPInt0', TPErrMsg )
                end if
                if ( i > TPMNH1 ) then
                        Q = 0.0d+00
                        U = 0.0d+00
                        TPMSSTPInt0 = 0
                        return
                end if
                if ( X .le. - TPMSSTPXmax ) then
                        j  = 1
                        XX = - TPMSSTPXmax
                else if ( X .ge. TPMSSTPXmax ) then
                        j  = TPMNX1
                        XX = TPMSSTPXmax
                else
                        XX = X
                end if
                Q = CalcLinFun2_0 ( i, j, H, XX, TPMNH, TPMNX, TPMSSTPH, TPMSSTPX, TPMSSTPG ) 
                U = CalcSpline2_0 ( i, j, H, XX, TPMNH, TPMNX, TPMSSTPH, TPMSSTPX, TPMSSTPF, TPMSSTPFxx, TPMSSTPFyy, TPMSSTPFxxyy ) 
                TPMSSTPInt0 = 1
        end function TPMSSTPInt0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSSTPInt0S ( Q, U, H, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the transfer function Q and potential U for the SSTP potential
        ! calculated with interpolation in the table and switch to the case of zero H
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U
        real*8, intent(in)      :: H, X
        !-------------------------------------------------------------------------------------------
        integer*4               :: IntSign
        real*8                  :: t, W, Qa, Ua
        !-------------------------------------------------------------------------------------------
                if ( TPMHSwitch == 0 ) then
                        TPMSSTPInt0S = TPMSSTPInt0 ( Q, U, H, X )
                else
                        if ( H > TPMHS ) then
                                TPMSSTPInt0S = TPMSSTPInt0 ( Q, U, H, X )
                        else
                                t = H / TPMHS
                                W = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                TPMSSTPInt0S = TPMSSTPInt0 ( Q, U, H, X )
                                IntSign = TPMSSTPInt0 ( Qa, Ua, 0.0d+00, X )
                                Q = W * Qa + ( 1.0d+00 - W ) * Q
                                U = W * Ua + ( 1.0d+00 - W ) * U
                        end if
                end if
        end function TPMSSTPInt0S !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMSSTPInt1 ( Q, U, Uh, Ux, H, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the transfer function Q, potential U, and derivarives Uh=dU/dH and
        ! Ux=dU/dX for the SSTP potential calculated with interpolation in the table without switch
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U, Uh, Ux
        real*8, intent(in)      :: H, X
        !-------------------------------------------------------------------------------------------
        integer*4               :: i, j
        real*8                  :: XX
        !-------------------------------------------------------------------------------------------
                i = 1 + int ( H / TPMDH )
                j = 1 + int ( ( X + TPMSSTPXMax ) / TPMSSTPDX )
                if ( ( i < TPMSSTPNH ) .and. ( j > TPMSSTPNX ) ) then
                        !call PrintTPErrMsg ()        
                        !!call PrintStdLogMsg ( TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10)' ) 'Tubes intersect each other: H=', H, ', X=', X
                        !call Error ( 'TPMSSTPInt1', TPErrMsg )
                end if
                if ( i > TPMNH1 ) then
                        Q = 0.0d+00
                        U = 0.0d+00
                        Uh = 0.0d+00
                        Ux = 0.0d+00
                        TPMSSTPInt1 = 0
                        return
                end if
                if ( X .le. - TPMSSTPXmax ) then
                        j  = 1
                        XX = - TPMSSTPXmax
                else if ( X .ge. TPMSSTPXmax ) then
                        j  = TPMNX1
                        XX = TPMSSTPXmax
                else
                        XX = X
                end if
                Q = CalcLinFun2_0 ( i, j, H, XX, TPMNH, TPMNX, TPMSSTPH, TPMSSTPX, TPMSSTPG ) 
                call CalcSpline2_1 ( U, Uh, Ux, i, j, H, XX, TPMNH, TPMNX, TPMSSTPH, TPMSSTPX, TPMSSTPF, TPMSSTPFxx, TPMSSTPFyy, TPMSSTPFxxyy ) 
                TPMSSTPInt1 = 1
        end function TPMSSTPInt1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMSSTPInt1S ( Q, U, Uh, Ux, H, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the transfer function Q, potential U, and derivarives Uh=dU/dH and
        ! Ux=dU/dX for the SSTP potential calculated with interpolation in the table and switch to 
        ! the case of zero H.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U, Uh, Ux
        real*8, intent(in)      :: H, X
        !-------------------------------------------------------------------------------------------
        integer*4               :: IntSign
        real*8                  :: t, W, W1, dWdH, Qa, Ua, Uha, Uxa
        !-------------------------------------------------------------------------------------------
                if ( TPMHSwitch == 0 ) then
                        TPMSSTPInt1S = TPMSSTPInt1 ( Q, U, Uh, Ux, H, X )
                else
                        if ( H > TPMHS ) then
                                TPMSSTPInt1S = TPMSSTPInt1 ( Q, U, Uh, Ux, H, X )
                        else
                                t = H / TPMHS
                                W = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                dWdH = 6.0d+00 * t * ( t - 1.0d+00 ) / TPMHS
                                TPMSSTPInt1S = TPMSSTPInt1 ( Q, U, Uh, Ux, H, X )
                                IntSign = TPMSSTPInt1 ( Qa, Ua, Uha, Uxa, 0.0d+00, X )
                                W1 = 1.0d+00 - W
                                Q = W * Qa + W1 * Q
                                U = W * Ua + W1 * U
                                Uh = W1 * Uh + ( Ua - U ) * dWdH
                                Ux = W * Uxa + W1 * Ux
                        end if
                end if
        end function TPMSSTPInt1S !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine TPMSSTPWrite () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function writes the table of the SSTP potential to the disk file
        !-------------------------------------------------------------------------------------------
        integer*4 :: Fuid, i, j
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( TPMSSTPFile, 'wt', '' )
                write ( unit = Fuid, fmt = '(4i8)' ) TPMChiIndM, TPMChiIndN, TPMNH1, TPMNX1
                do i = 0, TPMNH1
                        do j = 0, TPMNX1
                                if ( ( i .ge. TPMSSTPNH ) .or. ( j .le. TPMSSTPNX ) ) write ( unit = Fuid, fmt = '(2e26.17)' ) TPMSSTPG(i,j), TPMSSTPF(i,j)
                        end do
                end do
                call CloseFile ( Fuid )
        end subroutine TPMSSTPWrite !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPMSSTPRead () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function reads the table of the SSTP potential from the disk file
        !-------------------------------------------------------------------------------------------
        integer*4       :: Fuid, i, j
        integer*4       :: iTPMChiIndM, iTPMChiIndN, iTPMNH1, iTPMNX1
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( TPMSSTPFile, 'rt', '' )
                read ( unit = Fuid, fmt = '(4i8)' ) iTPMChiIndM, iTPMChiIndN, iTPMNH1, iTPMNX1
                if ( iTPMChiIndM .NE. TPMChiIndM .OR. iTPMChiIndN .NE. TPMChiIndN ) then
                        print *, 'ERROR in [TPMSSTPRead]: iTPMChiIndM .NE. TPMChiIndM .OR. iTPMChiIndN .NE. TPMChiIndN'
                        stop
                end if
                if ( iTPMNH1 .NE. TPMNH1 .OR. iTPMNX1 .NE. TPMNX1 ) then
                        print *, 'ERROR in [TPMSSTPRead]: iTPMNH1 .NE. TPMNH1 .OR. iTPMNX1 .NE. TPMNX1'
                        stop
                end if
                do i = 0, TPMNH1
                        do j = 0, TPMNX1
                                if ( ( i .ge. TPMSSTPNH ) .or. ( j .le. TPMSSTPNX ) ) read ( unit = Fuid, fmt = '(2e26.17)' ) TPMSSTPG(i,j), TPMSSTPF(i,j)
                        end do
                end do
                call CloseFile ( Fuid )
        end subroutine TPMSSTPRead !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine TPMSSTPInit () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function calculates the table of the SSTP potential
        !-------------------------------------------------------------------------------------------
        integer*4       :: i, j
        real*8          :: E
        character*512   :: Msg
        real*8, dimension(0:TPMNMAX-1) :: FF, DD, MM, K0, K1, K2
        !-------------------------------------------------------------------------------------------
                TPMDE = M_2PI / TPMNE
                E = 0.0d+00
                do i = 0, TPMNE1
                        TPMCE(i) = cos ( E ) 
                        TPMSE(i) = sin ( E )
                        E = E + TPMDE
                end do
                do i = 0, TPMNH1
                        TPMSSTPH(i) = TPMDH * i
                end do
                TPMSSTPX1 = - 2.0d+00 * TPMSSTPDelta
                TPMSSTPXmax = TPBRcutoff + TPMSSTPDelta
                TPMSSTPDX = 2.0 * TPMSSTPXmax / TPMNX1
                do j = 0, TPMNX1
                        TPMSSTPX(j) = - TPMSSTPXmax + TPMSSTPDX * j
                end do
                TPMSSTPNH = 1 + int ( ( TPMR1 + TPMR2 + TPMSSTPDelta ) / TPMDH )
                TPMSSTPNX = int ( ( TPMSSTPXMax + TPMSSTPX1 ) / TPMSSTPDX ) - 1
                if ( TPMStartMode == 0 ) then
                        do i = 0, TPMNH1
                                do j = 0, TPMNX1
                                        if ( ( i .ge. TPMSSTPNH ) .or. ( j .le. TPMSSTPNX ) ) then
                                                call TPMSSTPIntegrator ( TPMSSTPG(i,j), TPMSSTPF(i,j), TPMSSTPH(i), TPMSSTPX(j) )
                                                print '(2i5,a,e20.10,a,e20.10,a,e20.10,a,e20.10)', i, j, ' H=', TPMSSTPH(i), ', X=', TPMSSTPX(j), ', Q=', TPMSSTPG(i,j), ', U=', TPMSSTPF(i,j)
                                        end if
                                end do
                        end do
                        call TPMSSTPWrite ()
                else
                        call TPMSSTPRead ()
                end if
                call CreateSpline2Ext ( 3, 3, 3, 3, TPMNH, TPMSSTPNH, TPMNX, TPMSSTPNX, TPMNMAX, TPMSSTPH, TPMSSTPX, TPMSSTPF, TPMSSTPFxx, TPMSSTPFyy, TPMSSTPFxxyy, FF, MM, DD, K0, K1, K2 )                
        end subroutine TPMSSTPInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! STP Potential for an infinite tube interacting with a parallel segment. No actual initialization 
! is necessary for this potential, since the data are taken from the table for SSTP potenrials.
!---------------------------------------------------------------------------------------------------

        integer*4 function TPMSTPInt0 ( Q, U, H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the transfer function Q and potential U for the STP potential
        ! calculated with interpolation in the table
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
        integer*4               :: i
        !-------------------------------------------------------------------------------------------
                i = 1 + int ( H / TPMDH )
                if ( i < TPMSSTPNH ) then
                        !call PrintTPErrMsg ()        
                        !!call PrintStdLogMsg ( TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10)' ) 'Tubes intersect each other: H=', H
                        !call Error ( 'TPMSTPInt0', TPErrMsg )
                end if
                if ( H > TPMHmax ) then
                        Q = 0.0d+00
                        U = 0.0d+00
                        TPMSTPInt0 = 0
                        return
                end if
                if ( i == TPMNH ) i = TPMNH - 1
                Q = CalcLinFun1_0 ( i, H, TPMNH, TPMSSTPH, TPMSTPG ) 
                U = CalcSpline1_0 ( i, H, TPMNH, TPMSSTPH, TPMSTPF, TPMSTPFxx ) 
                TPMSTPInt0 = 1
        end function TPMSTPInt0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSTPInt1 ( Q, U, dUdH, H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the transfer function Q, potential U, and derivative dUdH for 
        ! the STP potential calculated with interpolation in the table
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U, dUdH
        real*8, intent(in)      :: H
        integer*4               :: i
        !-------------------------------------------------------------------------------------------
                i = 1 + int ( H / TPMDH )
                if ( i < TPMSSTPNH ) then
                        !call PrintTPErrMsg ()        
                        !!call PrintStdLogMsg ( TPErrMsg )
                        !write ( TPErrMsg, '(a,e20.10)' ) 'Tubes intersect each other: H=', H
                        !call Error ( 'TPMSTPInt0', TPErrMsg )
                end if
                if ( H > TPMHmax ) then
                        Q = 0.0d+00
                        U = 0.0d+00
                        dUdH = 0.0d+00
                        TPMSTPInt1 = 0
                        return
                end if
                Q = CalcLinFun1_0 ( i, H, TPMNH, TPMSSTPH, TPMSTPG ) 
                call CalcSpline1_1 ( U, dUdH, i, H, TPMNH, TPMSSTPH, TPMSTPF, TPMSTPFxx ) 
                TPMSTPInt1 = 1
        end function TPMSTPInt1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine TPMSTPInit () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function initializes the table of the STP potential
        !-------------------------------------------------------------------------------------------
                TPMSTPG(0:TPMNH1)   = TPMSSTPG(0:TPMNH1,TPMNX1)
                TPMSTPF(0:TPMNH1)   = TPMSSTPF(0:TPMNH1,TPMNX1)
                TPMSTPFxx(0:TPMNH1) = TPMSSTPFxx(0:TPMNH1,TPMNX1)
        end subroutine TPMSTPInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!---------------------------------------------------------------------------------------------------
! Fitting functions for SST and ST potential.
! This correction functions are choosen empirically to improve accuracy of SST and ST potentials.
!---------------------------------------------------------------------------------------------------

        subroutine TPMAInit ( X1_1, X1_2, X2_1, X2_2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)              :: X1_1, X1_2, X2_1, X2_2
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2)          :: R1_1, R1_2, R2_1, R2_2
        real*8, dimension(0:2)          :: Fa, Ma
        real*8                          :: Qa, Ua, Qb, Ub, X, H, HH, Ucoeff, Uamin, Ubmin
        integer*4                       :: i, j, IntSign, Fuid
        real*8, dimension(0:TPMNHMAX-1) :: D, K0, K1, K2
        integer*4                       :: iTPMChiIndM, iTPMChiIndN, iTPMAN
        !-------------------------------------------------------------------------------------------
                TPMAHmin = TPMR1 + TPMR2 + TPMSTDelta
                TPMAHmax = TPMR1 + TPMR2 + 0.95d+00 * TPBRcutoff
                TPMADH = ( TPMAHmax - TPMAHmin ) / ( TPMAN - 1 )
                if ( TPMStartMode == 1 ) then
                        Fuid = OpenFile ( TPMAFile, 'rt', '' )
                        read ( unit = Fuid, fmt = '(4i8)' ) iTPMChiIndM, iTPMChiIndN, iTPMAN
                        if ( iTPMChiIndM .NE. TPMChiIndM .OR. iTPMChiIndN .NE. TPMChiIndN ) then
                                print *, 'ERROR in [TPMAInit]: iTPMChiIndM .NE. TPMChiIndM .OR. iTPMChiIndN .NE. TPMChiIndN'
                                stop
                        end if
                        if ( iTPMAN .NE. TPMAN ) then
                                print *, 'ERROR in [TPMAInit]: iTPMAN .NE. TPMAN'
                                stop
                        end if
                        do i = 0, TPMAN - 1
                                TPMAH(i) = TPMAHmin + i * TPMADH
                                read ( unit = Fuid, fmt = * ) TPMAF(i)
                        end do
                        call CloseFile ( Fuid )
                        call CreateSpline1 ( 3, 3, TPMAN, TPMAH, TPMAF, TPMAFxx, D, K0, K1, K2 )
                        return
                end if
                call TPTInit ( TPMR1, TPMR2, TPTNXMAX, TPTNEMAX )
                do i = 0, TPMAN - 1
                        TPMAH(i) = TPMAHmin + i * TPMADH
                        call TPTGetEnds ( R1_1, R1_2, R2_1, R2_2, X1_1, X1_2, X2_1, X2_2, TPMAH(i), M_PI_2 )
                        X = - ( X1_2 - X1_1 ) / 2.0d+00
                        do j = 0, ( TPTNXMAX - 1 ) / 2
                                HH = sqrt ( TPMAH(i) * TPMAH(i) + sqr ( X * sin ( M_PI_2 ) ) )
                                IntSign = TPMSTPInt0 ( Qb, Ub, HH ) 
                                call TPTSetSegPosition2 ( TPTSeg1, R1_1, R1_2 )
                                call TPTSetSegPosition2 ( TPTSeg2, R2_1, R2_2 )
                                IntSign = TPTSectionPotential ( Qa, Ua, Fa, Ma, TPTSeg1, j, TPTSeg2 ) 
                                if ( j == 0 ) then
                                        Uamin = Ua
                                        Ubmin = Ub
                                else 
                                        if ( Ua < Uamin ) then
                                                Uamin = Ua
                                        end if
                                        if ( Ub < Ubmin ) then
                                                Ubmin = Ub
                                        end if
                                end if
                                X = X + TPTSeg1%DX
                        end do
                        TPMAF(i) = Uamin / Ubmin
                end do
                Fuid = OpenFile ( TPMAFile, 'wt', '' )
                write ( unit = Fuid, fmt = '(4i8)' ) TPMChiIndM, TPMChiIndN, TPMAN
                do i = 0, TPMAN - 1
                        write ( unit = Fuid, fmt = * ) TPMAF(i)
                end do
                call CloseFile ( Fuid )
                call CreateSpline1 ( 3, 3, TPMAN, TPMAH, TPMAF, TPMAFxx, D, K0, K1, K2 )
        end subroutine TPMAInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real*8 function TPMA0 ( H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
        integer*4               :: i
        real*8                  :: A0, t, S
        !-------------------------------------------------------------------------------------------
                if ( H > TPMAHmax ) then
                        TPMA0 = 1.0d+00
                        return
                else if ( H < TPMAHmin ) then
                        if ( H < TPMAHmin0 ) then
                                TPMA0 = 1.5d+00
                                return
                        end if
                        A0   = 1.5d+00
                        t    = ( H - TPMAHmin0 ) / TPMAHmin
                        S    = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                        TPMA0 = ( 1.0d+00 - S ) * CalcSpline1_0 ( 1, H, TPMAN, TPMAH, TPMAF, TPMAFxx ) + A0 * S
                        return
                end if
                i = 1 + int ( ( H - TPMAHmin ) / TPMADH )
                TPMA0 = CalcSpline1_0 ( i, H, TPMAN, TPMAH, TPMAF, TPMAFxx ) 
        end function TPMA0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPMA1 ( A, Ah, H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)     :: A, Ah
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
        integer*4               :: i
        real*8                  :: A0, t, S, dSdH
        !-------------------------------------------------------------------------------------------
                if ( H > TPMAHmax ) then
                        A = 1.0d+00
                        Ah = 0.0d+00
                        return
                else if ( H < TPMAHmin ) then
                        if ( H < TPMAHmin0 ) then
                                A = 1.5d+00
                                Ah = 0.0d+00
                                return
                        end if
                        A0   = 1.5d+00
                        t    = ( H - TPMAHmin0 ) / TPMAHmin
                        S    = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                        dSdH = 6.0d+00 * t * ( t - 1.0d+00 ) / TPMAHmin
                        call CalcSpline1_1 ( A, Ah, 1, H, TPMAN, TPMAH, TPMAF, TPMAFxx )
                        Ah = Ah * ( 1.0d+00 - S ) + dSdH * ( A0 - A )
                        A  = A * ( 1.0d+00 - S ) + A0 * S
                        return
                end if
                i = 1 + int ( ( H - TPMAHmin ) / TPMADH )
                call CalcSpline1_1 ( A, Ah, i, H, TPMAN, TPMAH, TPMAF, TPMAFxx )
        end subroutine TPMA1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        real*8 function TPMCu0 ( H, cosA, sinA ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the correction function for the magnitude of the potential. 
        !-------------------------------------------------------------------------------------------
        real*8, intent(in)      :: H, cosA, sinA
        !-------------------------------------------------------------------------------------------
                TPMCu0 = 1.0d+00 + ( TPMA0 ( H ) - 1.0d+00 ) * sqr ( sinA )
        end function TPMCu0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPMCu1 ( Cu, CuH, CuA, H, cosA, sinA ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Thi subroutine calculates the correction function Cu for magnitude of the potential and 
        ! its derivatives CuH, CuA.
        !-------------------------------------------------------------------------------------------
        real*8, intent(ouT)     :: Cu, CuH, CuA
        real*8, intent(in)      :: H, cosA, sinA
        real*8                  :: AA, AAh, D
        !-------------------------------------------------------------------------------------------
                call TPMA1 ( AA, AAh, H )
                D = sqr ( sinA ) 
                AA = AA - 1.0d+00
                Cu = 1.0d+00 + AA * D
                CuH = AAh * D
                CuA = AA * 2.0d+0 * cosA * sinA
        end subroutine TPMCu1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real*8 function TPMCa0 ( cosA, sinA ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the correction function for the argument of the potential. 
        ! If correction is not necessary, it should return sinA.
        !-------------------------------------------------------------------------------------------
        real*8, intent(in)      :: cosA, sinA
        !-------------------------------------------------------------------------------------------
                TPMCa0 = sinA / ( 1.0d+00 - TPMCaA * sqr ( sinA ) )
        end function TPMCa0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPMCa1 ( Ca, CaA, Ka, KaA, cosA, sinA ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This subroutine calculates the correction function Cu for the depth of the potential well 
        ! and its derivatives CuH, CuA. If correction is not necessary, it should return Ca = sinA 
        ! and CaA = cosA.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Ca, CaA, Ka, KaA
        real*8, intent(in)      :: cosA, sinA
        !-------------------------------------------------------------------------------------------
                Ka  = 1.0d+00 / ( 1.0d+00 - TPMCaA * sqr ( sinA ) )
                Ca = sinA * Ka 
                KaA = 2.0d+00 * TPMCaA * sinA * cosA * sqr ( Ka )
                CaA = cosA * Ka + sinA * KaA
        end subroutine TPMCa1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real*8 function TPMCe0 ( sinA ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the correction function for the argument of the potential. 
        ! If correction is not necessary, it should return sinA.
        !-------------------------------------------------------------------------------------------
        real*8, intent(in)      :: sinA
        !-------------------------------------------------------------------------------------------
                TPMCe0 = 1.0d+00 - TPMCeA * sinA * sinA
        end function TPMCe0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPMCe1 ( Ce, CeA, Ke, cosA, sinA ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! If correction is not necessary, it should return Ce = 1 and CeA = 0.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Ce, CeA, Ke
        real*8, intent(in)      :: cosA, sinA
        !-------------------------------------------------------------------------------------------
                Ce = 1.0d+00 - TPMCeA * sinA * sinA
                CeA = - 2.0d+00 * TPMCeA * sinA * cosA
                Ke = TPMCeA
        end subroutine TPMCe1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!---------------------------------------------------------------------------------------------------
! SST Potential for the semi-infinite tube interacting with segment. 
! This potential does not need any initialization. All necessry data is taken from tables of the 
! SSTP potential.
!---------------------------------------------------------------------------------------------------

        integer*4 function TPMSSTPotential ( Q, U, X1, X2, H, cosA, D, N ) !!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function calculates the transfer function Q and potenial U applyed to a segment 
        ! from asemi-infinte tube based on numerical integration (trapesond rule) along the segment 
        ! axis for non-parallel objects. 
        ! Relative position of the nanotube and segment are given by axial positions of the segment
        ! ends X1 and X2, height H, cosA= cos(A), where A is the cross-axis angle, and displacement 
        ! D of the nanotube end.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)     :: Q, U
        real*8, intent(in)      :: X1, X2, H, cosA, D
        integer*4, intent(in)   :: N ! Number of nodes for numerical integration
        real*8                  :: sinA, Qs, Us, DX, X, XX, HH, Cu, Ca, Ce
        integer*4               :: i
        !-------------------------------------------------------------------------------------------
                Q = 0.0d+00
                U = 0.0d+00
                DX = ( X2 - X1 ) / ( N - 1 )
                X = X1
                sinA = dsqrt ( 1.0d+00 - cosA * cosA )
                Cu = TPMCu0 ( H, cosA, sinA )
                Ca = TPMCa0 ( cosA, sinA )
                Ce = TPMCe0 ( sinA )
                TPMSSTPotential = 0
                do i = 0, N - 1
                        XX = X * cosA - Ce * D
                        HH = sqrt ( H * H + sqr ( X * Ca ) )
                        if ( TPMSSTPInt0S ( Qs, Us, HH, XX ) > 0 ) TPMSSTPotential = 1
                        if ( i == 0 .or. i == N - 1 ) then
                                Q = Q + 0.5d+00 * Qs
                                U = U + 0.5d+00 * Us
                        else
                                Q = Q + Qs
                                U = U + Us
                        end if
                        X = X + DX
                end do
                Q = Cu * Q * DX
                U = Cu * U * DX
        end function TPMSSTPotential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSSTPotentialPar ( Q, U, R1_1, Laxis1, R2_1, Laxis2, L1, N ) !!!!!!!!!! 
        ! Potential applyed to the segment from the semi-infinte tube is calculated by numerical 
        ! integration (trapesond rule) along the segment axis for parallel objects.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)                     :: Q, U
        real*8, dimension(0:2), intent(in)      :: R1_1, Laxis1, R2_1, Laxis2
        real*8, intent(in)                      :: L1
        integer*4, intent(in)                   :: N ! Number of nodes for numerical integration
        !-------------------------------------------------------------------------------------------
        real*8                                  :: Qs, Us, DX, X, S, H
        real*8, dimension(0:2)                  :: R1, L12
        integer*4                               :: i
        !-------------------------------------------------------------------------------------------
                DX = L1 / ( N - 1 )
                X = 0.0d+00
                Q = 0.0d+00
                U = 0.0d+00
                TPMSSTPotentialPar = 0
                do i = 0, N - 1
                        R1 = R1_1 + X * Laxis1
                        call LinePoint ( S, L12, R2_1, Laxis2, R1 )
                        L12 = L12 - R1
                        call ApplyPeriodicBC ( L12 )
                        H = S_V3norm3 ( L12 )
                        if ( TPMSSTPInt0S ( Qs, Us,  H, S ) > 0 ) then
                                TPMSSTPotentialPar = 1
                                if ( i == 0 .or. i == TPMNN - 1 ) then
                                        Q = Q + 0.5d+00 * Qs
                                        U = U + 0.5d+00 * Us
                                else
                                        Q = Q + Qs
                                        U = U + Us
                                end if
                                X = X + DX
                        end if
                end do
                Q = Q * DX
                U = U * DX
        end function TPMSSTPotentialPar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMSSTForces ( Q, U, F1, F2, Fd, X1, X2, H, cosA, D, N ) !!!!!!!!!!!!!!!!
        ! Potential and forces applyed to the segment from the semi-infinte  tube are calculated 
        ! by numerical integration (trapesond rule) along the segment axis.
        ! Non-parallel case.
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)                     :: Q, U, Fd
        real*8, dimension(0:2), intent(out)     :: F1, F2
        real*8, intent(in)                      :: X1, X2, H, cosA, D
        integer*4, intent(in)                   :: N ! Number of nodes for numerical integration
        !-------------------------------------------------------------------------------------------
        real*8                                  :: DX, sinA
        real*8                                  :: Qs, Us, Ush, Usx, Fx, Fy, Fz
        real*8                                  :: C, C1, C2, I0, Ih, Ih1, Ih2, Ix, Ix1, X, XX, HH
        real*8                                  :: Ca, CaA, Ka, KaA, Cu, CuH, CuA, Ce, CeA, Ke, Uh, Ua
        integer*4                               :: IntSign, i
        !-------------------------------------------------------------------------------------------
                I0   = 0.0d+00
                Ih   = 0.0d+00
                Ih1  = 0.0d+00
                Ih2  = 0.0d+00
                Ix   = 0.0d+00
                Ix1  = 0.0d+00
                Q    = 0.0d+00
                U    = 0.0d+00
                F1   = 0.0d+00
                F2   = 0.0d+00
                Fd   = 0.0d+00
                sinA = dsqrt ( 1.0d+00 - cosA * cosA )
                X    = X1
                DX   = ( X2 - X1 ) / ( N - 1 )
                TPMSSTForces = 0
                call TPMCa1 ( Ca, CaA, Ka, KaA, cosA, sinA )
                call TPMCu1 ( Cu, CuH, CuA, H, cosA, sinA )
                call TPMCe1 ( Ce, CeA, Ke, cosA, sinA )
                do i = 0, N - 1
                        XX = X * cosA - Ce * D
                        HH = sqrt ( H * H + sqr ( X * Ca ) )
                        if ( TPMSSTPInt1S ( Qs, Us, Ush, Usx, HH, XX ) > 0 ) TPMSSTForces = 1
                        if ( i == 0 .or. i == N - 1 ) then
                                Qs  = 0.5d+0 * Qs
                                Us  = 0.5d+0 * Us
                                Ush = 0.5d+0 * Ush
                                Usx = 0.5d+0 * Usx
                        end if
                        if ( HH .le. TPGeomPrec ) then
                                Ush = 0.0d+00
                        else
                                Ush = Ush / HH
                        end if
                        Q   = Q + Qs
                        I0  = I0 + Us
                        Ih  = Ih + Ush
                        Ih1 = Ih1 + X * Ush
                        Ih2 = Ih2 + X * X * Ush
                        Ix  = Ix + Usx
                        Ix1 = Ix1 + X * Usx
                        X   = X + DX
                end do
                if ( TPMSSTForces == 0 ) return
                
                C   = DX * Cu
                I0  = I0 * C
                Ih  = Ih * C
                Ih1 = Ih1 * C
                Ih2 = Ih2 * C
                Ix  = Ix * C
                Ix1 = Ix1 * C

                DX  = X2 - X1 

                Q  = Q * C
                U  = I0
                Uh = ( CuH * U / Cu + h * Ih ) / DX
                Ua = ( CuA * I0 / Cu + Ca * CaA * Ih2 - sinA * Ix1 - CeA * D * Ix ) / DX

                C1 = Ka * Ka * Ih1
                C  = h * ( C1 + cosA * Ke * Ix ) / DX
                F1(0) = X2 * Uh - C 
                F2(0) = C - X1 * Uh
                
                C = ( cosA * sinA * C1 + ( Ke * sinA - sinA ) * Ix ) / DX
                F1(1) = Ua - X2 * C 
                F2(1) = X1 * C - Ua
                
                C1 = Ca * Ca * Ih1 + cosA * Ix
                C2 = Ca * Ca * Ih2 + cosA * Ix1
                F1(2) = ( U - X2 * C1 + C2 ) / DX
                F2(2) = ( X1 * C1 - C2 - U ) / DX
                
                Fd = Ce * Ix
        end function TPMSSTForces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSSTForcesPar ( Q, U, F1, F2, Fd, R1_1, Laxis1, R2_1, Laxis2, L1, N ) !
        ! Potential and forces applyed to the segment from the semi-infinte tube are calculated by 
        ! numerical integration (trapesond rule) along the segment axis.
        ! Non-parallel case
        !-------------------------------------------------------------------------------------------
        real*8, intent(out)                     :: Q, U, Fd
        real*8, dimension(0:2), intent(out)     :: F1, F2
        real*8, dimension(0:2), intent(in)      :: R1_1, Laxis1, R2_1, Laxis2
        real*8, intent(in)                      :: L1
        integer*4, intent(in)                   :: N ! Number of nodes for numerical integration
        !-------------------------------------------------------------------------------------------
        real*8                                  :: Qs, Us, Ush, Usx, DX, X, S, H, Beta, Gamma
        real*8, dimension(0:2)                  :: R1, L12, Fs
        integer*4                               :: i, N1
        !-------------------------------------------------------------------------------------------
                Q  = 0.0d+00
                U  = 0.0d+00
                F1 = 0.0d+00
                F2 = 0.0d+00
                Fd = 0.0d+00
                X  = 0.0d+00
                N1 = N - 1
                DX = L1 / N1
                TPMSSTForcesPar = 0
                do i = 0, N1
                        R1 = R1_1 + X * Laxis1
                        call LinePoint ( S, L12, R2_1, Laxis2, R1 )
                        L12 = L12 - R1
                        call ApplyPeriodicBC ( L12 )
                        H = S_V3norm3 ( L12 )
                        if ( TPMSSTPInt1S ( Qs, Us, Ush, Usx, H, S ) > 0 ) then
                                TPMSSTForcesPar = 1
                                if ( H .ge. TPGeomPrec ) then
                                        Fs = Ush * L12 / H -  Usx * Laxis2
                                else
                                        Fs = - Usx * Laxis2
                                end if
                                Beta = real ( i ) / N1
                                Gamma = 1.0d+00 - Beta
                                if ( i == 0 .or. i == N1 ) then
                                        Q = Q + 0.5d+00 * Qs
                                        U = U + 0.5d+00 * Us
                                        Fd = Fd + 0.5d+00 * Usx
                                        Gamma = 0.5d+00 * Gamma
                                        Beta = 0.5d+000 * Beta
                                else
                                        Q  = Q + Qs
                                        U  = U + Us
                                        Fd = Fd + Usx
                                end if
                                F1 = F1 + Gamma * Fs 
                                F2 = F2 + Beta * Fs 
                        end if
                        X = X + DX
                end do
                Q = Q * DX
                U = U * DX
                Fd = Fd * DX
                Fs = U * Laxis1 / L1
                F1 = DX * F1 + Fs
                F2 = DX * F2 - Fs
        end function TPMSSTForcesPar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!---------------------------------------------------------------------------------------------------
! ST: Potential for the infinite tube interacting with segment
!--------------------------------------------------------------------------------------------------
        
        !
        ! These functions are used to smooth boundaries in (H,X) domain for ST potential
        !
        
        real*8 function TPMSTXMin0 ( H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
        real*8                  :: X
        !-------------------------------------------------------------------------------------------
                if ( H < TPMSTH1 ) then
                        TPMSTXMin0 = sqrt ( TPMSTH2 * TPMSTH2 - H * H )
                        return
                else if ( H > TPMSTH2 ) then
                        TPMSTXMin0 = 0.0d+00
                        return
                end if
                X = ( H - TPMSTH1 ) / TPMSTDH12
                TPMSTXMin0 = sqrt ( TPMSTH2 * TPMSTH2 - H * H ) * ( 1.0d+00 - X * X * X * ( 3.0d+00 * X * ( 2.0d+00 * X - 5.0d+00 ) + 10.0d+00 ) )
        end function TPMSTXMin0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        real*8 function TPMSTXMax0 ( H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
                TPMSTXMax0 = sqrt ( TPMSTXMax * TPMSTXMax - H * H )
        end function TPMSTXMax0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine TPMSTXMin1 ( XMin, dXMindH, H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)     :: XMin, dXMindH
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
        real*8                  :: X, F, dFdX
        !-------------------------------------------------------------------------------------------
                if ( H < TPMSTH1 ) then
                        XMin = sqrt ( TPMSTH2 * TPMSTH2 - H * H )
                        dXMindH = - H / XMin
                        return
                else if ( H > TPMSTH2 ) then
                        XMin     = 0.0d+00
                        dXMindH = 0.0d+00
                        return
                end if
                X = ( H - TPMSTH1 ) / TPMSTDH12
                F = 1.0d+00 - X * X * X * ( 3.0d+00 * X * ( 2.0d+00 * X - 5.0d+00 ) + 10.0d+00 )
                X = X * ( X - 1.0d+00 )
                dFdX = - 30.0d+00 * X * X
                XMin = sqrt ( TPMSTH2 * TPMSTH2 - H * H ) 
                dXMindH = - H * F / XMin + XMin * dFdX / TPMSTDH12
                XMin = XMin * F
        end subroutine TPMSTXMin1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine TPMSTXMax1 ( XMax, dXMaxdH, H ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)     :: XMax, dXMaxdH
        real*8, intent(in)      :: H
        !-------------------------------------------------------------------------------------------
                XMax    = sqrt ( TPMSTXMax * TPMSTXMax - H * H )
                dXMaxdH = - H / XMax
        end subroutine TPMSTXMax1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !
        ! ST (segment-tube) table
        !
        
        subroutine TPMSTIntegrator ( G, F, Q, U, H, X, DX ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(inout)   :: G, F, Q, U
        real*8, intent(in)      :: H, X, DX
        !-------------------------------------------------------------------------------------------
        real*8                  :: FFx, HH, DDX
        integer*4               :: IntSign
        !-------------------------------------------------------------------------------------------
                DDX = 0.5 * DX
                G = G + Q * DDX
                F = F + U * DDX
                Q = 0.0d+00
                U = 0.0d+00
                HH  = dsqrt ( sqr ( H ) + sqr ( X ) )
                if ( HH > TPMHmax ) return
                IntSign = TPMSTPInt0 ( Q, U, HH ) 
                if ( IntSign == 1 ) then
                        G = G + Q * DDX
                        F = F + U * DDX
                end if
        end subroutine TPMSTIntegrator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSTInt0 ( G, F, H, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)     :: G, F
        real*8, intent(in)      :: H, X
        !-------------------------------------------------------------------------------------------
        integer*4               :: i, j
        real*8                  :: S, XA, XXX, XXXX, XMin, XMax
        !-------------------------------------------------------------------------------------------
                if ( H > TPMHmax ) then
                        G = 0.0d+00
                        F = 0.0d+00
                        TPMSTInt0 = 0
                        return
                else if ( H < 0.0d+00 ) then
                        G = 0.0d+00
                        F = 0.0d+00
                        TPMSTInt0 = 2
                        !call PrintTPErrMsg ()        
                        !!call PrintStdLogMsg ( TPErrMsg )
                        !all Error ( 'TPMSTInt0', 'H < 0' )
                        !!return 
                end if
                S = signum ( X )
                XA = dabs ( X )
                i = 1 + int ( H / TPMDH )
                if ( i > TPMNH1 ) i = TPMNH1
                XMin = TPMSTXMin0 ( H )
                XMax = TPMSTXMax0 ( H )
                XXX  = ( XA - XMin ) / ( XMax - XMin )
                if ( XXX < 0.0d+00 ) then
                        j  = 1
                        XXXX = 0.0d+00
                        !call PrintTPErrMsg ()        
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10,a,e20.10)' ) 'X < XMin, X=', XA, ', XMin=', XMin, ', H=', H
                        !call Error ( 'TPMSTInt0', TPErrMsg )
                else if ( XXX .ge. 1.0d+00 ) then
                        j  = TPMNX1
                        XXXX = 1.0d+00
                else
                        XXXX = XXX
                        j = 1 + int ( XXXX * TPMNX1 )
                end if
                G = S * CalcLinFun2_0 ( i, j, H, XXXX, TPMNH, TPMNX, TPMSTH, TPMSTX, TPMSTG ) 
                F = S * CalcSpline2_0 ( i, j, H, XXXX, TPMNH, TPMNX, TPMSTH, TPMSTX, TPMSTF, TPMSTFxx, TPMSTFyy, TPMSTFxxyy ) 
                TPMSTInt0 = 1
        end function TPMSTInt0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSTInt1 ( G, F, Fh, Fx, H, X ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(inout)   :: G, F, Fh, Fx
        real*8, intent(in)      :: H, X
        !-------------------------------------------------------------------------------------------
        integer*4               :: i, j
        real*8                  :: S, XA, DX, XXX, XXXX, XMin, XMax, dXMindH, dXMaxdH
        !-------------------------------------------------------------------------------------------
                if ( H > TPMHmax ) then
                        G  = 0.0d+00
                        F  = 0.0d+00
                        Fh = 0.0d+00
                        Fx = 0.0d+00
                        TPMSTInt1 = 0
                        return
                else if ( H < 0.0d+00 ) then
                        G  = 0.0d+00
                        F  = 0.0d+00
                        Fh = 0.0d+00
                        Fx = 0.0d+00
                        TPMSTInt1 = 2
                        !call PrintTPErrMsg ()        
                        !!call PrintStdLogMsg ( trim ( TPErrMsg ) )
                        !call Error ( 'TPMSTInt1', 'H < 0' )
                        !!return 
                end if
                S = signum ( X )
                XA = dabs ( X )
                i = 1 + int ( H / TPMDH )
                if ( i > TPMNH1 ) i = TPMNH1
                call TPMSTXMin1 ( XMin, dXMindH, H )
                call TPMSTXMax1 ( XMax, dXMaxdH, H )
                DX = XMax - XMin
                XXX  = ( XA - XMin ) / DX
                if ( XXX < 0.0d+00 ) then
                        j  = 1
                        XXX = 0.0d+00
                        XXXX = 0.0d+00
                        !call PrintTPErrMsg ()        
                        !write ( TPErrMsg, '(a,e20.10,a,e20.10,a,e20.10)' ) 'X < XMin, X=', XA, ', XMin=', XMin, ', H=', H
                        !call Error ( 'TPMSTInt', TPErrMsg )
                else if ( XXX .ge. 1.0d+00 ) then
                        j  = TPMNX1
                        XXX = 1.0d+00
                        XXXX = 1.0d+00
                else
                        XXXX = XXX 
                        j = 1 + int ( XXXX * TPMNX1 )
                end if
                G = S * CalcLinFun2_0 ( i, j, H, XXXX, TPMNH, TPMNX, TPMSTH, TPMSTX, TPMSTG ) 
                call CalcSpline2_1 ( F, Fh, Fx, i, j, H, XXXX, TPMNH, TPMNX, TPMSTH, TPMSTX, TPMSTF, TPMSTFxx, TPMSTFyy, TPMSTFxxyy ) 
                Fx = Fx / DX
                Fh = Fh - Fx * ( dXMaxdH * XXX + dXMindH * ( 1.0d+00 - XXX ) )
                F = F * S
                Fh = Fh * S
                TPMSTInt1 = 1
        end function TPMSTInt1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMSTPotential ( Q, U, X1, X2, H, cosA, CaseID ) !!!!!!!!!!!!!!!!!!!!!!!!
        real*8, intent(out)     :: Q, U
        real*8, intent(in)      :: X1, X2, H, cosA
        integer*4, intent(in)   :: CaseID
        !-------------------------------------------------------------------------------------------
        real*8                  :: sinA, GG1, GG2, FF1, FF2, Ca, Cu
        !-------------------------------------------------------------------------------------------
                if ( CaseID == MD_LINES_PAR ) then
                        TPMSTPotential = TPMSTPInt0 ( Q, U, H )
                        U = U * ( X2 - X1 )
                        return 
                end if
                TPMSTPotential = 0
                sinA = dsqrt ( 1.0d+00 - cosA * cosA )
                Cu = TPMCu0 ( H, cosA, sinA )
                Ca = TPMCa0 ( cosA, sinA )
                if ( TPMSTInt0  ( GG1, FF1, H, X1 * Ca ) > 0 ) TPMSTPotential = 1
                if ( TPMSTInt0  ( GG2, FF2, H, X2 * Ca ) > 0 ) TPMSTPotential = 1
                Q = Cu * ( GG2 - GG1 ) / Ca
                U = Cu * ( FF2 - FF1 ) / Ca
        end function TPMSTPotential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSTForces ( Q, U, F1, F2, X1, X2, H, cosA, CaseID ) !!!!!!!!!!!!!!!!!!!
        real*8, intent(out)                     :: Q, U
        real*8, dimension(0:2), intent(out)     :: F1, F2
        real*8, intent(in)                      :: X1, X2, H, cosA
        integer*4, intent(in)                   :: CaseID
        !-------------------------------------------------------------------------------------------
        real*8                                  :: DX, sinA
        real*8                                  :: GG1, GG2, FF1, FF2, Fh1, Fh2, Fx1, Fx2
        real*8                                  :: B, C, D
        real*8                                  :: Ca, CaA, Ka, KaA, Cu, CuH, CuA
        integer*4                               :: IntSign1, IntSign2
        !-------------------------------------------------------------------------------------------
                DX = X2 - X1
                if ( CaseID == MD_LINES_PAR ) then
                        F1 = 0.0d+00
                        F2 = 0.0d+00
                        TPMSTForces = TPMSTPInt1 ( Q, U, F1(0), H )
                        F1(0) = F1(0) * 0.5 * DX
                        F2(0) = F1(0)
                        F1(2) = U
                        F2(2) = - U
                        Q     = Q * DX
                        U     = U * DX
                        return 
                end if

                sinA = dsqrt ( 1.0d+00 - cosA * cosA )
                call TPMCa1 ( Ca, CaA, Ka, KaA, cosA, sinA )
                IntSign1 = TPMSTInt1  ( GG1, FF1, Fh1, Fx1, H, X1 * Ca )               
                IntSign2 = TPMSTInt1  ( GG2, FF2, Fh2, Fx2, H, X2 * Ca )
                if ( ( IntSign1 .ne. 1 ) .and. ( IntSign2 .ne. 1 ) ) then
                        Q = 0.0d+00
                        U = 0.0d+00
                        F1 = 0.0d+00
                        F2 = 0.0d+00
                        TPMSTForces = 0
                        return 
                end if
                
                call TPMCu1 ( Cu, CuH, CuA, H, cosA, sinA )
                
                Q = Cu * ( GG2 - GG1 ) / Ca
                U = Cu * ( FF2 - FF1 ) / Ca

                B = Cu * ( Fx2 - Fx1 ) / sinA
                C = H * B / sinA
                D = CuH * U / Cu + Cu * ( Fh2 - Fh1 ) / Ca
                F1(0) = ( X2 * D - C ) / DX
                F2(0) = ( C - X1 * D ) / DX
                
                C = cosA * B
                D = ( CuA / Cu - CaA / Ca ) * U + CaA * Cu * ( X2 * Fx2 - X1 * Fx1 ) / Ca
                F1(1) = ( D - X2 * C ) / DX
                F2(1) = ( X1 * C - D ) / DX
                
                F1(2) = Cu * Fx1
                F2(2) = - Cu * Fx2
                
                TPMSTForces = 1
        end function TPMSTForces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMSTForceTorque( Qi, Ui, Fi, Ti, Q, U, F, T, Psi, PsiA, Cap, L, H, cosA, CaseID ) 
        real*8, intent(out)                     :: Qi, Ui, Fi, Ti, Q, U, F, T, Psi, PsiA, Cap
        real*8, intent(in)                      :: L, H, cosA
        integer*4, intent(in)                   :: CaseID
        !-------------------------------------------------------------------------------------------
        real*8                                  :: L2, sinA
        real*8                                  :: GG, FF, Fh, Fx, GGi, FFi, Fhi, Fxi
        real*8                                  :: B, C, D
        real*8                                  :: Ca, CaA, Ka, KaA, Cu, CuH, CuA
        integer*4                               :: IntSign
        !-------------------------------------------------------------------------------------------
                if ( CaseID == MD_LINES_PAR ) then
                        TPMSTForceTorque = TPMSTPInt1 ( Q, U, F, H )
                        Q  = Q * L
                        U  = U * L
                        F  = - F * L
                        T  = 0.0d+00
                        Qi = 0.0d+00
                        Ui = 0.0d+00
                        Fi = 0.0d+00
                        Ti = 0.0d+00
                        Psi = 0.0d+00
                        PsiA = 0.0d+00
                        Cap = 0.0d+00
                        return 
                end if
                
                L2 = 0.5d+00 * L
                sinA = dsqrt ( 1.0d+00 - cosA * cosA )
                call TPMCa1 ( Ca, CaA, Ka, KaA, cosA, sinA )
                IntSign = TPMSTInt1  ( GG, FF, Fh, Fx, H, L2 * Ca )               
                IntSign = TPMSTInt1  ( GGi, FFi, Fhi, Fxi, H, TPMSTXmax )               
                if ( IntSign .ne. 1 ) then
                        Qi = 0.0d+00
                        Ui = 0.0d+00
                        Fi = 0.0d+00
                        Ti = 0.0d+00
                        Q  = 0.0d+00
                        U  = 0.0d+00
                        F  = 0.0d+00
                        T  = 0.0d+00
                        Psi  = 0.0d+00
                        PsiA = 0.0d+00
                        Cap  = 0.0d+00
                        TPMSTForceTorque = 0
                        return 
                end if
                
                call TPMCu1 ( Cu, CuH, CuA, H, cosA, sinA )
                
                Psi = Cu / Ca
                PsiA = ( CuA * Ca - Cu * CaA ) / Ca / Ca
                Cap = CuA / Cu - KaA / Ka - cosA / sinA
                Qi = 2.0d+00 * Psi * GGi
                Ui = 2.0d+00 * Psi * FFi
                Fi = - 2.0d+00 * ( CuH * FFi / Ca + Psi * Fhi )
                Ti = - Cap * Ui
                
                Q = 2.0d+00 * Cu * GG / Ca
                U = 2.0d+00 * Cu * FF / Ca
                F = - 2.0d+00 * ( CuH * FF / Ca + Psi * Fh )
                T = - 2.0d+00 * ( ( CuA * Ka - Cu * KaA ) / ( Ka * Ka * sinA ) - Cu * cosA / ( Ka * sinA * sinA ) ) * FF &
                    - 2.0d+00 * Cu / ( Ka * sinA ) * Fx * L2 * ( KaA * sinA + Ka * cosA )
                
                TPMSTForceTorque = 1
        end function TPMSTForceTorque !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine TPMSTInit () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8                          :: X, Q, U, DX, DDX, XMin, XMax
        integer*4                       :: i, j, k
        real*8, dimension(0:TPMNMAX-1)  :: FF, DD, MM, K0, K1, K2
        !-------------------------------------------------------------------------------------------
                TPMSTH1   = TPMR1 + TPMR2
                TPMSTH2   = TPMSTH1 + TPMSTDelta
                TPMSTDH12 = TPMSTH2 - TPMSTH1
                TPMSTXmax = TPMHMax + TPMSTDelta
                DX = 1.0 / TPMNX1
                do j = 0, TPMNX1  
                        TPMSTX(j) = DX * j
                end do
                do i = 0, TPMNH1
                        TPMSTH(i) = TPMDH * i
                        XMin = TPMSTXmin0 ( TPMSTH(i) )
                        XMax = TPMSTXMax0 ( TPMSTH(i) )
                        Q  = 0.0d+00
                        U  = 0.0d+00
                        DX = ( XMax - XMin ) * TPMSTX(1) / TPMSTNXS
                        X  = XMin
                        call TPMSTIntegrator ( TPMSTG(i,0), TPMSTF(i,0), Q, U, TPMSTH(i), X, DX )
                        TPMSTG(i,0) = 0.0d+00
                        TPMSTF(i,0) = 0.0d+00
                        TPMSTFyy(i,0) = U
                        TPMSTFyy(i,TPMNX1) = 0.0d+00
                        do j = 1, TPMNX1  
                                TPMSTG(i,j) = TPMSTG(i,j-1)
                                TPMSTF(i,j) = TPMSTF(i,j-1)
                                do k = 0, TPMSTNXS - 1
                                        X = X + DX
                                        call TPMSTIntegrator ( TPMSTG(i,j), TPMSTF(i,j), Q, U, TPMSTH(i), X, DX )
                                end do
                                if ( j < TPMNX1 ) DX  = ( XMax - XMin ) * ( TPMSTX(j+1) - TPMSTX(j) ) / TPMSTNXS
                        end do
                end do
                call CreateSpline2 ( 3, 3, 3, 3, TPMNH, TPMNX, TPMNMAX, TPMSTH, TPMSTX, TPMSTF, TPMSTFxx, TPMSTFyy, TPMSTFxxyy, FF, MM, DD, K0, K1, K2 )                
        end subroutine TPMSTInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Interaction functions: They can be used for calculation of the potential and forces between a 
! segment and infinte or semi-infinite nanotube.
!---------------------------------------------------------------------------------------------------

        subroutine TPMSegmentForces ( F2_1, F2_2, F1_1, F1_2, R1_1, R1_2, R2, Laxis2, L2 ) !!!!!!!!!
        real*8, dimension(0:2), intent(out)     :: F2_1, F2_2
        real*8, dimension(0:2), intent(in)      :: F1_1, F1_2, R1_1, R1_2, R2, Laxis2
        real*8, intent(in)                      :: L2 
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2)                  :: F, M, RR
        !-------------------------------------------------------------------------------------------
                RR = R1_1 - R2
                ! Taking into account periodic boundaries
                call ApplyPeriodicBC ( RR )
                call V3_V3xxV3 ( M, RR, F1_1 )
                RR = R1_2 - R2
                ! Taking into account periodic boundaries
                call ApplyPeriodicBC ( RR )
                call V3_V3xxV3 ( F, RR, F1_2 )
                M = - ( M + F )
                F = - ( F1_1 + F1_2 )
                call TPBSegmentForces ( F2_1, F2_2, F, M, Laxis2, L2 )
        end subroutine TPMSegmentForces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Interaction of a segment with semi-infinite or infinite tube
        !

        integer*4 function TPMInteractionF ( Q, U, F1_1, F1_2, F2_1, F2_2, Fd, R1_1, R1_2, R2_1, R2_2, SType2 )
        ! SType2 in the type of the second segment:
        !   SType2 == 0, internal segment
        !   Stype2 == 1, point R2_1 is the end of the tube
        !   Stype2 == 2, point R2_2 in the end of the tube
        !-------------------------------------------------------------------------------------------
        real*8, intent(inout)                   :: Q, U, Fd
        real*8, dimension(0:2), intent(inout)   :: F1_1, F1_2, F2_1, F2_2
        real*8, dimension(0:2), intent(in)      :: R1_1, R1_2, R2_1, R2_2
        !-------------------------------------------------------------------------------------------
        integer*4                               :: SType2
        real*8, dimension(0:2)                  :: R1, R2, Laxis1, Laxis2, F1, F2, L12, Ly, DR, F1_1a, F1_2a, F1_1b, F1_2b
        real*8                                  :: H, cosA, D1, D2, L1, L2, cosA2, t, W, W1, dWdt, Qa, Ua, Qb, Ub, Fda, Fdb, FF
        integer*4                               :: GeomID, SwitchID, S, IntSigna, IntSignb
        !-------------------------------------------------------------------------------------------
                R1 = 0.5d+00 * ( R1_1 + R1_2 )
                R2 = 0.5d+00 * ( R2_1 + R2_2 )
                Laxis1 = R1_2 - R1_1
                Laxis2 = R2_2 - R2_1
                L1 = S_V3norm3 ( Laxis1 )
                L2 = S_V3norm3 ( Laxis2 )
                Laxis1 = Laxis1 / L1
                Laxis2 = Laxis2 / L2
                L1 = 0.5d+00 * L1 
                L2 = 0.5d+00 * L2
                if ( SType2 == 2 ) Laxis2 = - Laxis2
                GeomID = LineLine ( H, cosA, D1, D2, L12, R1, Laxis1, R2, Laxis2, TPGeomPrec )
                
                ! Angle switch
                if ( TPMASwitch == 0 ) then
                        if ( GeomID == MD_LINES_PAR ) then
                                SwitchID = 2
                        else
                                SwitchID = 0
                        end if
                else
                        cosA2 = cosA * cosA
                        if ( cosA2 .ge. TPMASMax .or. GeomID == MD_LINES_PAR ) then
                                SwitchID = 2
                        else if ( cosA2 .le. TPMASMin ) then
                                SwitchID = 0
                        else
                                t = ( cosA2 - TPMASMin ) / TPMASDelta
                                W = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                dWdt = 6.0d+00 * t * ( t - 1.0d+00 ) / TPMASDelta
                                SwitchID = 1
                        end if
                end if
                
                if ( SwitchID < 2 ) then
                        D2 = D2 - L2
                        if ( SType2 == 0 ) then
                                IntSigna = TPMSTForces  ( Qa, Ua, F1, F2, D1 - L1, D1 + L1, H, cosA, MD_LINES_NONPAR )
                                Fda = 0.0d+00
                        else
                                IntSigna = TPMSSTForces  ( Qa, Ua, F1, F2, Fda, D1 - L1, D1 + L1, H, cosA, D2, TPMNN )
                        end if
                        call V3_V3xxV3 ( Ly, Laxis1, Laxis2 )
                        S = signum ( S_V3xV3 ( Ly, L12 ) )
                        call V3_V3xxV3 ( Ly, Laxis1, L12 )
                        Ly = Ly * S
                        if ( IntSigna > 0 ) then
                                F1_1a = F1(0) * L12 + F1(1) * Ly + F1(2) * Laxis1
                                F1_2a = F2(0) * L12 + F2(1) * Ly + F2(2) * Laxis1 
                        else
                                F1_1a = 0.0d+00
                                F1_2a = 0.0d+00
                        end if
                end if
                
                if ( SwitchID > 0 ) then 
                        if ( SType2 == 0 ) then
                                call LinePoint ( H, L12, R2, Laxis2, R1 )
                                L12 = L12 - R1
                                call ApplyPeriodicBC ( L12 )
                                H = S_V3norm3 ( L12 )
                                IntSignb = TPMSTForces  ( Qb, Ub, F1, F2, - L1, L1, H, cosA, MD_LINES_PAR )
                                Fdb = 0.0d+00
                                if ( IntSignb > 0 ) then
                                        if ( H .le. TPGeomPrec ) then
                                                F1_1b = F1(2) * Laxis1
                                                F1_2b = F2(2) * Laxis1
                                        else
                                                L12 = L12 / H
                                                F1_1b = F1(0) * L12 + F1(2) * Laxis1
                                                F1_2b = F2(0) * L12 + F2(2) * Laxis1 
                                        end if
                                else
                                        F1_1b = 0.0d+00
                                        F1_2b = 0.0d+00
                                end if
                        else if ( Stype2 == 1 ) then
                                IntSignb = TPMSSTForcesPar ( Qb, Ub, F1_1b, F1_2b, Fdb, R1_1, Laxis1, R2_1, Laxis2, 2.0d+00 * L1, TPMNN )
                        else
                                IntSignb = TPMSSTForcesPar ( Qb, Ub, F1_1b, F1_2b, Fdb, R1_1, Laxis1, R2_2, Laxis2, 2.0d+00 * L1, TPMNN )
                        end if
                end if

                if ( SwitchID == 0 ) then
                        Q    = Qa
                        U    = Ua
                        F1_1 = F1_1a
                        F1_2 = F1_2a
                        Fd   = Fda
                        TPMInteractionF = IntSigna
                else if ( SwitchID == 2 ) then
                        Q    = Qb
                        U    = Ub
                        F1_1 = F1_1b
                        F1_2 = F1_2b
                        Fd   = Fdb
                        TPMInteractionF = IntSignb
                else 
                        W1   = 1.0d+00 - W
                        Q    = W * Qa + W1 * Qb
                        U    = W * Ua + W1 * Ub
                        Ly   = Ly * ( Ua - Ub ) * dWdt * cosA * sqrt ( 1.0d+00 - sqr ( cosA ) ) / L1
                        F1_1 = W * F1_1a + W1 * F1_1b - Ly
                        F1_2 = W * F1_2a + W1 * F1_2b + Ly
                        Fd   = W * Fda + W1 * Fdb
                        TPMInteractionF = 0
                        if ( IntSigna > 0 .or. IntSignb > 0 ) TPMInteractionF = 1
                end if      

                ! Calculation of forces for the comlimentary tube
                if ( SType2 == 2 ) Laxis2 = - Laxis2
                call TPMSegmentForces ( F2_1, F2_2, F1_1, F1_2, R1_1, R1_2, R2, Laxis2, 2.0d+00 * L2 )
                ! After the previous subroutine F2_1*Laxis2 = F2_2*Laxis2, but this is not true for the semi-infinite tube.
                ! The force along the tube sould be applied to the end of the tube, while for the
                ! another point corresponding force is equal to zero.
                if ( SType2 == 1 ) then
                        FF = S_V3xV3 ( F2_1, Laxis2 )
                        DR = ( Fd - FF ) * Laxis2
                        F2_1 = F2_1 + DR
                        F2_2 = F2_2 - DR
                else if ( SType2 == 2 ) then
                        FF = S_V3xV3 ( F2_2, Laxis2 )
                        DR = ( - Fd - FF ) * Laxis2
                        F2_2 = F2_2 + DR
                        F2_1 = F2_1 - DR
                end if
        end function TPMInteractionF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer*4 function TPMInteractionU ( Q, U, R1_1, R1_2, R2_1, R2_2, SType2 ) !!!!!!!!!!!!!!!!
        real*8, intent(inout)                   :: Q, U
        real*8, dimension(0:2), intent(in)      :: R1_1, R1_2, R2_1, R2_2
        integer*4, intent(in)                   :: SType2
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2)                  :: R1, R2, Laxis1, Laxis2, F1, F2, L12, DR
        real*8                                  :: H, cosA, D1, D2, L1, L2, cosA2, t, W, Qa, Ua, Qb, Ub
        integer*4                               :: GeomID, SwitchID, IntSigna, IntSignb
        !-------------------------------------------------------------------------------------------
                R1 = 0.5d+00 * ( R1_1 + R1_2 )
                R2 = 0.5d+00 * ( R2_1 + R2_2 )
                Laxis1 = R1_2 - R1_1
                Laxis2 = R2_2 - R2_1
                L1 = S_V3norm3 ( Laxis1 )
                L2 = S_V3norm3 ( Laxis2 )
                Laxis1 = Laxis1 / L1
                Laxis2 = Laxis2 / L2
                if ( SType2 == 2 ) Laxis2 = - Laxis2
                GeomID = LineLine ( H, cosA, D1, D2, L12, R1, Laxis1, R2, Laxis2, TPGeomPrec )
                L1 = 0.5d+00 * L1
                L2 = 0.5d+00 * L2

                ! Angle switch
                if ( TPMASwitch == 0 ) then
                        if ( GeomID == MD_LINES_PAR ) then
                                SwitchID = 2
                        else
                                SwitchID = 0
                        end if
                else
                        cosA2 = cosA * cosA
                        if ( cosA2 .ge. TPMASMax .or. GeomID == MD_LINES_PAR ) then
                                SwitchID = 2
                        else if ( cosA2 .le. TPMASMin ) then
                                SwitchID = 0
                        else
                                t = ( cosA2 - TPMASMin ) / TPMASDelta
                                W = 1.0d+00 - t * t * ( 3.0d+00 - 2.0d+00 * t )
                                SwitchID = 1
                        end if
                end if

                if ( SwitchID < 2 ) then
                        if ( Stype2 == 0 ) then
                                IntSigna = TPMSTPotential  ( Qa, Ua, D1 - L1, D1 + L1, H, cosA, MD_LINES_NONPAR )
                        else
                                IntSigna = TPMSSTPotential  ( Qa, Ua, D1 - L1, D1 + L1, H, cosA, D2 - L2, TPMNN )
                        end if
                end if
                
                if ( SwitchID > 0 ) then
                        if ( Stype2 == 0 ) then
                                call LinePoint ( H, L12, R2, Laxis2, R1 )
                                L12 = L12 - R1
                                call ApplyPeriodicBC ( L12 )
                                IntSignb = TPMSTPotential  ( Qb, Ub, - L1, L1, S_V3norm3 ( L12 ), cosA, MD_LINES_PAR )
                        else if ( Stype2 == 1 ) then
                                IntSignb = TPMSSTPotentialPar ( Qb, Ub, R1_1, Laxis1, R2_1, Laxis2, 2.0d+00 * L1, TPMNN )
                        else
                                IntSignb = TPMSSTPotentialPar ( Qb, Ub, R1_1, Laxis1, R2_2, Laxis2, 2.0d+00 * L1, TPMNN )
                        end if
                end if
                
                if ( SwitchID == 0 ) then
                        Q = Qa
                        U = Ua
                        TPMInteractionU = IntSigna
                else if ( SwitchID == 2 ) then
                        Q = Qb
                        U = Ub
                        TPMInteractionU = IntSignb
                else
                        Q = W * Qa + ( 1.0d+00 - W ) * Qb
                        U = W * Ua + ( 1.0d+00 - W ) * Ub
                        TPMInteractionU = 0
                        if ( IntSigna > 0 .or. IntSignb > 0 ) TPMInteractionU = 1
                end if
        end function TPMInteractionU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer*4 function TPMInteractionFNum ( Q, U, F1_1, F1_2, F2_1, F2_2, R1_1, R1_2, R2_1, R2_2, Stype2, Delta )
        real*8, intent(inout)                   :: Q, U
        real*8, dimension(0:2), intent(inout)   :: F1_1, F1_2, F2_1, F2_2
        real*8, dimension(0:2), intent(in)      :: R1_1, R1_2, R2_1, R2_2
        integer*4, intent(in)                   :: SType2
        real*8, intent(in)                      :: Delta
        !-------------------------------------------------------------------------------------------
        integer*4                               :: i, j, IntSign
        real*8                                  :: QQ, DD, D2
        real*8, dimension(0:1,0:2)              :: U1_1, U1_2, U2_1, U2_2
        real*8, dimension(0:2)                  :: RR
        !-------------------------------------------------------------------------------------------
                U = 0.0d+00
                F1_1 = 0.0d+00
                F1_2 = 0.0d+00
                F2_1 = 0.0d+00
                F2_2 = 0.0d+00
                TPMInteractionFNum = TPMInteractionU ( Q, U, R1_1, R1_2, R2_1, R2_2, SType2 )
                if ( TPMInteractionFNum == 0 ) return
                D2 = 2.0d+00 * Delta
                do i = 0, 2
                        DD = - Delta
                        do j = 0 , 1  
                                RR = R1_1 
                                RR(i) = RR(i) + DD
                                IntSign = TPMInteractionU ( QQ, U1_1(j,i), RR, R1_2, R2_1, R2_2, SType2 )
                                RR = R1_2 
                                RR(i) = RR(i) + DD
                                IntSign = TPMInteractionU ( QQ, U1_2(j,i), R1_1, RR, R2_1, R2_2, SType2 )
                                RR = R2_1 
                                RR(i) = RR(i) + DD;
                                IntSign = TPMInteractionU ( QQ, U2_1(j,i), R1_1, R1_2, RR, R2_2, SType2 )
                                RR = R2_2
                                RR(i) = RR(i) + DD
                                IntSign = TPMInteractionU ( QQ, U2_2(j,i), R1_1, R1_2, R2_1, RR, SType2 )
                                DD = DD + D2
                        end do
                end do
                do i = 0, 2
                        F1_1(i) = ( U1_1(0,i) - U1_1(1,i) ) / D2
                        F1_2(i) = ( U1_2(0,i) - U1_2(1,i) ) / D2
                        F2_1(i) = ( U2_1(0,i) - U2_1(1,i) ) / D2
                        F2_2(i) = ( U2_2(0,i) - U2_2(1,i) ) / D2
                end do
        end function TPMInteractionFNum !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!---------------------------------------------------------------------------------------------------
! Initialization
!---------------------------------------------------------------------------------------------------

        subroutine TPMInit ( ChiIndM, ChiIndN ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer*4, intent(in)   :: ChiIndM, ChiIndN
        real*8                  :: RT, DX
        !-------------------------------------------------------------------------------------------
                TPPotType = TP_POT_MONO_R
        
                ! Here we calculate the radius of nanotubes
                RT = TPBAcc * sqrt ( 3.0d+00 * ( ChiIndM * ChiIndM + ChiIndN * ChiIndN + ChiIndM * ChiIndN ) ) / M_2PI;
                print *, '(a,i3,a,i3,a,e18.10,a)', 'TPM is iniatized for (', ChiIndM, ',', ChiIndN, ') CNTs, RT = ', RT, ' A'

                TPMChiIndM = ChiIndM
                TPMChiIndN = ChiIndN
                TPMR1 = RT
                TPMR2 = RT

                TPMCaA = 0.275d+00 * ( 1.0d+00 - 1.0d+00 / ( 1.0d+00 + 0.59d+00 * RT ) )
                TPMCeA = 0.35d+00 + 0.0226d+00 * ( RT - 6.785d+00 )
                TPMAHmin0 = 10.0d+00 * ( RT / 6.785d+00 )**1.5

                TPMHmax = TPMR1 + TPMR2 + TPBRcutoff
                TPMDH   = TPMHmax / TPMNH1
                
                ! Parameters of the angle switch
                TPMASMin = sqr ( cos ( rad ( TPMAS ) ) )
                TPMASMax = 1.0d+00 - TPGeomPrec
                TPMASDelta = TPMASMax - TPMASMin
                
                call TPMSSTPInit ()
                
                call TPMSTPInit ()
                
                DX = TPMR1 + TPMR2 + TPBRcutoff
                call TPMAInit ( - DX, DX, - DX, DX )                
                
                call TPMSTInit ()
                
        end subroutine TPMInit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module TubePotMono !****************************************************************************

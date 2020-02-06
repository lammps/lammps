program TMDPotGen !*********************************************************************************
!
! Stand-alone generator of tabulate files with tubular potential data for single-walled CNTs.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
!
!***************************************************************************************************

use TubePotMono

implicit none

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------

        integer*4       :: ChiIndM = 10         ! Chirality index m of nanotubes
        integer*4       :: ChiIndN = 10         ! Chirality index n of nanotubes

!        real*8          :: RT                   ! CNT radius (Angstrom)
        
!---------------------------------------------------------------------------------------------------
! Body
!---------------------------------------------------------------------------------------------------
        
!                BC_X = 0
!                BC_Y = 0
!                BC_Z = 0
                
                TPMStartMode = 0
!                TPMNN = 100
!                TPMHSwitch   = 0       
!                TPMHS        = 3.0d+00
!                TPMASwitch   = 0       
!                TPMAS        = 3.0d+00
                

                ! Reading and printing of governing parameters
                call LoadGoverningParameters () 
                call PrintGoverningParameters () 

                ! Here we calculate the radius of nanotubes
!                RT = TPBAcc * sqrt ( 3.0d+00 * ( ChiIndM * ChiIndM + ChiIndN * ChiIndN + ChiIndM * ChiIndN ) ) / M_2PI;

                call TPBInit ()
                call TPMInit ( ChiIndM, ChiIndN )

contains !------------------------------------------------------------------------------------------

        subroutine LoadGoverningParameters () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function reads governing parameters from xdt file
        !-------------------------------------------------------------------------------------------
        integer*4       :: Fuid, i
        character*512   :: Msg
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDPotGen.xdt', 'rt', '' )
                read ( unit = Fuid, fmt = '(i22)' ) ChiIndM
                read ( unit = Fuid, fmt = '(i22)' ) ChiIndN
                call CloseFile ( Fuid )
        end subroutine LoadGoverningParameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine PrintGoverningParameters () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function prints governing parameters to xlg file
        !-------------------------------------------------------------------------------------------
        integer*4       :: Fuid, i
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDPotGen.xlg', 'wt', '' )
                write ( unit = Fuid, fmt = '(i22,a)' ) ChiIndM, ' : ChiIndM'
                write ( unit = Fuid, fmt = '(i22,a)' ) ChiIndN, ' : ChiIndN'
                call CloseFile ( Fuid )
        end subroutine PrintGoverningParameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program TMDPotGen !*****************************************************************************

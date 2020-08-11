module TMDGenData !*********************************************************************************
!
! Common data for TMDGen
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, Version 13.00, 2020
!
!---------------------------------------------------------------------------------------------------

use TPMGeom

implicit none

!---------------------------------------------------------------------------------------------------
! Constants
!---------------------------------------------------------------------------------------------------

        integer*4, parameter                    :: MAX_TUBE = 1000000 ! Maximum number of tubes in 3D

        real*8, parameter                       :: K_MDDU = K_MDMU / K_MDLU / K_MDLU / K_MDLU ! MD density unit (kg/m**3)

        !
        ! These parameters are specific for carbon nanotubes and taken from module TubePotBase
        !

        real*8, parameter                       :: TPbConstD    = 5.196152422706632d+00 ! = 3.0**1.5
        ! Mass of C atom
        real*8, parameter                       :: TPBM         = 12.0107d+00 ! (a.m.u.)
        ! Lattice parameter and numerical density of atoms for a graphene sheet, see Dresselhaus et al, Carbon 33(7), 1995
        real*8, parameter                       :: TPBA         = 1.421d+00 ! (Angstrom)
        real*8, parameter                       :: TPBD         = 4.0d+00 / ( TPBConstD * TPBA * TPBA ) ! (1/Angstrom^2)
        ! Specific heat of carbon nanotubes
        real*8, parameter                       :: TPBSH        = 600.0d+00 / K_MDCU ! (eV/(Da*K))

!---------------------------------------------------------------------------------------------------
! Governing parameters
!---------------------------------------------------------------------------------------------------
        
        ! Parameters of the sample
        
        real*8                                  :: LS0 = 4000.0         ! Sample size in x, y-directions (Angstrom)
        real*8                                  :: HS0 = 4000.0         ! Sample size in z-direction (Angstrom)
        real*8                                  :: DS0 = 0.01           ! Density (g/cm**3)
        integer*4                               :: BC_X0 = 1            ! Boundary conditions in x-direction: 0, free; 1, periodic
        integer*4                               :: BC_Y0 = 1            ! Boundary conditions in y-direction: 0, free; 1, periodic
        integer*4                               :: BC_Z0 = 1            ! Boundary conditions in z-direction: 0, free; 1, periodic

        ! Parameters of tubes
        integer*4                               :: ChiIndM = 10         ! Chirality index m of nanotubes
        integer*4                               :: ChiIndN = 10         ! Chirality index n of nanotubes
        real*8                                  :: LT0  = 2000.0        ! Characterstic length of tubes (Angstrom)
        integer*4                               :: SegType = 0          ! 0, number of segments per tube is fixed
                                                                        ! 1, rounded length of segments is fixed 
        integer*4                               :: NSeg0 = 100          ! Number of segments per tube
        real*8                                  :: LSeg0 = 20.0d+00     ! Length of the segment (Angstrom)


        ! Parameters controlling the sample structure
        
        real*8                                  :: DeltaT = 3.0         ! Minimal distance between tubes (Angstrom)
        integer*4                               :: NAmax = 50000        ! Maximal number of attempts (for SampleType = 4 it is used as an input paramtere for number of tubes)
        real*8                                  :: GeomPrec = 1.0d-06   ! Geometrical precision

!---------------------------------------------------------------------------------------------------
! Computed data
!---------------------------------------------------------------------------------------------------

        real*8                                  :: RT0  = 6.785         ! Radius of tubes (Angstrom)
        
        real*8                                  :: VS0                  ! Desired volume of the sample, Angstrom**3
        real*8                                  :: MS0                  ! Desired mass of the sample, Da (For SampleType = 4 it is the defined fixed mass- definition is given in TMDGen7T)

        real*8					:: CTCD			! Center to center distance between any surrounding tube and center tube (used for SampleType == 4 only) 
        
        integer*4                               :: NT                   ! Real number of tubes
        real*8, dimension(0:MAX_TUBE-1)         :: RT                   ! Radii of tubes, Angstrom
        real*8, dimension(0:MAX_TUBE-1)         :: LT                   ! Lengths of tubes, Angstrom
        real*8, dimension(0:MAX_TUBE-1,0:2)     :: CT                   ! Coordinates of tubes' centers, Angstrom
        real*8, dimension(0:MAX_TUBE-1,0:2)     :: DT                   ! Directions of tubes
        integer*4, dimension(0:MAX_TUBE-1)      :: AT                   ! Parent axes of tubes. It is used only in GeneratorBundle ()
        
contains !******************************************************************************************

!---------------------------------------------------------------------------------------------------
! Pseudo-random number generator
!---------------------------------------------------------------------------------------------------

        real*8 function randnumber () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns a pseudo-random number with uniform distribution in [0,1]
        !-------------------------------------------------------------------------------------------
                call random_number ( randnumber )
        end function randnumber !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine SetRandomSeed () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This subroutine sets random seed for the pseudo-random number generator
        !-------------------------------------------------------------------------------------------
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
        !-------------------------------------------------------------------------------------------          
            call RANDOM_SEED ( size = n )
            allocate ( seed(n) )
            call SYSTEM_CLOCK ( COUNT = clock )
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call RANDOM_SEED ( PUT = seed )
            deallocate ( seed )
        end subroutine SetRandomSeed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Generators for (random) properties of nanotubes
!---------------------------------------------------------------------------------------------------

        real*8 function TubeMass ( i ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the mass of the tube in Da       
        !-------------------------------------------------------------------------------------------
        integer*4, intent(in)   :: i
        !-------------------------------------------------------------------------------------------
                TubeMass = M_2PI * RT(i) * LT(i) * TPBM * TPBD
        end function TubeMass !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real*8 function TubeSpecificHeat ( i ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns the specific heat of the tube
        !-------------------------------------------------------------------------------------------
        integer*4, intent(in)   :: i
        !-------------------------------------------------------------------------------------------
                TubeSpecificHeat = TPBSH
        end function TubeSpecificHeat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Reading and printing of input parameters
!---------------------------------------------------------------------------------------------------
        
        subroutine LoadGoverningParameters () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function reads governing parameters from xdt file
        !-------------------------------------------------------------------------------------------
        integer*4       :: Fuid, i
        character*512   :: Msg
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDGen.xdt', 'rt', '' )
                read ( unit = Fuid, fmt = '(e22.12)' ) LS0
                read ( unit = Fuid, fmt = '(e22.12)' ) HS0
                read ( unit = Fuid, fmt = '(e22.12)' ) DS0
                read ( unit = Fuid, fmt = '(i22)' ) BC_X0
                read ( unit = Fuid, fmt = '(i22)' ) BC_Y0
                read ( unit = Fuid, fmt = '(i22)' ) BC_Z0
                read ( unit = Fuid, fmt = '(i22)' ) ChiIndM
                read ( unit = Fuid, fmt = '(i22)' ) ChiIndN
                read ( unit = Fuid, fmt = '(e22.12)' ) LT0
                read ( unit = Fuid, fmt = '(i22)' ) SegType
                read ( unit = Fuid, fmt = '(i22)' ) NSeg0
                read ( unit = Fuid, fmt = '(e22.12)' ) LSeg0
                read ( unit = Fuid, fmt = '(e22.12)' ) DeltaT
                read ( unit = Fuid, fmt = '(i22)' ) NAmax
                read ( unit = Fuid, fmt = '(e22.12)' ) GeomPrec
                call CloseFile ( Fuid )
        end subroutine LoadGoverningParameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine PrintGoverningParameters () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function prints governing parameters to xlg file
        !-------------------------------------------------------------------------------------------
        integer*4       :: Fuid, i
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDGen.xlg', 'wt', '' )
                write ( unit = Fuid, fmt = '(e22.12,a)' ) LS0, ' : LS0, Angstrom'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) HS0, ' : HS0, Angstrom'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) DS0, ' : DS0, g/cm**3'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) DS0, ' : SC0, 1/A**2'
                write ( unit = Fuid, fmt = '(i22,a)' ) BC_X0, ' : BC_X0'
                write ( unit = Fuid, fmt = '(i22,a)' ) BC_Y0, ' : BC_Y0'
                write ( unit = Fuid, fmt = '(i22,a)' ) BC_Z0, ' : BC_Z0'
                write ( unit = Fuid, fmt = '(i22,a)' ) ChiIndM, ' : ChiIndM'
                write ( unit = Fuid, fmt = '(i22,a)' ) ChiIndN, ' : ChiIndN'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) LT0, ' : LT0, Angstrom'
                write ( unit = Fuid, fmt = '(i22,a)' ) SegType, ' : SegType'
                write ( unit = Fuid, fmt = '(i22,a)' ) NSeg0, ' : NSeg0'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) LSeg0, ' : LSeg0, Angstrom'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) DeltaT, ' : DeltaT'
                write ( unit = Fuid, fmt = '(i22,a)' ) NAmax, ' : NAmax'
                write ( unit = Fuid, fmt = '(e22.12,a)' ) GeomPrec, ' : GeomPrec'
                call CloseFile ( Fuid )
        end subroutine PrintGoverningParameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Printing of sample parameters
!---------------------------------------------------------------------------------------------------
        
        subroutine PrintSampleParameters ( ParType ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function prints the most imprtant parameters of the sample.
        ! In the code, it used twice to print parameters of the desired and really generated samples.
        !-------------------------------------------------------------------------------------------
        character*(*), intent(in)       :: ParType
        real*8                          :: MP, M, V
        !-------------------------------------------------------------------------------------------
                print '(a,a,a)', '*** ', trim(ParType), ' properties of the sample'
                print '(a34,a,f15.4,a)', 'L', ' : ', LS0, ' A'
                print '(a34,a,f15.4,a)', 'H', ' : ', HS0, ' A'
                print '(a34,a,f15.4,a)', 'Density', ' : ', DS0, ' g/cm**3'
                print '(a34,a,e15.8,a)', 'Volume', ' : ', VS0, ' A*3'
                print '(a34,a,e15.8,a)', 'Mass', ' : ', MS0, ' Da'
                print '(a34,a,i10)', 'BC_X', ' : ', BC_X0
                print '(a34,a,i10)', 'BC_Y', ' : ', BC_Y0
                print '(a34,a,i10)', 'BC_Z', ' : ', BC_Z0
        end subroutine PrintSampleParameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! Initializing of basic geometrical parameters of the generated sample
!---------------------------------------------------------------------------------------------------        
        
        subroutine InitSample () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function initializes the geometrical parameters of the sample (sizes, etc.) 
        !-------------------------------------------------------------------------------------------

                BC_X = BC_X0
                BC_Y = BC_Y0
                BC_Z = BC_Z0
                DomXmin = - LS0 / 2.0d+00
                DomXmax = LS0 / 2.0d+00
                DomYmin = - LS0 / 2.0d+00
                DomYmax = LS0 / 2.0d+00
                DomZmin = - HS0 / 2.0d+00
                DomZmax = HS0 / 2.0d+00
                
                if ( BC_X0 == 0 ) then
                        DomXmin = 0.0d+00
                        DomXmax = LS0
                end if
                if ( BC_Y0 == 0 ) then
                        DomYmin = 0.0d+00
                        DomYmax = LS0
                end if
                if ( BC_Z0 == 0 ) then
                        DomZmin = 0.0d+00
                        DomZmax = HS0
                end if
        
                DomLX      = DomXmax - DomXmin
                DomLY      = DomYmax - DomYmin
                DomLZ      = DomZmax - DomZmin
                DomLXHalf  = 0.5d+00 * DomLX 
                DomLYHalf  = 0.5d+00 * DomLY
                DomLZHalf  = 0.5d+00 * DomLZ 
                
                DS0 = DS0 / ( K_MDDU / 1.0d+03 )
                VS0 = LS0 * LS0 * HS0
                MS0 = DS0 * VS0
        end subroutine InitSample !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------------
! A few auxiliary functions
!---------------------------------------------------------------------------------------------------

        subroutine GetTubeEnds ( X0, X1, i ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function calculates coordinates of two ends of nanotube i
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2), intent(out)     :: X0, X1
        integer*4, intent(in)                   :: i
        !-------------------------------------------------------------------------------------------
        real*8                                  :: LT2
        !-------------------------------------------------------------------------------------------
                LT2 = 0.5d+00 * LT(i)
                X0 = CT(i,0:2) - LT2 * DT(i,0:2)
                X1 = CT(i,0:2) + LT2 * DT(i,0:2)
        end subroutine GetTubeEnds !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        logical function IsTubeInside ( i ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function returns true if nanotube i lies inside the sample. Otherwise it returns false.
        !-------------------------------------------------------------------------------------------
        integer*4, intent(in)   :: i
        !-------------------------------------------------------------------------------------------
        integer*4               :: n
        real*8, dimension(0:2)  :: X0, X1, Xmin, Xmax
        !-------------------------------------------------------------------------------------------
                IsTubeInside = .true.
                if ( BC_X == 1 .and. BC_Y == 1 .and. BC_Z == 1 ) return
                call GetTubeEnds ( X0, X1, i )
                do n = 0, 2
                        Xmin(n) = min ( X0(n), X1(n) )
                        Xmax(n) = max ( X0(n), X1(n) )
                end do
                IsTubeInside = .false.
                if ( BC_X == 0 .and. ( Xmin(0) < DomXmin .or. Xmax(0) > DomXmax ) ) return
                if ( BC_Y == 0 .and. ( Xmin(1) < DomYmin .or. Xmax(1) > DomYmax ) ) return
                if ( BC_Z == 0 .and. ( Xmin(2) < DomZmin .or. Xmax(2) > DomZmax ) ) return
                IsTubeInside = .true.
        end function IsTubeInside !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TMDGenData !*****************************************************************************

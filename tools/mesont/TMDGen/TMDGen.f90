program TMDGen !************************************************************************************
!
! Stand-alone generator of 3D CNT samples.
!
!---------------------------------------------------------------------------------------------------
!
! Intel Fortran
!
! Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
!
!***************************************************************************************************

use TMDGen3D

implicit none

!---------------------------------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------------------------------

        integer*4       :: Nseg, Nnode
        real*8          :: DS00

!---------------------------------------------------------------------------------------------------
! Body
!---------------------------------------------------------------------------------------------------
        
                print *, 'TMD generator of 3D CNT samples, v. 13.00'
                print '(a34,a,i10)', 'Maximum number of nanotubes', ' : ', MAX_TUBE
                
                call SetRandomSeed ()

                ! Reading and printing of governing parameters
                call LoadGoverningParameters () 
                call PrintGoverningParameters () 

                ! Here we calculate the radius of nanotubes
                RT0 = TPBA * sqrt ( 3.0d+00 * ( ChiIndM * ChiIndM + ChiIndN * ChiIndN + ChiIndM * ChiIndN ) ) / M_2PI;

                ! Here we calculate parameters of the desired sample
                call InitSample ()
                DS0  = DS0 * ( K_MDDU / 1.0d+03 )
                call PrintSampleParameters ( 'Desired' )
                DS00 = DS0
                DS0  = DS0 / ( K_MDDU / 1.0d+03 )
                
                call Generator3D ()

                ! Here we write the major output file with the sample
                !call WriteOutputFile_old_format ()
                !call WriteOutputFile ()
        
                ! Here we write an auxiliary Tecplot file to visualize the initial sample
                !PrintTecplotFile ()    
                call WriteLAMMPSFile()

                ! Here we print parameters of the final sample
                call PrintSampleParameters ( 'Final' )
                print '(a34,a,f15.4,a)', 'Nanotube radius ', ' : ', RT0, ' a'
                print '(a34,a,f15.4,a)', 'Nanotube length ', ' : ', LT0, ' a'
                print '(a34,a,f15.4,a)', 'Nanotube mass ', ' : ', M_2PI * RT0 * LT0 * TPBM * TPBD, ' Da'
                if ( SegType == 0 ) then
                        LSeg0 = LT0 / NSeg0
                else
                        NSeg0 = int ( LT0 / LSeg0 ) + 1
                        LSeg0 = LT0 / NSeg0
                end if
                print '(a34,a,f15.4,a)', 'Nanotube segment length ', ' : ', LSeg0, ' a'
                print '(a34,a,f15.4,a)', 'Nanotube segment mass ', ' : ', M_2PI * RT0 * LSeg0 * TPBM * TPBD, ' Da'
                print '(a34,a,f15.4)', 'Desired / Real densities ', ' : ', DS00 / DS0
                print '(a34,a,i10)', 'Real number of tubes', ' : ', NT
                print '(a34,a,i10)', 'Real number of segments', ' : ', Nseg
                print '(a34,a,i10)', 'Real number of nodes', ' : ', Nnode

contains !******************************************************************************************

        subroutine DiscretizeTube ( X0, DL, NS, i ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function calculaats geometrical parameters that are necessary to represent straight 
        ! tube i as a sequence of segments.
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2), intent(out)     :: X0
        real*8, intent(out)                     :: DL
        integer*4, intent(out)                  :: NS
        integer*4, intent(in)                   :: i
        !-------------------------------------------------------------------------------------------
        real*8, dimension(0:2)                  :: X1
        !-------------------------------------------------------------------------------------------
                call GetTubeEnds ( X0, X1, i )
                if ( SegType == 0 ) then
                        NS = NSeg0
                else
                        NS = int ( LT(i) / LSeg0 ) + 1
                end if
                DL = LT(i) / NS
        end subroutine DiscretizeTube !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine WriteOutputFile_old_format () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function writes a dat file (version 2) with the initial nanotube sample.
        ! This file is used by TMD/TMDMPI to start a new simulation.
        !-------------------------------------------------------------------------------------------
        integer*4               :: Fuid, i, j, NTS, Prop
        real*8                  :: DL, L, L00, M00, I00, J00, C00, LL00, MM00, II00, JJ00, CC00
        real*8, dimension(0:2)  :: X, X0
        logical*4               :: PrintNode                
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDGen_old.dat', "wt", "" )
                write ( unit = Fuid, fmt = '(i12)' )  3
                write ( unit = Fuid, fmt = '(2i4,4e20.12)' )  ChiIndM, ChiIndN, RT0, TPBA, TPBD, TPBM
                write ( unit = Fuid, fmt = '(3e20.12)' )  DomXmin, DomYmin, DomZmin
                write ( unit = Fuid, fmt = '(3e20.12)' )  DomXmax, DomYmax, DomZmax
                write ( unit = Fuid, fmt = '(3i12)' ) BC_X, BC_Y, BC_Z 
                write ( unit = Fuid, fmt = '(i12)' ) NT
                Nseg = 0
                Nnode = 0
                do i = 0, NT - 1
                        call DiscretizeTube ( X0, DL, NTS, i )
                        L00 = LT(i) / NTS
                        M00 = TubeMass ( i ) / NTS
                        I00 = 0.0d+00
                        J00 = M00 * sqr ( RT(i) )
                        C00 = M00 * TubeSpecificHeat ( i )
                        Nseg = Nseg + NTS
                        write ( unit = Fuid, fmt = '(i12)' ) NTS + 1 
                        Nnode = Nnode + NTS + 1
                        L = 0.0d+00
                        do j = 0, NTS 
                                X = X0 + L * DT(i,0:2)
                                MM00 = M00
                                II00 = I00
                                JJ00 = J00
                                CC00 = C00
                                LL00 = L00
                                if ( j == 0 .or. j == NTS ) then
                                        MM00 = 0.5d+00 * M00
                                        II00 = 0.5d+00 * I00
                                        JJ00 = 0.5d+00 * J00
                                        CC00 = 0.5d+00 * C00
                                end if
                                if ( j == NTS ) LL00 = 0.0d+00 
                                Prop = 0
                                write ( unit = Fuid, fmt = '(i2,6e20.12)' ) Prop, RT(0), LL00, MM00, II00, JJ00, CC00
                                write ( unit = Fuid, fmt = '(6e20.12)' ) X, RT(i), 0.0d+00, 300.0d+00
                                L = L + DL
                        end do
                end do
                write ( unit = Fuid, fmt = '(i12)' ) 0
                write ( unit = Fuid, fmt = '(i12)' ) 0
                call CloseFile ( Fuid )
        end subroutine WriteOutputFile_old_format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine WriteOutputFile () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function writes a dat file (version 2) with the initial nanotube sample.
        ! This file is used by TMD/TMDMPI to start a new simulation.
        !-------------------------------------------------------------------------------------------
        integer*4               :: Fuid, i, j, NTS
        real*8                  :: DL, L, L00, M00, LL00, MM00
        real*8, dimension(0:2)  :: X, X0
        logical*4               :: PrintNode                
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDGen.dat', "wt", "" )
                write ( unit = Fuid, fmt = '(2i4,4e20.12)' )  ChiIndM, ChiIndN, RT0, TPBA, TPBD, TPBM
                write ( unit = Fuid, fmt = '(3e20.12)' )  DomXmin, DomYmin, DomZmin
                write ( unit = Fuid, fmt = '(3e20.12)' )  DomXmax, DomYmax, DomZmax
                write ( unit = Fuid, fmt = '(3i12)' ) BC_X, BC_Y, BC_Z 
                write ( unit = Fuid, fmt = '(i12)' ) NT
                Nseg = 0
                Nnode = 0
                do i = 0, NT - 1
                        call DiscretizeTube ( X0, DL, NTS, i )
                        L00 = LT(i) / NTS
                        M00 = TubeMass ( i ) / NTS
                        Nseg = Nseg + NTS
                        write ( unit = Fuid, fmt = '(i12)' ) NTS + 1 
                        Nnode = Nnode + NTS + 1
                        L = 0.0d+00
                        do j = 0, NTS 
                                X = X0 + L * DT(i,0:2)
                                MM00 = M00
                                LL00 = L00
                                if ( j == 0 .or. j == NTS ) MM00 = 0.5d+00 * M00
                                if ( j == NTS ) LL00 = 0.0d+00 
                                write ( unit = Fuid, fmt = '(5e20.12)' ) X, LL00, MM00
                                L = L + DL
                        end do
                end do
                call CloseFile ( Fuid )
        end subroutine WriteOutputFile !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine PrintTecplotFile () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function prints Tecplot file to visualize the generated sample
        !-------------------------------------------------------------------------------------------
        integer*4                       :: Fuid, i
        real*8                          :: LT2
        !-------------------------------------------------------------------------------------------
                Fuid = OpenFile ( 'TMDGen.plt', "wt", "" )
                write ( unit = Fuid, fmt = '(a)' ) 'VARIABLES="X" "Y" "Z"'
                do i = 0, NT - 1
                        write ( unit = Fuid, fmt = '(a,i,a)' ) 'ZONE T="T', i, '"'
                        LT2 = 0.5d+00 * LT(i)
                        write ( unit = Fuid, fmt = '(3e20.12)' ) CT(i,0:2) - LT2 * DT(i,0:2)
                        write ( unit = Fuid, fmt = '(3e20.12)' ) CT(i,0:2) + LT2 * DT(i,0:2)
                end do
                call CloseFile ( Fuid )
        end subroutine PrintTecplotFile !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        subroutine WriteLAMMPSFile () !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This function writes a dat file (version 2) with the initial nanotube sample.
        ! This file is used by TMD/TMDMPI to start a new simulation.
        !-------------------------------------------------------------------------------------------
            integer*4               :: file_id, i, j, NTS, node_id, b1, b2
            real*8                  :: DL, L, L00, M00, LL00, MM00
            real*8, dimension(0:2)  :: X, X0
            logical*4               :: PrintNode                
            !-------------------------------------------------------------------------------------------
            open(newunit = file_id, file = 'TMDSample.init')
            write(file_id,*)
            write(file_id,*)
            !count the number of nodes and segments
            Nseg = 0
            Nnode = 0
            do i = 0, NT - 1
                call DiscretizeTube (X0, DL, NTS, i)
                Nseg = Nseg + NTS
                Nnode = Nnode + NTS + 1
            enddo
            write(file_id,'(i9,a)') Nnode, " atoms"
            write(file_id,*) 
            write(file_id,*) "1 atom types"
            write(file_id,*) 
            write(file_id,'(2e20.12,2a)') DomXmin, DomXmax, " xlo xhi"
            write(file_id,'(2e20.12,2a)') DomYmin, DomYmax, " ylo yhi"
            write(file_id,'(2e20.12,2a)') DomZmin, DomZmax, " zlo zhi"
            write(file_id,*) 
            write(file_id,*) "Masses"
            write(file_id,*) 
            write(file_id,*) "1 1.0"
            write(file_id,*) 
            write(file_id,*) "Atoms"
            write(file_id,*) 
            
            node_id = 1
            do i = 0, NT - 1
                call DiscretizeTube(X0, DL, NTS, i)
                L00 = LT(i) / NTS
                M00 = TubeMass (i) / NTS
                
                b1 = -1
                L = 0.0d+00
                do j = 0, NTS 
                    b2 = node_id + 1
                    if (j == NTS) b2 = -1
                    MM00 = M00
                    LL00 = L00
                    if (j == 0 .or. j == NTS) MM00 = 0.5d+00 * M00
                    if (j == NTS) LL00 = 0.0d+00 
                    X = X0 + L * DT(i,0:2)
                    write(file_id,'(2i9,a,2i9,3e14.7,a,3e20.12,a)') node_id, i, " 1 ", b1, b2, MM00, RT(i), LL00, " 0 ", X, " 0 0 0"
                    b1 = node_id
                    node_id = node_id + 1
                    L = L + DL
                enddo
            enddo
            close(file_id)
        end subroutine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end program TMDGen !********************************************************************************

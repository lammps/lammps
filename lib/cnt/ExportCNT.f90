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
!   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
!------------------------------------------------------------------------- 

module ExportCNT !*******************************************************************************
    use iso_c_binding
    use CNTPot
    use TPMLib
    use TubePotMono
    use TPMForceField
    implicit none

contains
    subroutine InitCNTPotModule_(STRModel, STRParams, YMType, BNDModel, Rref) &
        bind(c, name = "InitCNTPotModule")
        integer*4, intent(in)   :: STRModel, STRParams, YMType, BNDModel
        real*8, intent(in)      :: Rref
        
        call InitCNTPotModule(STRModel, STRParams, YMType, BNDModel, Rref)
    endsubroutine
    
    subroutine TPBInit_() &
        bind(c, name = "TPBInit")
        
        call TPBInit()
    endsubroutine

    subroutine TPMInit_(M, N) &
        bind(c, name = "TPMInit")
        integer*4, intent(in)      :: M, N
        
        call TPMInit(M, N)
    endsubroutine
    
    subroutine SetTablePath_(TPMSSTPFile_, N1, TPMAFile_, N2) &
        bind(c, name = "SetTablePath")
        integer*4, intent(in)                           :: N1, N2
        character, intent(in), dimension(N1)            :: TPMSSTPFile_
        character, intent(in), dimension(N2)            :: TPMAFile_
        integer                                         :: i
        
        do i = 1, len(TPMSSTPFile)
            if (i <= N1) then
                TPMSSTPFile(i:i) = TPMSSTPFile_(i)
            else 
                TPMSSTPFile(i:i) = ' '
            endif
        enddo
        do i = 1, len(TPMAFile)
            if (i <= N2) then
                TPMAFile(i:i) = TPMAFile_(i)
            else 
                TPMAFile(i:i) = ' '
            endif
        enddo
    endsubroutine
    
    function get_R_ () &
        bind(c, name = "get_R")
        real*8 :: get_R_
        get_R_ = TPMR1
        return
    endfunction
        

    subroutine TubeStretchingForceField_(U1, U2, F1, F2, S1, S2, X1, X2, R12, L12) &
        bind(c, name = "TubeStretchingForceField")
        real*8, intent(inout)                           :: U1, U2       ! Interaction energies associated with nodes X1 and X2
        real*8, intent(inout), dimension(0:2)           :: F1, F2       ! Forces exerted on nodes X1 and X2
        real*8, intent(inout), dimension(0:2,0:2)       :: S1, S2       ! Contributions of nodes X1 and X2 to the virial stress tensor
        real*8, intent(in), dimension(0:2)              :: X1, X2       ! Coordinates of the segmnet nodes 
        real*8, intent(in)                              :: R12          ! Radius of nanotube the segment (X1,X2) belongs to
        real*8, intent(in)                              :: L12          ! Equilubrium length of segment (X1,X2)
        
        call TubeStretchingForceField(U1, U2, F1, F2, S1, S2, X1, X2, R12, L12)
    endsubroutine

    subroutine TubeBendingForceField_(U1, U2, U3, F1, F2, F3, S1, S2, S3, X1, X2, X3, R123, L123, BBF2) &
        bind(c, name = "TubeBendingForceField")
        real*8, intent(inout)                           :: U1, U2, U3   ! Interaction energies associated with nodes X1, X2, and X3
        real*8, intent(inout), dimension(0:2)           :: F1, F2, F3   ! Forces exerted on nodes X1, X2, and X3
        real*8, intent(inout), dimension(0:2,0:2)       :: S1, S2, S3   ! Contributions of nodes X1, X2, and X3 to the virial stress tensor
        real*8, intent(in), dimension(0:2)              :: X1, X2, X3   ! Coordinates of nodes 
        real*8, intent(in)                              :: R123         ! Radius of nanotube the segment (X1,X2) belongs to
        real*8, intent(in)                              :: L123         ! Equilubrium length of segment (X1,X2) and (X2,X3) (It is assumed to be the same for both segments)
        integer*4, intent(inout)                        :: BBF2
        
        call TubeBendingForceField(U1, U2, U3, F1, F2, F3, S1, S2, S3, X1, X2, X3, R123, L123, BBF2 )
    endsubroutine

    subroutine SegmentTubeForceField_(U1, U2, U, F1, F2, F, Fe, S1, S2, S, Se, X1, X2, R12, N, X, Xe, BBF, R, E1, E2, Ee, TPMType) &
        bind(c, name = "SegmentTubeForceField")
        integer*4, intent(in)                           :: N            ! Number of nodes in array X
        real*8, intent(inout)                           :: U1, U2       ! Interaction energies associated with nodes X1 and X2
        real*8, intent(inout), dimension(0:N-1)         :: U            ! Interaction energies associated with nodes X
        real*8, intent(inout), dimension(0:2)           :: F1, F2       ! Forces exerted on nodes X1 and X2
        real*8, intent(inout), dimension(0:2,0:N-1)     :: F            ! Forces exerted on nodes X
        real*8, intent(inout), dimension(0:2)           :: Fe           ! Force exerted on node Xe (can be updated only if Ee > 0)
        real*8, intent(inout), dimension(0:2,0:2)       :: S1, S2       ! Contributions of nodes X1 and X2 to the virial stress tensor
        real*8, intent(inout), dimension(0:2,0:2,0:N-1) :: S            ! Contributions of nodes X to the virial stress tensor
        real*8, intent(inout), dimension(0:2,0:2)       :: Se           ! Contributions of node Xe to the virial stress tensor (can be updated only if Ee > 0)
        real*8, intent(in), dimension(0:2)              :: X1, X2       ! Coordinates of the segmnet nodes 
        real*8, intent(in)                              :: R12          ! Radius of nanotube the segment (X1,X2) belongs to
        real*8, intent(in), dimension(0:2,0:N-1)        :: X            ! Coordinates of the nanotube nodes
        real*8, intent(in), dimension(0:2)              :: Xe           ! Additional node of the extended chain if Ee > 0
        integer*4, intent(in), dimension(0:N-1)         :: BBF          ! Bending buckling flags (BBF(i) = 1 in a case of buckling in node i)
        real*8, intent(in)                              :: R            ! Radius of nanotube X
        integer*4, intent(in)                           :: E1, E2       ! E1 = 1 if the chnane node 0 is a CNT end; E2 = 1 if the chnane node N-1 is a CNT end;
        integer*4, intent(in)                           :: Ee           ! Parameter defining the type of the extended chain (0,1,2)
        integer*4, intent(in)                           :: TPMType      ! Type of the tubular potential (0 or 1)
        
        call SegmentTubeForceField(U1, U2, U, F1, F2, F, Fe, S1, S2, S, Se, X1, X2, R12, N, X, Xe, BBF, R, E1, E2, Ee, TPMType)
    endsubroutine
endmodule ExportCNT !**************************************************************************

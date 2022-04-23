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
!   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
!-------------------------------------------------------------------------

module ExportCNT !**********************************************************************************

    use iso_c_binding
    use CNTPot
    use TPMLib
    use TubePotMono
    use TPMForceField
    implicit none

contains

    subroutine InitCNTPotModule_(STRModel, STRParams, YMType, BNDModel, Rref) &
        bind(c, name = "mesont_lib_InitCNTPotModule")
        integer(c_int), intent(in)      :: STRModel, STRParams, YMType, BNDModel
        real(c_double), intent(in)      :: Rref

        call InitCNTPotModule(STRModel, STRParams, YMType, BNDModel, Rref)
    endsubroutine

    subroutine TPBInit_() &
        bind(c, name = "mesont_lib_TPBInit")

        call TPBInit()
    endsubroutine

    subroutine TPMInit_(M, N) &
        bind(c, name = "mesont_lib_TPMInit")
        integer(c_int), intent(in) :: M, N

        call TPMInit(M, N)
    endsubroutine

    subroutine SetTablePath_(TPMFile_, N) &
        bind(c, name = "mesont_lib_SetTablePath")
        integer(c_int), intent(in)                      :: N
        character(c_char), intent(in), dimension(N)     :: TPMFile_
        integer                                         :: i

        do i = 1, len(TPMFile)
            if (i <= N) then
                TPMFile(i:i) = TPMFile_(i)
            else
                TPMFile(i:i) = ' '
            endif
        enddo
    endsubroutine

    function get_R_ () &
        bind(c, name = "mesont_lib_get_R")
        real(c_double) :: get_R_
        get_R_ = TPMR1
        return
    endfunction


    subroutine TubeStretchingForceField_(U1, U2, F1, F2, S1, S2, X1, X2, R12, L12) &
        bind(c, name = "mesont_lib_TubeStretchingForceField")
        ! Interaction energies associated with nodes X1 and X2
        real(c_double), intent(inout)                           :: U1, U2
        ! Forces exerted on nodes X1 and X2
        real(c_double), intent(inout), dimension(0:2)           :: F1, F2
        ! Contributions of nodes X1 and X2 to the virial stress tensor
        real(c_double), intent(inout), dimension(0:2,0:2)       :: S1, S2
        ! Coordinates of the segment nodes
        real(c_double), intent(in), dimension(0:2)              :: X1, X2
        ! Radius of a nanotube the segment (X1,X2) belongs to
        real(c_double), intent(in)                              :: R12
        ! Equilibrium length of segment (X1,X2)
        real(c_double), intent(in)                              :: L12

        call TubeStretchingForceField(U1, U2, F1, F2, S1, S2, X1, X2, R12, L12)
    endsubroutine

    subroutine TubeBendingForceField_(U1, U2, U3, F1, F2, F3, S1, S2, S3, X1, X2, X3, R123, L123, BBF2) &
        bind(c, name = "mesont_lib_TubeBendingForceField")
        ! Interaction energies associated with nodes X1, X2, and X3
        real(c_double), intent(inout)                           :: U1, U2, U3
        ! Forces exerted on nodes X1, X2, and X3
        real(c_double), intent(inout), dimension(0:2)           :: F1, F2, F3
        ! Contributions of nodes X1, X2, and X3 to the virial stress tensor
        real(c_double), intent(inout), dimension(0:2,0:2)       :: S1, S2, S3
        ! Coordinates of nodes
        real(c_double), intent(in), dimension(0:2)              :: X1, X2, X3
        ! Radius of nanotube the segment (X1,X2) belongs to
        real(c_double), intent(in)                              :: R123
        ! Equilibrium length of segment (X1,X2) and (X2,X3) (It is assumed to be the same for both segments)
        real(c_double), intent(in)                              :: L123
        integer(c_int), intent(inout)                           :: BBF2

        call TubeBendingForceField(U1, U2, U3, F1, F2, F3, S1, S2, S3, X1, X2, X3, R123, L123, BBF2 )
    endsubroutine

    subroutine SegmentTubeForceField_(U1,U2,U,F1,F2,F,Fe,S1,S2,S,Se,X1,X2,R12,N,X,Xe,BBF,R,E1,E2,Ee,TPMType)&
        bind(c, name = "mesont_lib_SegmentTubeForceField")
        ! Number of nodes in array X
        integer(c_int), intent(in)                              :: N
        ! Interaction energies associated with nodes X1 and X2
        real(c_double), intent(inout)                           :: U1, U2
        ! Interaction energies associated with nodes X
        real(c_double), intent(inout), dimension(0:N-1)         :: U
        ! Forces exerted on nodes X1 and X2
        real(c_double), intent(inout), dimension(0:2)           :: F1, F2
        ! Forces exerted on nodes X
        real(c_double), intent(inout), dimension(0:2,0:N-1)     :: F
        ! Force exerted on node Xe (can be updated only if Ee > 0)
        real(c_double), intent(inout), dimension(0:2)           :: Fe
        ! Contributions of nodes X1 and X2 to the virial stress tensor
        real(c_double), intent(inout), dimension(0:2,0:2)       :: S1, S2
        ! Contributions of nodes X to the virial stress tensor
        real(c_double), intent(inout), dimension(0:2,0:2,0:N-1) :: S
        ! Contributions of node Xe to the virial stress tensor (can be updated only if Ee > 0)
        real(c_double), intent(inout), dimension(0:2,0:2)       :: Se
        ! Coordinates of the segment nodes
        real(c_double), intent(in), dimension(0:2)              :: X1, X2
        ! Radius of nanotube the segment (X1,X2) belongs to
        real(c_double), intent(in)                              :: R12
        ! Coordinates of the nanotube nodes
        real(c_double), intent(in), dimension(0:2,0:N-1)        :: X
        ! Additional node of the extended chain if Ee > 0
        real(c_double), intent(in), dimension(0:2)              :: Xe
        ! Bending buckling flags (BBF(i) = 1 in a case of buckling in node i)
        integer(c_int), intent(in), dimension(0:N-1)            :: BBF
        ! Radius of nanotube X
        real(c_double), intent(in)                              :: R
        ! E1 = 1 if the chain node 0 is a CNT end; E2 = 1 if the chain node N-1 is a CNT end;
        integer(c_int), intent(in)                           :: E1, E2
        ! Parameter defining the type of the extended chain (0,1,2)
        integer(c_int), intent(in)                           :: Ee
        ! Type of the tubular potential (0 or 1)
        integer(c_int), intent(in)                           :: TPMType

        call SegmentTubeForceField(U1, U2, U, F1, F2, F, Fe, S1, S2, S, Se, X1, X2, R12, N, X, Xe, BBF, R, E1, E2, Ee, TPMType)
    endsubroutine

endmodule ExportCNT !*******************************************************************************

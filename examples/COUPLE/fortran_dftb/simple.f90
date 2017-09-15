   module callback
     implicit none
     contains
       subroutine fortran_callback(lmp, timestep, nlocal, ids, c_pos, c_fext) & 
     & bind(C, name='f_callback')
       use, intrinsic :: ISO_C_binding
       use LAMMPS
       implicit none
       type (C_ptr), value :: lmp
       integer(C_int64_t), intent(in), value :: timestep
       integer(C_int), intent(in), value :: nlocal
       real (C_double), dimension(:,:), pointer :: x
       type(c_ptr) :: c_pos, c_fext, c_ids
       double precision, pointer :: fext(:,:), pos(:,:)
       integer, intent(in) :: ids(nlocal)
       real(C_double) :: virial(6)
       real (C_double) :: etot
       real(C_double), pointer :: ts_lmp
       double precision :: stress(3,3), ts_dftb
       integer :: natom , i
       real (C_double), parameter :: econv = 627.4947284155114 ! converts from Ha to
       double precision, parameter :: fconv = 1185.793095983065 ! converts from Ha/bohr to
       double precision, parameter :: autoatm = 2.9037166638E8
       double precision lx, ly, lz
       real (C_double), pointer :: boxxlo, boxxhi
       real (C_double), pointer :: boxylo, boxyhi
       real (C_double), pointer :: boxzlo, boxzhi
       double precision, parameter :: nktv2p = 68568.4149999999935972
       double precision :: volume
       type (C_ptr) :: Cptr
       type (C_ptr), pointer, dimension(:) :: Catom

       call c_f_pointer(c_pos, pos, [3,nlocal])
       call c_f_pointer(c_fext, fext, [3,nlocal])
       call lammps_extract_global(boxxlo, lmp, 'boxxlo')
       call lammps_extract_global(boxxhi, lmp, 'boxxhi')
       call lammps_extract_global(boxylo, lmp, 'boxylo')
       call lammps_extract_global(boxyhi, lmp, 'boxyhi')
       call lammps_extract_global(boxzlo, lmp, 'boxzlo')
       call lammps_extract_global(boxzhi, lmp, 'boxzhi')
       lx = boxxhi - boxxlo
       ly = boxyhi - boxylo
       lz = boxzhi - boxzlo
       volume = lx*ly*lz
       open (unit = 10, status = 'replace', action = 'write', file='lammps.gen')
       write(10,*)nlocal,"S"
       write(10,*) "C"
       do i = 1, nlocal
         write(10,'(2I,3F15.6)')i,1,pos(:,ids(i))
       enddo
       write(10,*)"0.0 0.0 0.0"
       write(10,*)lx,0,0
       write(10,*)0,ly,0
       write(10,*)0,0,lz
       close(10)
       call system("./dftb+ > dftb.out")
       open (unit = 10, status = 'old', file = 'results.out')
       read(10,*)etot
       read(10,*)ts_dftb
       do i = 1, 3
         read(10,*)stress(i,:)
       enddo
       stress (:,:) = stress(:,:)*autoatm
       virial(1) = stress(1,1)/(nktv2p/volume)
       virial(2) = stress(2,2)/(nktv2p/volume)
       virial(3) = stress(3,3)/(nktv2p/volume)
       virial(4) = stress(1,2)/(nktv2p/volume)
       virial(5) = stress(1,3)/(nktv2p/volume)
       virial(6) = stress(2,3)/(nktv2p/volume)
       etot = etot*econv
       call lammps_set_external_vector(lmp,1,ts_dftb*econv)
       do i = 1, nlocal
         read(10,*)fext(:,ids(i))
         fext(:,ids(i)) = fext(:,ids(i))*fconv
       enddo
       close(10)
       call lammps_set_user_energy (lmp, etot)
       call lammps_set_user_virial (lmp, virial)

       end  subroutine
    end module callback


program simple_fortran_callback

   use MPI
   use LAMMPS
   use callback
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int, C_FUNPTR
   implicit none
   type (C_ptr) :: lmp
   integer :: error, narg, me, nprocs

   call MPI_Init (error)
   call MPI_Comm_rank (MPI_COMM_WORLD, me, error)
   call MPI_Comm_size (MPI_COMM_WORLD, nprocs, error)

   call lammps_open_no_mpi ('lmp -log log.simple', lmp)
   call lammps_file (lmp, 'in.simple')
   call lammps_set_callback(lmp)
   call lammps_set_external_vector_length(lmp,2)

   call lammps_command (lmp, 'run 10')
   call lammps_close (lmp)
   call MPI_Finalize (error)


end program simple_fortran_callback



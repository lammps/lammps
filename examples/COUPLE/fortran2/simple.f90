program simple

   use MPI
   use LAMMPS

   ! The following line is unnecessary, as I have included these three entities
   ! with the LAMMPS module, but I leave them in anyway to remind people where
   ! they came from
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

   ! Notes:
   !  * If LAMMPS returns a scalar that is allocated by the library interface
   !     (see library.cpp), then that memory is deallocated automatically and
   !     the argument to lammps_extract_fix must be a SCALAR.
   !  * If LAMMPS returns a pointer to an array, consisting of internal LAMMPS
   !     data, then the argument must be an interoperable Fortran pointer.
   !     Interoperable means it is of type INTEGER (C_INT) or of type
   !     REAL (C_DOUBLE) in this context.
   !  * Pointers should NEVER be deallocated, as that would deallocate internal
   !     LAMMPS data!
   !  * Note that just because you can read the values of, say, a compute at
   !     any time does not mean those values represent the "correct" values.
   !     LAMMPS will abort you if you try to grab a pointer to a non-current
   !     entity, but once it's bound, it's your responsibility to check that
   !     it's current before evaluating.
   !  * IMPORTANT:  Two-dimensional arrays (such as 'x' from extract_atom)
   !     will be transposed from what they might look like in C++.  This is
   !     because of different bookkeeping conventions between Fortran and C
   !     that date back to about 1970 or so (when C was written).
   !  * Arrays start from 1, EXCEPT for mass from extract_atom, which
   !     starts from 0.  This is because the C array actually has a blank
   !     first element (and thus mass[1] corresponds to the mass of type 1)

   type (C_ptr) :: lmp
   real (C_double), pointer :: compute => NULL()
   real (C_double) :: fix, fix2
   real (C_double), dimension(:), pointer :: compute_v => NULL()
   real (C_double), dimension(:,:), pointer :: x => NULL()
   real (C_double), dimension(:), pointer :: mass => NULL()
   integer, dimension(:), allocatable :: types
   double precision, dimension(:), allocatable :: r
   integer :: error, narg, me, nprocs
   character (len=1024) :: command_line

   call MPI_Init (error)
   call MPI_Comm_rank (MPI_COMM_WORLD, me, error)
   call MPI_Comm_size (MPI_COMM_WORLD, nprocs, error)

   ! You are free to pass any string you like to lammps_open or
   ! lammps_open_no_mpi; here is how you pass it the command line
   !call get_command (command_line)
   !call lammps_open (command_line, MPI_COMM_WORLD, lmp)

   ! And here's how to to it with a string constant of your choice
   call lammps_open_no_mpi ('lmp -log log.simple', lmp)

   call lammps_file (lmp, 'in.simple')
   call lammps_command (lmp, 'run 500')

   ! This extracts f_2 as a scalar (the last two arguments can be arbitrary)
   call lammps_extract_fix (fix, lmp, '2', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, 1, 1)
   print *, 'Fix is ', fix

   ! This extracts f_4[1][1] as a scalar
   call lammps_extract_fix (fix2, lmp, '4', LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY, 1, 1)
   print *, 'Fix 2 is ', fix2

   ! This extracts the scalar compute of compute thermo_temp
   call lammps_extract_compute (compute, lmp, 'thermo_temp', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
   print *, 'Compute is ', compute

   ! This extracts the vector compute of compute thermo_temp
   call lammps_extract_compute (compute_v, lmp, 'thermo_temp', LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)
   print *, 'Vector is ', compute_v

   ! This extracts the masses
   call lammps_extract_atom (mass, lmp, 'mass')
   print *, 'Mass is ', mass(1:)

   ! Extracts a pointer to the arrays of positions for all atoms
   call lammps_extract_atom (x, lmp, 'x')
   if ( .not. associated (x) ) print *, 'x is not associated'
   print *, 'x is ', x(:,1)  ! Prints x, y, z for atom 1

   ! Extracts pointer to atom types
   call lammps_gather_atoms (lmp, 'type', 1, types)
   print *, 'types is ', types(1:3)

   ! Allocates an array and assigns all positions to it
   call lammps_gather_atoms (lmp, 'x', 3, r)
   print *, 'natoms = ', int(lammps_get_natoms(lmp))
   print *, 'size(r) = ', size(r)
   print *, 'r is ', r(1:6)

   ! Puts those position data back
   call lammps_scatter_atoms (lmp, 'x', r)

   call lammps_command (lmp, 'run 1')
   print *, 'x is ', x(:,1)  ! Note that the position updates!
   print *, 'Compute is ', compute  ! This did only because "temp" is part of
                                    ! the thermo output; the vector part did
                                    ! not, and won't until we give LAMMPS a
                                    ! thermo output or other command that
                                    ! requires its value

   call lammps_close (lmp)

   call MPI_Finalize (error)

end program simple

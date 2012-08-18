program simple

   use LAMMPS
   implicit none

   type (lammps_instance) :: lmp
   double precision :: compute, fix, fix2
   double precision, dimension(:), allocatable :: compute_v, mass, r
   double precision, dimension(:,:), allocatable :: x
   real, dimension(:,:), allocatable :: x_r

   call lammps_open_no_mpi ('',lmp)
   call lammps_file (lmp, 'in.simple')
   call lammps_command (lmp, 'run 500')

   call lammps_extract_fix (fix, lmp, '2', 0, 1, 1, 1)
   print *, 'Fix is ', fix

   call lammps_extract_fix (fix2, lmp, '4', 0, 2, 1, 1)
   print *, 'Fix 2 is ', fix2

   call lammps_extract_compute (compute, lmp, 'thermo_temp', 0, 0)
   print *, 'Compute is ', compute

   call lammps_extract_compute (compute_v, lmp, 'thermo_temp', 0, 1)
   print *, 'Vector is ', compute_v

   call lammps_extract_atom (mass, lmp, 'mass')
   print *, 'Mass is ', mass

   call lammps_extract_atom (x, lmp, 'x')
   if ( .not. allocated (x) ) print *, 'x is not allocated'
   print *, 'x is ', x(1,:)

   call lammps_extract_atom (x_r, lmp, 'x')
   if ( .not. allocated (x_r) ) print *, 'x is not allocated'
   print *, 'x_r is ', x_r(1,:)

   call lammps_get_coords (lmp, r)
   print *, 'r is ', r(1:3)

   call lammps_close (lmp)

end program simple

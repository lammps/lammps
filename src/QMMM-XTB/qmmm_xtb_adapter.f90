! LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
! https://www.lammps.org/, Sandia National Laboratories
! LAMMPS development team: developers@lammps.org
!
! This file is distributed under the GNU General Public License.
!
! This adapter uses the LGPL-3.0-or-later libxtb Fortran module interface.
! It is an independent implementation of the SCC callback required for a
! sander-compatible periodic QM/MM Hamiltonian.
module lammps_qmmm_xtb_adapter
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use xtb_mctc_accuracy, only : wp
   use xtb_scf, only : scf
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init_environment => init
   use xtb_type_molecule, only : TMolecule, init_molecule => init
   use xtb_type_restart, only : TRestart
   use xtb_type_solvation, only : TSolvation
   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   implicit none
   private

   public :: lammps_qmmm_xtb_create
   public :: lammps_qmmm_xtb_calculate
   public :: lammps_qmmm_xtb_destroy

   ! The shift is V_MM^corr + A_image q.  getEnergy supplies
   ! q.V_MM^corr + 1/2 q.A_image.q, matching the variational energy used by
   ! sander.  Coordinate derivatives are evaluated by LAMMPS after SCC so
   ! addGradient deliberately contributes no force here.
   type, extends(TSolvation) :: periodic_shift_type
      real(wp), allocatable :: mm_shift(:)
      real(wp), allocatable :: image_response(:, :)
   contains
      procedure :: update => periodic_shift_update
      procedure :: addShift => periodic_shift_add
      procedure :: getEnergy => periodic_shift_energy
      procedure :: addGradient => periodic_shift_gradient
   end type periodic_shift_type

   type :: adapter_state_type
      logical :: active = .false.
      logical :: have_result = .false.
      integer :: nqm = 0
      type(TEnvironment) :: env
      type(TMolecule) :: mol
      type(TxTBCalculator) :: calc
      type(TRestart) :: restart
      class(TSolvation), allocatable :: shift
   end type adapter_state_type

   type(adapter_state_type), save :: state
   logical, save :: xtb_runtime_initialized = .false.

   interface
      subroutine mctc_init(progname, ntimer, verbose)
         character(len=*), intent(in) :: progname
         integer, intent(in) :: ntimer
         logical, intent(in) :: verbose
      end subroutine mctc_init
   end interface

contains

subroutine periodic_shift_update(self, env, num, xyz)
   class(periodic_shift_type), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: num(:)
   real(wp), intent(in) :: xyz(:, :)
end subroutine periodic_shift_update

subroutine periodic_shift_add(self, env, qat, qsh, atomicShift, shellShift)
   class(periodic_shift_type), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(inout) :: atomicShift(:)
   real(wp), intent(inout) :: shellShift(:)

   atomicShift = atomicShift + self%mm_shift + matmul(self%image_response, qat)
end subroutine periodic_shift_add

subroutine periodic_shift_energy(self, env, qat, qsh, energy)
   class(periodic_shift_type), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(out) :: energy
   real(wp), allocatable :: image_potential(:)

   image_potential = matmul(self%image_response, qat)
   energy = dot_product(qat, self%mm_shift) + 0.5_wp*dot_product(qat, image_potential)
end subroutine periodic_shift_energy

subroutine periodic_shift_gradient(self, env, num, xyz, qat, qsh, gradient)
   class(periodic_shift_type), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: num(:)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(inout) :: gradient(:, :)
end subroutine periodic_shift_gradient

integer(c_int) function lammps_qmmm_xtb_create(nqm, atomic_numbers, qm_xyz_bohr, method, charge, uhf, &
      accuracy, maxiter, electronic_temperature) bind(C) result(status)
   integer(c_int), value :: nqm, method, charge, uhf, maxiter
   integer(c_int), intent(in) :: atomic_numbers(*)
   real(c_double), intent(in) :: qm_xyz_bohr(*)
   real(c_double), value :: accuracy, electronic_temperature
   integer, allocatable :: numbers(:)
   real(wp), allocatable :: xyz(:, :)
   logical :: failed
   integer :: i

   status = 1_c_int
   call lammps_qmmm_xtb_destroy()
   if (nqm <= 0) return

   ! xTB's parameter reader uses its process-global persistent environment,
   ! even when a separate TEnvironment is passed to the calculator.  The
   ! public xTB drivers and sander initialize that runtime exactly once.
   if (.not.xtb_runtime_initialized) then
      call mctc_init('xtb/lammps', 11, .false.)
      xtb_runtime_initialized = .true.
   end if

   allocate(numbers(nqm))
   allocate(xyz(3, nqm))
   do i = 1, nqm
      numbers(i) = atomic_numbers(i)
      xyz(:, i) = real(qm_xyz_bohr(3*(i-1)+1:3*i), wp)
   end do

   ! Keep recoverable xTB warnings (for example an EEQ-to-SAD fallback) from
   ! being promoted to fatal errors in a long-running MD calculation.
   call init_environment(state%env, .false.)
   call init_molecule(state%mol, numbers, xyz, chrg=real(charge, wp), uhf=int(uhf))
   ! Pass the parameter filename explicitly.  Besides making the selected
   ! Hamiltonian unambiguous, this avoids relying on the compiler ABI for an
   ! omitted optional CHARACTER argument ahead of METHOD and ACCURACY.
   select case(method)
   case(1)
      call newXTBCalculator(state%env, state%mol, state%calc, 'param_gfn1-xtb.txt', &
         method=1, accuracy=real(accuracy, wp))
   case(2)
      call newXTBCalculator(state%env, state%mol, state%calc, 'param_gfn2-xtb.txt', &
         method=2, accuracy=real(accuracy, wp))
   case default
      return
   end select
   state%calc%maxiter = int(maxiter)
   state%calc%etemp = real(electronic_temperature, wp)
   call state%env%check(failed)
   if (failed) then
      call state%env%show('LAMMPS could not initialize the xTB calculator')
      return
   end if

   ! The EEQ wavefunction guess depends on geometry.  Initialize it from the
   ! actual QM coordinates instead of a coincident placeholder geometry.
   call newWavefunction(state%env, state%mol, state%calc, state%restart)
   call state%env%check(failed)
   if (failed) then
      call state%env%show('LAMMPS could not initialize the xTB wavefunction')
      return
   end if

   allocate(periodic_shift_type :: state%shift)
   select type(shift => state%shift)
   type is(periodic_shift_type)
      allocate(shift%mm_shift(nqm), shift%image_response(nqm, nqm))
      shift%mm_shift = 0.0_wp
      shift%image_response = 0.0_wp
   end select

   state%nqm = nqm
   state%active = .true.
   state%have_result = .false.
   status = 0_c_int
end function lammps_qmmm_xtb_create

integer(c_int) function lammps_qmmm_xtb_calculate(nqm, qm_xyz_bohr, npoint, point_xyz_bohr, &
      point_charge, point_atomic_numbers, mm_hardness, mm_shift_hartree, &
      image_response_hartree, energy_hartree, &
      qm_gradient_hartree_bohr, mulliken_charge, point_gradient_hartree_bohr) &
      bind(C) result(status)
   integer(c_int), value :: nqm, npoint
   integer(c_int), intent(in) :: point_atomic_numbers(*)
   real(c_double), intent(in) :: qm_xyz_bohr(*), point_xyz_bohr(*)
   real(c_double), intent(in) :: point_charge(*)
   real(c_double), value :: mm_hardness
   real(c_double), intent(in) :: mm_shift_hartree(*), image_response_hartree(*)
   real(c_double), intent(out) :: energy_hartree
   real(c_double), intent(out) :: qm_gradient_hartree_bohr(*)
   real(c_double), intent(out) :: mulliken_charge(*)
   real(c_double), intent(out) :: point_gradient_hartree_bohr(*)
   type(scc_results) :: results
   real(wp), allocatable :: gradient(:, :)
   real(wp) :: energy, gap
   logical :: failed
   integer :: i, j, index, z

   status = 1_c_int
   energy_hartree = 0.0_c_double
   if (.not.state%active .or. nqm /= state%nqm .or. npoint < 0) return

   do i = 1, nqm
      do j = 1, 3
         state%mol%xyz(j, i) = real(qm_xyz_bohr(3*(i-1)+j), wp)
      end do
   end do

   call state%calc%pcem%allocate(int(npoint))
   do i = 1, npoint
      do j = 1, 3
         state%calc%pcem%xyz(j, i) = real(point_xyz_bohr(3*(i-1)+j), wp)
      end do
      state%calc%pcem%q(i) = real(point_charge(i), wp)
      ! Match sander's mmhardness branches exactly.  A zero atomic number is
      ! an extra point site and therefore uses hydrogen's hardness.
      if (mm_hardness > 1.0e-6_c_double) then
         state%calc%pcem%gam(i) = real(mm_hardness, wp)
      else if (mm_hardness > -1.0e-6_c_double) then
         state%calc%pcem%gam(i) = state%calc%xtbData%coulomb%chemicalHardness(1)
      else
         z = max(1, int(point_atomic_numbers(i)))
         if (z > size(state%calc%xtbData%coulomb%chemicalHardness)) then
            call state%env%error('MM element is outside the selected xTB parameter range', &
               'qmmm_xtb_adapter.f90')
            call state%calc%pcem%deallocate()
            return
         end if
         state%calc%pcem%gam(i) = abs(real(mm_hardness, wp)) * &
            state%calc%xtbData%coulomb%chemicalHardness(z)
      end if
   end do
   state%calc%pcem%grd = 0.0_wp

   select type(shift => state%shift)
   type is(periodic_shift_type)
      do i = 1, nqm
         shift%mm_shift(i) = real(mm_shift_hartree(i), wp)
      end do
      do j = 1, nqm
         do i = 1, nqm
            index = (j-1)*nqm+i
            shift%image_response(i, j) = real(image_response_hartree(index), wp)
         end do
      end do
   end select

   allocate(gradient(3, nqm), source=0.0_wp)
   energy = 0.0_wp
   gap = 0.0_wp
   call state%mol%update
   call scf(state%env, state%mol, state%restart%wfn, state%calc%basis, state%calc%pcem, &
      state%calc%xtbData, state%shift, gap, state%calc%etemp, state%calc%maxiter, 0, &
      state%have_result, .true., state%calc%accuracy, energy, gradient, results)

   call state%env%check(failed)
   if (failed .or. .not.results%converged) then
      call state%env%show('LAMMPS xTB SCC calculation failed')
      call state%calc%pcem%deallocate()
      return
   end if

   energy_hartree = real(energy, c_double)
   do i = 1, nqm
      mulliken_charge(i) = real(state%restart%wfn%q(i), c_double)
      do j = 1, 3
         qm_gradient_hartree_bohr(3*(i-1)+j) = real(gradient(j, i), c_double)
      end do
   end do
   do i = 1, npoint
      do j = 1, 3
         point_gradient_hartree_bohr(3*(i-1)+j) = &
            real(state%calc%pcem%grd(j, i), c_double)
      end do
   end do

   state%have_result = .true.
   call state%calc%pcem%deallocate()
   status = 0_c_int
end function lammps_qmmm_xtb_calculate

subroutine lammps_qmmm_xtb_destroy() bind(C)
   if (allocated(state%shift)) deallocate(state%shift)
   call state%calc%pcem%deallocate()
   if (allocated(state%calc%basis)) then
      call state%calc%basis%deallocate()
      deallocate(state%calc%basis)
   end if
   if (allocated(state%calc%xtbData)) deallocate(state%calc%xtbData)
   if (allocated(state%restart%wfn%q)) call state%restart%wfn%deallocate()
   if (allocated(state%mol%xyz)) call state%mol%deallocate()
   state%active = .false.
   state%have_result = .false.
   state%nqm = 0
end subroutine lammps_qmmm_xtb_destroy

end module lammps_qmmm_xtb_adapter

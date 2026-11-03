// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_ke_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeKEKokkos<DeviceType>::ComputeKEKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeKE(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeKEKokkos<DeviceType>::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  atomKK->sync(execution_space,datamask_read);

  auto d_v = atomKK->k_v.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  const int nlocal = atom->nlocal;
  const int groupbit_kk = groupbit;

  double ke = 0.0;

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal),
      KOKKOS_LAMBDA(const int i, double &l_ke) {
        if (d_mask(i) & groupbit_kk)
          l_ke += static_cast<double>(d_rmass(i)) *
            static_cast<double>(d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2));
      }, ke);

  } else {

    atomKK->k_mass.template sync<DeviceType>();
    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal),
      KOKKOS_LAMBDA(const int i, double &l_ke) {
        if (d_mask(i) & groupbit_kk)
          l_ke += static_cast<double>(d_mass(d_type(i))) *
            static_cast<double>(d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2));
      }, ke);

  }

  MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeKEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeKEKokkos<LMPHostType>;
#endif
}

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

#include "compute_uf3_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeUF3Kokkos<DeviceType>::ComputeUF3Kokkos(LAMMPS *lmp, int narg, char **arg) :
    ComputeUF3(lmp, narg, arg), atomKK(nullptr)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | TYPE_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> ComputeUF3Kokkos<DeviceType>::~ComputeUF3Kokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void ComputeUF3Kokkos<DeviceType>::init()
{
  ComputeUF3::init();

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void ComputeUF3Kokkos<DeviceType>::compute_array()
{
  atomKK->sync(execution_space, datamask_read);

  if (execution_space != HostKK)
    error->all(FLERR,
                "compute uf3/kk: only Kokkos host execution is implemented; use compute ... uf3/kk/host");

  ComputeUF3::compute_array();
}

namespace LAMMPS_NS {
template class ComputeUF3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeUF3Kokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS

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

#include "fix_store_force_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "memory_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixStoreForceKokkos<DeviceType>::FixStoreForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixStoreForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = F_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;

  // the base constructor allocated foriginal with Memory; redo it as a
  // dual view so the kernel can fill it on the device

  memory->destroy(foriginal);
  foriginal = nullptr;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixStoreForceKokkos<DeviceType>::~FixStoreForceKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_foriginal,foriginal);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixStoreForceKokkos<DeviceType>::init()
{
  FixStoreForce::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix store/force/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixStoreForceKokkos<DeviceType>::post_force(int /*vflag*/)
{
  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_foriginal,foriginal);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_foriginal,foriginal,nmax,3,"store/force:foriginal");
    array_atom = foriginal;
    d_foriginal = k_foriginal.template view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);

  d_f = atomKK->k_f.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  const int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixStoreForce>(0,nlocal),*this);
  copymode = 0;

  // array_atom is read on the host by dumps, computes and variables

  k_foriginal.template modify<DeviceType>();
  k_foriginal.sync_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixStoreForceKokkos<DeviceType>::operator()(TagFixStoreForce, const int &i) const
{
  if (d_mask[i] & groupbit) {
    d_foriginal(i,0) = static_cast<KK_FLOAT>(d_f(i,0));
    d_foriginal(i,1) = static_cast<KK_FLOAT>(d_f(i,1));
    d_foriginal(i,2) = static_cast<KK_FLOAT>(d_f(i,2));
  } else {
    d_foriginal(i,0) = d_foriginal(i,1) = d_foriginal(i,2) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixStoreForceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixStoreForceKokkos<LMPHostType>;
#endif
}

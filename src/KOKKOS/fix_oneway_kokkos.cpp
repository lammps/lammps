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

#include "fix_oneway_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "kokkos_base.h"
#include "region.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOneWayKokkos<DeviceType>::FixOneWayKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixOneWay(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | MASK_MASK;
  datamask_modify = V_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOneWayKokkos<DeviceType>::~FixOneWayKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOneWayKokkos<DeviceType>::init()
{
  FixOneWay::init();

  if (!(utils::strmatch(region->style, "^block") || utils::strmatch(region->style, "^sphere")))
    error->all(FLERR, "Cannot (yet) use {}-style region with fix oneway/kk", region->style);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOneWayKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(execution_space, X_MASK | V_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  // build region match array

  region->prematch();
  DAT::tdual_int_1d k_match = DAT::tdual_int_1d("oneway:k_match", atom->nlocal);
  KokkosBase *regionKKBase = dynamic_cast<KokkosBase *>(region);
  regionKKBase->match_all_kokkos(groupbit, k_match);
  k_match.template sync<DeviceType>();
  d_match = k_match.template view<DeviceType>();

  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixOneWay>(0, nlocal), *this);
  copymode = 0;

  atomKK->modified(execution_space, V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixOneWayKokkos<DeviceType>::operator()(TagFixOneWay, const int &i) const
{
  if ((mask[i] & groupbit) && d_match[i]) {
    // bits 0-1 = coordinate index (0=x, 1=y, 2=z), bit 2 = minus direction
    const int idx = direction & 3;
    if (direction & 4) {
      if (v(i,idx) > 0.0) v(i,idx) = -v(i,idx);
    } else {
      if (v(i,idx) < 0.0) v(i,idx) = -v(i,idx);
    }
  }
}

namespace LAMMPS_NS {
template class FixOneWayKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixOneWayKokkos<LMPHostType>;
#endif
}

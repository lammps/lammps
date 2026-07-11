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

#include "fix_nve_noforce_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVENoforceKokkos<DeviceType>::FixNVENoforceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVENoforce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | MASK_MASK;
  datamask_modify = X_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVENoforceKokkos<DeviceType>::~FixNVENoforceKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVENoforceKokkos<DeviceType>::init()
{
  FixNVENoforce::init();

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Cannot (yet) use respa with fix nve/noforce/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVENoforceKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | V_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVENoforce>(0, nlocal), *this);
  copymode = 0;

  atomKK->modified(execution_space, X_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNVENoforceKokkos<DeviceType>::operator()(TagFixNVENoforce, const int &i) const
{
  if (mask[i] & groupbit) {
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

namespace LAMMPS_NS {
template class FixNVENoforceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVENoforceKokkos<LMPHostType>;
#endif
}

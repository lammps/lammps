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

#include "fix_efield_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEfieldKokkos<DeviceType>::FixEfieldKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEfield(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | Q_MASK | F_MASK | RMASS_MASK | MASK_MASK | TYPE_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEfieldKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // update efield due to variables

  update_efield_variables();

  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;

  eflag = 0;

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), *this);
  
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::operator()(const int i) const
{
  if (mask[i] & groupbit) {
    double qi = q[i];
    f(i,0) += qi*ex;
    f(i,1) += qi*ey;
    f(i,2) += qi*ez;
  }
}

namespace LAMMPS_NS {
template class FixEfieldKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixEfieldKokkos<LMPHostType>;
#endif
}

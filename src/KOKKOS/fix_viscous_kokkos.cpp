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

#include "fix_viscous_kokkos.h"

#include "atom_kokkos.h"
#include "update.h"
#include "modify.h"
#include "input.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "kokkos_base.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixViscousKokkos<DeviceType>::FixViscousKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixViscous(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixViscousKokkos<DeviceType>::~FixViscousKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixViscousKokkos<DeviceType>::init()
{
  FixViscous::init();

  k_gamma = Kokkos::DualView<double*, Kokkos::LayoutRight, DeviceType>("FixViscousKokkos:gamma",atom->ntypes+1);

  for (int i = 1; i <= atom->ntypes; i++) k_gamma.h_view(i) = gamma[i];

  k_gamma.template modify<LMPHostType>();
  k_gamma.template sync<DeviceType>();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixViscousKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, V_MASK | F_MASK | MASK_MASK | TYPE_MASK);

  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixViscous>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixViscousKokkos<DeviceType>::operator()(TagFixViscous, const int &i) const {
  if (mask[i] & groupbit) {
    double drag = k_gamma.d_view(type[i]);
    f(i,0) -= drag*v(i,0);
    f(i,1) -= drag*v(i,1);
    f(i,2) -= drag*v(i,2);
  }
}

namespace LAMMPS_NS {
template class FixViscousKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixViscousKokkos<LMPHostType>;
#endif
}


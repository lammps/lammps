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

#include "fix_viscous_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// type of scaling (mirrors enum in fix_viscous_sphere.cpp)
enum { VS_NONE, VS_TYPE, VS_VARIABLE };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixViscousSphereKokkos<DeviceType>::FixViscousSphereKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixViscousSphere(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixViscousSphereKokkos<DeviceType>::~FixViscousSphereKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixViscousSphereKokkos<DeviceType>::init()
{
  FixViscousSphere::init();

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Cannot (yet) use respa with fix viscous/sphere/kk");

  if (scalestyle == VS_VARIABLE)
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot (yet) use variable style with fix viscous/sphere/kk");

  if (scalestyle == VS_TYPE) {
    k_eff_gamma = Kokkos::DualView<KK_FLOAT*, Kokkos::LayoutRight, DeviceType>(
        "FixViscousSphereKokkos:eff_gamma", atom->ntypes + 1);
    for (int i = 1; i <= atom->ntypes; i++)
      k_eff_gamma.view_host()(i) = (KK_FLOAT)(gamma * scalegamma[i]);
    k_eff_gamma.modify_host();
    k_eff_gamma.template sync<DeviceType>();
  } else {
    m_gamma = (KK_FLOAT) gamma;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixViscousSphereKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, OMEGA_MASK | TORQUE_MASK | MASK_MASK | TYPE_MASK);

  omega = atomKK->k_omega.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixViscousSphere>(0, nlocal), *this);
  copymode = 0;

  atomKK->modified(execution_space, TORQUE_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixViscousSphereKokkos<DeviceType>::operator()(TagFixViscousSphere, const int &i) const
{
  if (mask[i] & groupbit) {
    KK_FLOAT drag;
    if (scalestyle == VS_TYPE)
      drag = k_eff_gamma.view_device()(type[i]);
    else
      drag = m_gamma;
    torque(i,0) -= drag * omega(i,0);
    torque(i,1) -= drag * omega(i,1);
    torque(i,2) -= drag * omega(i,2);
  }
}

namespace LAMMPS_NS {
template class FixViscousSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixViscousSphereKokkos<LMPHostType>;
#endif
}

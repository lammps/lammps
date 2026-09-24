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

#include "fix_viscous_nonlinear_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "math_const.h"
#include "update.h"

using namespace LAMMPS_NS;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixViscousNonlinearKokkos<DeviceType>::FixViscousNonlinearKokkos(LAMMPS *lmp, int narg,
                                                                char **arg) :
  FixViscousNonlinear(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = V_MASK | F_MASK | RADIUS_MASK | MASK_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixViscousNonlinearKokkos<DeviceType>::init()
{
  FixViscousNonlinear::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix viscous/nonlinear/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixViscousNonlinearKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space,datamask_read);

  d_f = atomKK->k_f.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  for (int i = 0; i < 3; i++) m_v_fluid[i] = static_cast<KK_FLOAT>(v_fluid[i]);

  const int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixViscousNonlinear>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixViscousNonlinearKokkos<DeviceType>::operator()(TagFixViscousNonlinear,
                                                       const int &i) const
{
  if (d_mask[i] & groupbit) {

    // apply Schiller-Naumann drag relative to the (uniform) fluid velocity

    const KK_FLOAT vrel0 = d_v(i,0) - m_v_fluid[0];
    const KK_FLOAT vrel1 = d_v(i,1) - m_v_fluid[1];
    const KK_FLOAT vrel2 = d_v(i,2) - m_v_fluid[2];
    const KK_FLOAT vmag = Kokkos::sqrt(vrel0*vrel0 + vrel1*vrel1 + vrel2*vrel2);
    if (vmag == static_cast<KK_FLOAT>(0.0)) return;

    const KK_FLOAT rho_fluid_kk = static_cast<KK_FLOAT>(rho_fluid);
    const KK_FLOAT mu_fluid_kk = static_cast<KK_FLOAT>(mu_fluid);

    const KK_FLOAT r = d_radius[i];
    const KK_FLOAT re = rho_fluid_kk * vmag * (static_cast<KK_FLOAT>(2.0) * r) / mu_fluid_kk;
    const KK_FLOAT cd = (static_cast<KK_FLOAT>(24.0) / re) * (static_cast<KK_FLOAT>(1.0) +
      static_cast<KK_FLOAT>(0.15) * Kokkos::pow(re, static_cast<KK_FLOAT>(0.687)));

    // F = -1/2 Cd rho_g (pi r^2) |v_rel| v_rel

    const KK_FLOAT pref = static_cast<KK_FLOAT>(0.5) * cd * rho_fluid_kk *
      static_cast<KK_FLOAT>(MY_PI) * r * r * vmag;
    d_f(i,0) -= static_cast<KK_ACC_FLOAT>(pref * vrel0);
    d_f(i,1) -= static_cast<KK_ACC_FLOAT>(pref * vrel1);
    d_f(i,2) -= static_cast<KK_ACC_FLOAT>(pref * vrel2);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixViscousNonlinearKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixViscousNonlinearKokkos<LMPHostType>;
#endif
}

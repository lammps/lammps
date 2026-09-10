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

#include "fix_spring_rg_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "group.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringRGKokkos<DeviceType>::FixSpringRGKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSpringRG(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringRGKokkos<DeviceType>::init()
{
  // the base class init() computes the group mass and, on the first run
  // with the "rg0 NULL" option, the reference radius of gyration on the host

  atomKK->sync(Host, X_MASK | IMAGE_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  FixSpringRG::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix spring/rg/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringRGKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  l_rmass_flag = (atom->rmass != nullptr) ? 1 : 0;
  if (l_rmass_flag) {
    rmass = atomKK->k_rmass.view<DeviceType>();
  } else {
    atomKK->k_mass.sync<DeviceType>();
    mass = atomKK->k_mass.view<DeviceType>();
    type = atomKK->k_type.view<DeviceType>();
  }
  const int nlocal = atom->nlocal;

  prd = Few<double,3>(domain->prd);
  h = Few<double,6>(domain->h);
  triclinic = domain->triclinic;

  // total mass and center of mass of the group, reduced on the device
  // from unwrapped coordinates (mirrors Group::mass() + Group::xcm())

  double mlocal[4], mall[4];
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixSpringRGXcm>(0,nlocal),*this,mlocal);
  copymode = 0;
  MPI_Allreduce(mlocal,mall,4,MPI_DOUBLE,MPI_SUM,world);

  if (group->dynamic[igroup]) masstotal = mall[3];

  l_xcm[0] = l_xcm[1] = l_xcm[2] = 0.0;
  if (masstotal > 0.0) {
    l_xcm[0] = mall[0]/masstotal;
    l_xcm[1] = mall[1]/masstotal;
    l_xcm[2] = mall[2]/masstotal;
  }

  // radius of gyration about the center of mass, reduced on the device
  // (mirrors Group::gyration(); only slot 0 of the reduction is used)

  double rglocal[4], rgall[4];
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixSpringRGGyration>(0,nlocal),*this,rglocal);
  copymode = 0;
  MPI_Allreduce(rglocal,rgall,4,MPI_DOUBLE,MPI_SUM,world);

  double rg = 0.0;
  if (masstotal > 0.0) rg = sqrt(rgall[0]/masstotal);

  // rg == 0 means that either there are no atoms in the group or that
  //         they are exactly on top of each other. nothing to do then

  if ((rg == 0.0) || (masstotal == 0.0)) return;

  // apply restoring forces to atoms on device

  l_coeff = 2.0 * k * (1.0 - rg0/rg) / masstotal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringRGApply>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringRGKokkos<DeviceType>::operator()(TagFixSpringRGXcm, const int &i, value_type result) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double massone = l_rmass_flag ? (double)rmass[i] : (double)mass[type[i]];
    result[0] += massone * unwrap[0];
    result[1] += massone * unwrap[1];
    result[2] += massone * unwrap[2];
    result[3] += massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringRGKokkos<DeviceType>::operator()(TagFixSpringRGGyration, const int &i, value_type result) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double dx = unwrap[0] - l_xcm[0];
    const double dy = unwrap[1] - l_xcm[1];
    const double dz = unwrap[2] - l_xcm[2];
    const double massone = l_rmass_flag ? (double)rmass[i] : (double)mass[type[i]];
    result[0] += massone * (dx*dx + dy*dy + dz*dz);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringRGKokkos<DeviceType>::operator()(TagFixSpringRGApply, const int &i) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double massone = l_rmass_flag ? (double)rmass[i] : (double)mass[type[i]];
    f(i,0) -= static_cast<KK_ACC_FLOAT>(l_coeff * massone * (unwrap[0] - l_xcm[0]));
    f(i,1) -= static_cast<KK_ACC_FLOAT>(l_coeff * massone * (unwrap[1] - l_xcm[1]));
    f(i,2) -= static_cast<KK_ACC_FLOAT>(l_coeff * massone * (unwrap[2] - l_xcm[2]));
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixSpringRGKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSpringRGKokkos<LMPHostType>;
#endif
}

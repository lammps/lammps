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
  FixSpringRG::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix spring/rg/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringRGKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // compute center of mass and radius of gyration on host

  atomKK->sync(Host, X_MASK | IMAGE_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  double xcm[3];
  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);
  group->xcm(igroup, masstotal, xcm);
  double rg = group->gyration(igroup, masstotal, xcm);

  double term1 = 2.0 * k * (1.0 - rg0/rg);

  // apply restoring forces to atoms on device

  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  l_xcm[0] = xcm[0];
  l_xcm[1] = xcm[1];
  l_xcm[2] = xcm[2];
  l_coeff = (masstotal > 0.0) ? term1 / masstotal : 0.0;

  prd = Few<double,3>(domain->prd);
  h = Few<double,6>(domain->h);
  triclinic = domain->triclinic;

  copymode = 1;
  if (atom->rmass) {
    rmass = atomKK->k_rmass.view<DeviceType>();
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringRGRmass>(0,nlocal),*this);
  } else {
    mass = atomKK->k_mass.view<DeviceType>();
    type = atomKK->k_type.view<DeviceType>();
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringRG>(0,nlocal),*this);
  }
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringRGKokkos<DeviceType>::operator()(TagFixSpringRG, const int &i) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double massone = mass[type[i]];
    f(i,0) -= l_coeff * massone * (unwrap[0] - l_xcm[0]);
    f(i,1) -= l_coeff * massone * (unwrap[1] - l_xcm[1]);
    f(i,2) -= l_coeff * massone * (unwrap[2] - l_xcm[2]);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringRGKokkos<DeviceType>::operator()(TagFixSpringRGRmass, const int &i) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double massone = rmass[i];
    f(i,0) -= l_coeff * massone * (unwrap[0] - l_xcm[0]);
    f(i,1) -= l_coeff * massone * (unwrap[1] - l_xcm[1]);
    f(i,2) -= l_coeff * massone * (unwrap[2] - l_xcm[2]);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixSpringRGKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSpringRGKokkos<LMPHostType>;
#endif
}

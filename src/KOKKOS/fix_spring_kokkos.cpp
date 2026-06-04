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

#include "fix_spring_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "group.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double SMALL = 1.0e-10;
enum{TETHER,COUPLE};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringKokkos<DeviceType>::FixSpringKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSpring(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringKokkos<DeviceType>::init()
{
  FixSpring::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix spring/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringKokkos<DeviceType>::post_force(int /*vflag*/)
{
  double xcm[3], xcm2[3];
  double dx, dy, dz, fx, fy, fz, r, dr;

  if (styleflag == TETHER) {

    // compute center of mass of group on host

    atomKK->sync(Host, X_MASK | MASK_MASK);

    if (group->dynamic[igroup]) masstotal = group->mass(igroup);
    group->xcm(igroup, masstotal, xcm);

    // compute scalar forces from xcm displacement

    dx = xcm[0] - xc;
    dy = xcm[1] - yc;
    dz = xcm[2] - zc;
    if (!xflag) dx = 0.0;
    if (!yflag) dy = 0.0;
    if (!zflag) dz = 0.0;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    r = MAX(r, SMALL);
    dr = r - r0;

    fx = k_spring*dx*dr/r;
    fy = k_spring*dy*dr/r;
    fz = k_spring*dz*dr/r;
    ftotal[0] = -fx;
    ftotal[1] = -fy;
    ftotal[2] = -fz;
    ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
    if (dr < 0.0) ftotal[3] = -ftotal[3];
    espring = 0.5*k_spring * dr*dr;

    if (masstotal > 0.0) {
      l_fx = fx / masstotal;
      l_fy = fy / masstotal;
      l_fz = fz / masstotal;
    } else {
      l_fx = l_fy = l_fz = 0.0;
    }

    // apply forces to atoms on device

    atomKK->sync(execution_space, F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

    f = atomKK->k_f.view<DeviceType>();
    mask = atomKK->k_mask.view<DeviceType>();
    int nlocal = atom->nlocal;

    copymode = 1;
    if (atom->rmass) {
      rmass = atomKK->k_rmass.view<DeviceType>();
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringTetherRmass>(0,nlocal),*this);
    } else {
      mass = atomKK->k_mass.view<DeviceType>();
      type = atomKK->k_type.view<DeviceType>();
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringTether>(0,nlocal),*this);
    }
    copymode = 0;

  } else {    // COUPLE

    // compute both centers of mass on host

    atomKK->sync(Host, X_MASK | MASK_MASK);

    if (group->dynamic[igroup]) masstotal = group->mass(igroup);
    if (group->dynamic[igroup2]) masstotal2 = group->mass(igroup2);
    group->xcm(igroup, masstotal, xcm);
    group->xcm(igroup2, masstotal2, xcm2);

    // compute scalar forces from xcm displacement

    dx = xcm2[0] - xcm[0] - xc;
    dy = xcm2[1] - xcm[1] - yc;
    dz = xcm2[2] - xcm[2] - zc;
    if (!xflag) dx = 0.0;
    if (!yflag) dy = 0.0;
    if (!zflag) dz = 0.0;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    r = MAX(r, SMALL);
    dr = r - r0;

    fx = k_spring*dx*dr/r;
    fy = k_spring*dy*dr/r;
    fz = k_spring*dz*dr/r;
    ftotal[0] = fx;
    ftotal[1] = fy;
    ftotal[2] = fz;
    ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
    if (dr < 0.0) ftotal[3] = -ftotal[3];
    espring = 0.5*k_spring * dr*dr;

    if (masstotal2 > 0.0) {
      l_fx2 = fx / masstotal2;
      l_fy2 = fy / masstotal2;
      l_fz2 = fz / masstotal2;
    } else {
      l_fx2 = l_fy2 = l_fz2 = 0.0;
    }

    if (masstotal > 0.0) {
      l_fx = fx / masstotal;
      l_fy = fy / masstotal;
      l_fz = fz / masstotal;
    } else {
      l_fx = l_fy = l_fz = 0.0;
    }

    l_group2bit = group2bit;

    // apply forces to atoms on device

    atomKK->sync(execution_space, F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

    f = atomKK->k_f.view<DeviceType>();
    mask = atomKK->k_mask.view<DeviceType>();
    int nlocal = atom->nlocal;

    copymode = 1;
    if (atom->rmass) {
      rmass = atomKK->k_rmass.view<DeviceType>();
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringCoupleRmass>(0,nlocal),*this);
    } else {
      mass = atomKK->k_mass.view<DeviceType>();
      type = atomKK->k_type.view<DeviceType>();
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringCouple>(0,nlocal),*this);
    }
    copymode = 0;
  }

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringKokkos<DeviceType>::operator()(TagFixSpringTether, const int &i) const
{
  if (mask[i] & groupbit) {
    const double massone = mass[type[i]];
    f(i,0) -= l_fx * massone;
    f(i,1) -= l_fy * massone;
    f(i,2) -= l_fz * massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringKokkos<DeviceType>::operator()(TagFixSpringTetherRmass, const int &i) const
{
  if (mask[i] & groupbit) {
    const double massone = rmass[i];
    f(i,0) -= l_fx * massone;
    f(i,1) -= l_fy * massone;
    f(i,2) -= l_fz * massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringKokkos<DeviceType>::operator()(TagFixSpringCouple, const int &i) const
{
  const double massone = mass[type[i]];
  if (mask[i] & groupbit) {
    f(i,0) += l_fx * massone;
    f(i,1) += l_fy * massone;
    f(i,2) += l_fz * massone;
  }
  if (mask[i] & l_group2bit) {
    f(i,0) -= l_fx2 * massone;
    f(i,1) -= l_fy2 * massone;
    f(i,2) -= l_fz2 * massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringKokkos<DeviceType>::operator()(TagFixSpringCoupleRmass, const int &i) const
{
  const double massone = rmass[i];
  if (mask[i] & groupbit) {
    f(i,0) += l_fx * massone;
    f(i,1) += l_fy * massone;
    f(i,2) += l_fz * massone;
  }
  if (mask[i] & l_group2bit) {
    f(i,0) -= l_fx2 * massone;
    f(i,1) -= l_fy2 * massone;
    f(i,2) -= l_fz2 * massone;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixSpringKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSpringKokkos<LMPHostType>;
#endif
}

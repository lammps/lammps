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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy, alphataubio at gmail
------------------------------------------------------------------------- */

#include "fix_nve_limit_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "kokkos_type.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVELimitKokkos<DeviceType>::FixNVELimitKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVELimit(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK = (AtomKokkos *) atom;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNVELimitKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_v = atomKK->k_v.template view<DeviceType>();
  auto d_f = atomKK->k_f.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto l_groupbit = groupbit;
  auto l_dtf = dtf;
  auto l_dtv = dtv;
  auto l_vlimitsq = vlimitsq;

  int d_ncount;

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space, X_MASK|V_MASK|F_MASK|MASK_MASK|RMASS_MASK );

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), KOKKOS_LAMBDA(const int i, int &l_ncount) {
      if (d_mask[i] & l_groupbit) {
        const double dtfm = l_dtf / d_rmass[i];
        d_v(i,0) += dtfm * d_f(i,0);
        d_v(i,1) += dtfm * d_f(i,1);
        d_v(i,2) += dtfm * d_f(i,2);

        const double vsq = d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2);
        if (vsq > l_vlimitsq) {
          l_ncount++;
          const double scale = sqrt(l_vlimitsq/vsq);
          d_v(i,0) *= scale;
          d_v(i,1) *= scale;
          d_v(i,2) *= scale;
        }

        d_x(i,0) += l_dtv * d_v(i,0);
        d_x(i,1) += l_dtv * d_v(i,1);
        d_x(i,2) += l_dtv * d_v(i,2);
      }
    }, d_ncount);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    auto l_groupbit = groupbit;
    atomKK->sync(execution_space, X_MASK|V_MASK|F_MASK|MASK_MASK|TYPE_MASK );

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), KOKKOS_LAMBDA(const int i, int &l_ncount) {
      if (d_mask[i] & l_groupbit) {
        const double dtfm = l_dtf / d_mass[d_type[i]];
        d_v(i,0) += dtfm * d_f(i,0);
        d_v(i,1) += dtfm * d_f(i,1);
        d_v(i,2) += dtfm * d_f(i,2);

        const double vsq = d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2);
        if (vsq > l_vlimitsq) {
          l_ncount++;
          const double scale = sqrt(l_vlimitsq/vsq);
          d_v(i,0) *= scale;
          d_v(i,1) *= scale;
          d_v(i,2) *= scale;
        }

        d_x(i,0) += l_dtv * d_v(i,0);
        d_x(i,1) += l_dtv * d_v(i,1);
        d_x(i,2) += l_dtv * d_v(i,2);
      }
    }, d_ncount);
  }

  ncount += d_ncount;
  atomKK->modified(execution_space, X_MASK | V_MASK );
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVELimitKokkos<DeviceType>::final_integrate()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  auto d_v = atomKK->k_v.template view<DeviceType>();
  auto d_f = atomKK->k_f.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto l_groupbit = groupbit;
  auto l_dtf = dtf;
  auto l_vlimitsq = vlimitsq;

  int d_ncount;

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    atomKK->sync(execution_space, V_MASK|F_MASK|MASK_MASK|RMASS_MASK );

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), KOKKOS_LAMBDA(const int i, int &l_ncount) {
      if (d_mask[i] & l_groupbit) {
        const double dtfm = l_dtf / d_rmass[i];
        d_v(i,0) += dtfm * d_f(i,0);
        d_v(i,1) += dtfm * d_f(i,1);
        d_v(i,2) += dtfm * d_f(i,2);

        const double vsq = d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2);
        if (vsq > l_vlimitsq) {
          l_ncount++;
          const double scale = sqrt(l_vlimitsq/vsq);
          d_v(i,0) *= scale;
          d_v(i,1) *= scale;
          d_v(i,2) *= scale;
        }
      }
    }, d_ncount);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space, V_MASK|F_MASK|MASK_MASK|TYPE_MASK );

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), KOKKOS_LAMBDA(const int i, int &l_ncount) {
      if (d_mask[i] & l_groupbit) {
        const double dtfm = l_dtf / d_mass[d_type[i]];
        d_v(i,0) += dtfm * d_f(i,0);
        d_v(i,1) += dtfm * d_f(i,1);
        d_v(i,2) += dtfm * d_f(i,2);

        const double vsq = d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2);
        if (vsq > l_vlimitsq) {
          l_ncount++;
          const double scale = sqrt(l_vlimitsq/vsq);
          d_v(i,0) *= scale;
          d_v(i,1) *= scale;
          d_v(i,2) *= scale;
        }
      }
    }, d_ncount);
  }

  ncount += d_ncount;
  atomKK->modified(execution_space, V_MASK );
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNVELimitKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVELimitKokkos<LMPHostType>;
#endif
}


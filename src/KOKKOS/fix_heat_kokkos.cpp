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

#include "fix_heat_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "kokkos_base.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixHeatKokkos<DeviceType>::FixHeatKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixHeat(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixHeatKokkos<DeviceType>::init()
{
  // the base class init() sums the group mass on the host

  atomKK->sync(Host, MASK_MASK | TYPE_MASK | RMASS_MASK);

  FixHeat::init();

  if (region && !dynamic_cast<KokkosBase*>(region))
    error->all(FLERR, "Cannot use fix heat/kk with region style {} that has no KOKKOS support",
               region->style);

  // evaluating a per-atom heat variable is a host-only concept; supporting it
  // would force per-atom host/device transfers every heat step

  if (hstyle == ATOM)
    error->all(FLERR, "Fix heat/kk does not (yet) support an atom-style variable as heat flux");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixHeatKokkos<DeviceType>::end_of_step()
{
  // evaluate equal-style variable on the host

  if (hstyle == EQUAL) {
    modify->clearstep_compute();
    heat_input = input->variable->compute_equal(hvar);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  atomKK->sync(execution_space, V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);

  v = atomKK->k_v.view<DeviceType>();
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

  // with a region, additionally restrict all sums to atoms currently inside
  // it, and use the mass of those atoms in place of the static group mass

  DAT::tdual_int_1d k_match;
  l_region_flag = (region != nullptr) ? 1 : 0;
  if (l_region_flag) {
    region->prematch();
    KokkosBase *regionKKBase = dynamic_cast<KokkosBase*>(region);
    k_match = DAT::tdual_int_1d("heat:k_match", nlocal);
    regionKKBase->match_all_kokkos(groupbit, k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // group momentum, twice the unscaled kinetic energy sum, and group mass,
  // reduced on the device (mirrors Group::vcm() + Group::ke() + Group::mass())

  double result[5], all[5];
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixHeatKE>(0,nlocal),*this,result);
  copymode = 0;
  MPI_Allreduce(result,all,5,MPI_DOUBLE,MPI_SUM,world);

  if (l_region_flag) {
    masstotal = all[4];
    if (masstotal == 0.0) error->all(FLERR, "Fix heat group has no atoms");
  }

  double vcm[3];
  vcm[0] = all[0]/masstotal;
  vcm[1] = all[1]/masstotal;
  vcm[2] = all[2]/masstotal;
  const double ke = 0.5 * force->mvv2e * all[3] * force->ftm2v;
  const double vcmsq = vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2];

  // add heat via scale factor on velocities
  // scale = velocity scale factor to accomplish eflux change in energy
  // vsub = velocity subtracted from each atom to preserve momentum
  // overall KE cannot go negative

  const double heat = heat_input * nevery * update->dt * force->ftm2v;
  const double escale =
    (ke + heat - 0.5*vcmsq*masstotal) / (ke - 0.5*vcmsq*masstotal);
  if (escale < 0.0) error->all(FLERR, "Fix heat kinetic energy went negative");
  scale = sqrt(escale);

  l_scale = static_cast<KK_FLOAT>(scale);
  l_vsub[0] = static_cast<KK_FLOAT>((scale - 1.0)*vcm[0]);
  l_vsub[1] = static_cast<KK_FLOAT>((scale - 1.0)*vcm[1]);
  l_vsub[2] = static_cast<KK_FLOAT>((scale - 1.0)*vcm[2]);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixHeatApply>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space, V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixHeatKokkos<DeviceType>::operator()(TagFixHeatKE, const int &i, value_type result) const
{
  if ((mask[i] & groupbit) && (!l_region_flag || d_match[i])) {
    const KK_FLOAT massone = l_rmass_flag ? rmass[i] : mass[type[i]];
    const KK_FLOAT v0 = v(i,0);
    const KK_FLOAT v1 = v(i,1);
    const KK_FLOAT v2 = v(i,2);
    result[0] += static_cast<double>(massone * v0);
    result[1] += static_cast<double>(massone * v1);
    result[2] += static_cast<double>(massone * v2);
    result[3] += static_cast<double>(massone * (v0*v0 + v1*v1 + v2*v2));
    result[4] += static_cast<double>(massone);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixHeatKokkos<DeviceType>::operator()(TagFixHeatApply, const int &i) const
{
  if ((mask[i] & groupbit) && (!l_region_flag || d_match[i])) {
    v(i,0) = l_scale * v(i,0) - l_vsub[0];
    v(i,1) = l_scale * v(i,1) - l_vsub[1];
    v(i,2) = l_scale * v(i,2) - l_vsub[2];
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixHeatKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixHeatKokkos<LMPHostType>;
#endif
}

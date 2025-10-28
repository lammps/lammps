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

#include "fix_electron_stopping_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos_base.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "region.h"
#include "update.h"

#include <climits>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixElectronStoppingKokkos<DeviceType>::FixElectronStoppingKokkos(LAMMPS *lmp, int narg,
                                                                 char **arg) :
    FixElectronStopping(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | TYPE_MASK | TAG_MASK | MASK_MASK | RMASS_MASK;
  datamask_modify = F_MASK;

  const int ncol = atomKK->ntypes + 1;
  typename AT::tdual_double_2d k_elstop_ranges("k_elstop_ranges", ncol, table_entries);
  typename AT::tdual_double_2d::t_host h_elstop_ranges = k_elstop_ranges.view_host();
  for (int r = 0; r < table_entries; r++) {
    for (int c = 0; c < ncol; c++) { h_elstop_ranges(c, r) = elstop_ranges[c][r]; }
  }
  k_elstop_ranges.modify_host();
  k_elstop_ranges.template sync<DeviceType>();
  d_elstop_ranges = k_elstop_ranges.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void FixElectronStoppingKokkos<DeviceType>::init()
{
  FixElectronStopping::init();
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void FixElectronStoppingKokkos<DeviceType>::post_force(int /*vflag*/)
{
  SeLoss_sync_flag = 0;

  atomKK->sync(execution_space, datamask_read);
  neighbor->build_one(list);

  NeighListKokkos<DeviceType> *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  d_numneigh = k_list->d_numneigh;

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();

  int nlocal = atom->nlocal;

  dt = update->dt;
  mvv2e = force->mvv2e;

  // update region if necessary

  if (region) {
    if (!(utils::strmatch(region->style, "^block") || utils::strmatch(region->style, "^sphere")))
      error->all(FLERR, "Cannot (yet) use {}-style region with fix electron/stopping/kk",
                 region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("electron_stopping:k_match", nlocal);
    KokkosBase *regionKKBase = dynamic_cast<KokkosBase *>(region);
    regionKKBase->match_all_kokkos(groupbit, k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  copymode = 1;
  double seloss = SeLoss;

  // Min reduction used to print error message for the smallest local
  // index. Matches the behavior of the non-Kokkos version of this
  // fix.
  FixElectronStoppingErrorValue error_value = {nlocal, 0};

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixElectronStopping>(0, nlocal), *this,
                          Kokkos::Sum<double>(seloss),
                          Kokkos::Min<FixElectronStoppingErrorValue>(error_value));
  if (error_value.i < nlocal) {
    atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space, TAG_MASK);
    tagint error_tag = (atomKK->k_tag.view<LMPHostType>())(error_value.i);
    error->one(FLERR, "Fix electron/stopping/kk: kinetic energy too high for atom {}: {} vs {}",
               error_tag, error_value.energy, elstop_ranges[0][table_entries - 1]);
  }
  SeLoss += seloss;
  copymode = 0;

  atomKK->modified(execution_space, datamask_modify);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void
FixElectronStoppingKokkos<DeviceType>::operator()(TagFixElectronStopping, const int &i,
                                                  double &seloss,
                                                  FixElectronStoppingErrorValue &error_value) const
{
  if (d_mask(i) & groupbit) {
    if (d_numneigh(i) < minneigh) return;
    if (region && !d_match(i)) return;

    int itype = type(i);
    double massone = (d_rmass.data()) ? d_rmass(i) : d_mass(itype);
    double v2 = v(i, 0) * v(i, 0) + v(i, 1) * v(i, 1) + v(i, 2) * v(i, 2);
    double energy = 0.5 * mvv2e * massone * v2;

    if (energy < Ecut) return;
    if (energy < d_elstop_ranges(0, 0)) return;
    if (energy > d_elstop_ranges(0, table_entries - 1)) {
      error_value = {i, energy};
      return;
    }

    // Binary search to find correct energy range
    int iup = table_entries - 1;
    int idown = 0;
    while (true) {
      int ihalf = idown + (iup - idown) / 2;
      if (ihalf == idown) break;
      if (d_elstop_ranges(0, ihalf) < energy)
        idown = ihalf;
      else
        iup = ihalf;
    }

    double Se_lo = d_elstop_ranges(itype, idown);
    double Se_hi = d_elstop_ranges(itype, iup);
    double E_lo = d_elstop_ranges(0, idown);
    double E_hi = d_elstop_ranges(0, iup);

    // Get electronic stopping with a simple linear interpolation
    double Se = (Se_hi - Se_lo) / (E_hi - E_lo) * (energy - E_lo) + Se_lo;

    double vabs = Kokkos::sqrt(v2);
    double factor = -Se / vabs;

    f(i, 0) += v(i, 0) * factor;
    f(i, 1) += v(i, 1) * factor;
    f(i, 2) += v(i, 2) * factor;

    seloss += Se * vabs * dt;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixElectronStoppingKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixElectronStoppingKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS

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

#include "fix_aveforce_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "input.h"
#include "kokkos_base.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAveForceKokkos<DeviceType>::FixAveForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixAveForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAveForceKokkos<DeviceType>::~FixAveForceKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAveForceKokkos<DeviceType>::init()
{
  FixAveForce::init();

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Cannot (yet) use respa with fix aveforce/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAveForceKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // update region if necessary

  if (region) {
    if (!(utils::strmatch(region->style, "^block") || utils::strmatch(region->style, "^sphere")))
      error->all(FLERR, "Cannot (yet) use {}-style region with fix aveforce/kk", region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("aveforce:k_match", atom->nlocal);
    KokkosBase *regionKKBase = dynamic_cast<KokkosBase *>(region);
    regionKKBase->match_all_kokkos(groupbit, k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // evaluate variables on host before kernel

  if (varflag == EQUAL) {
    atomKK->sync(Host, ALL_MASK);
    modify->clearstep_compute();
    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  atomKK->sync(execution_space, F_MASK | MASK_MASK);

  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  // phase 1: parallel reduce to sum forces and count

  double result[4] = {0.0, 0.0, 0.0, 0.0};

  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAveForceReduce>(0, nlocal),
                          *this, result);
  copymode = 0;

  // MPI allreduce on host

  MPI_Allreduce(result, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);

  int ncount = static_cast<int>(foriginal_all[3]);
  if (ncount == 0) return;

  // compute average force

  m_fave[0] = foriginal_all[0] / ncount + xvalue;
  m_fave[1] = foriginal_all[1] / ncount + yvalue;
  m_fave[2] = foriginal_all[2] / ncount + zvalue;

  // phase 2: apply average force to all participating atoms

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveForceApply>(0, nlocal), *this);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixAveForceKokkos<DeviceType>::operator()(TagFixAveForceReduce, const int &i,
                                               value_type result) const
{
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    result[0] += f(i,0);
    result[1] += f(i,1);
    result[2] += f(i,2);
    result[3] += 1.0;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixAveForceKokkos<DeviceType>::operator()(TagFixAveForceApply, const int &i) const
{
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    if (xstyle) f(i,0) = m_fave[0];
    if (ystyle) f(i,1) = m_fave[1];
    if (zstyle) f(i,2) = m_fave[2];
  }
}

namespace LAMMPS_NS {
template class FixAveForceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixAveForceKokkos<LMPHostType>;
#endif
}

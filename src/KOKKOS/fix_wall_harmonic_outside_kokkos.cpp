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

#include "fix_wall_harmonic_outside_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallHarmonicOutsideKokkos<DeviceType>::FixWallHarmonicOutsideKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallHarmonicOutside(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | MASK_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallHarmonicOutsideKokkos<DeviceType>::~FixWallHarmonicOutsideKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vatom, vatom);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallHarmonicOutsideKokkos<DeviceType>::v_setup_peratom(int vflag)
{
  // the per-atom virial is accumulated into a dual view, so the plain
  // base-class vatom array must not be allocated here (alloc = 0)

  v_init(vflag,0);

  // reallocate the per-atom virial dual view if necessary.  This has to happen
  // here, after v_init() has set vflag_atom and maxvatom for this step, and
  // not in post_force() before the base class runs, where both are stale.

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "wall_harmonic_outside:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallHarmonicOutsideKokkos<DeviceType>::post_force(int vflag)
{
  FixWallHarmonicOutside::post_force(vflag);

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   recalling force applied if outside the control volume
   and within the interaction cutoff
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
------------------------------------------------------------------------- */

template <class DeviceType>
void FixWallHarmonicOutsideKokkos<DeviceType>::wall_particle(int m_in, int which, double coord_in)
{
  m = m_in;
  coord = static_cast<KK_FLOAT>(coord_in);

  atomKK->sync(execution_space, datamask_read);
  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  int nlocal = atomKK->nlocal;

  dim = which / 2;
  side = which % 2;
  if (side == 0) side = -1;

  double result[13] = {0.0};

  copymode = 1;
  Kokkos::parallel_reduce(nlocal, *this, result);
  copymode = 0;

  ewall[0] += result[0];
  ewall[m+1] += result[m+1];
  atomKK->modified(execution_space, datamask_modify);

  if (vflag_global) {
    virial[0] += result[7];
    virial[1] += result[8];
    virial[2] += result[9];
    virial[3] += result[10];
    virial[4] += result[11];
    virial[5] += result[12];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallHarmonicOutsideKokkos<DeviceType>::operator()(const int &i, value_type result) const
{
  if (d_mask(i) & groupbit) {
    KK_FLOAT dr;
    if (side < 0) dr = coord - d_x(i,dim);
    else dr = d_x(i,dim) - coord;

    // no force if above the interaction cutoff or inside the control volume

    if (dr >= (KK_FLOAT) cutoff[m]) return;
    if (dr <= static_cast<KK_FLOAT>(0.0)) return;

    KK_FLOAT fwall = (KK_FLOAT) side * static_cast<KK_FLOAT>(2.0) * (KK_FLOAT) epsilon[m] * dr;
    d_f(i,dim) -= static_cast<KK_ACC_FLOAT>(fwall);
    result[0] += static_cast<double>((KK_FLOAT) epsilon[m] * dr * dr);
    result[m+1] += static_cast<double>(fwall);

    if (evflag) {
      KK_FLOAT vn;
      if (side < 0)
        vn = -fwall * dr;
      else
        vn = fwall * dr;
      v_tally(result, dim, i, vn);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallHarmonicOutsideKokkos<DeviceType>::v_tally(value_type result, int n, int i,
                                                       KK_FLOAT vn) const
{
  if (vflag_global)
    result[n+7] += static_cast<double>(vn);

  if (vflag_atom)
    Kokkos::atomic_add(&(d_vatom(i,n)), static_cast<KK_ACC_FLOAT>(vn));
}

namespace LAMMPS_NS {
template class FixWallHarmonicOutsideKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallHarmonicOutsideKokkos<LMPHostType>;
#endif
}

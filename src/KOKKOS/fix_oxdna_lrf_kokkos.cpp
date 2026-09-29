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

#include "fix_oxdna_lrf_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "memory_kokkos.h"
#include "update.h"
#include "verlet_kokkos.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaLRFKokkos<DeviceType>::FixOxdnaLRFKokkos(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  kokkosable = 1;
  avecEllipKK = nullptr;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  // Since this fix is called pre-force, datamsk_read can contain all read parameters
  // needed for the oxdna styles. This means each oxdna
  // style only needs to resync f, torque, energy and virial (which can change between
  // pair/bond styles).
  datamask_read = MASK_MASK | ELLIPSOID_MASK | BONUS_MASK |
                  X_MASK | TYPE_MASK | TAG_MASK | CG_DNA_MASK;
  datamask_modify = EMPTY_MASK;

  MemKK::realloc_kokkos(k_nx, "FixOxdnaLRF:nx", atom->nmax);
  MemKK::realloc_kokkos(k_ny, "FixOxdnaLRF:ny", atom->nmax);
  MemKK::realloc_kokkos(k_nz, "FixOxdnaLRF:nz", atom->nmax);
  d_nx = k_nx.template view<DeviceType>();
  d_ny = k_ny.template view<DeviceType>();
  d_nz = k_nz.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaLRFKokkos<DeviceType>::~FixOxdnaLRFKokkos() = default;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaLRFKokkos<DeviceType>::init()
{
  avecEllipKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
  if (!avecEllipKK) error->all(FLERR, "Fix OXDNA/LRF/kk requires atom style ellipsoid/kk");

  // with rRESPA the local reference frames would never be computed

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "The oxDNA styles do not support run style respa");

  // offer to zero the device forces and torques in the frame kernel, which
  // runs over all owned and ghost atoms anyway; VerletKokkos::setup() decides

  if (std::is_same_v<DeviceType, LMPDeviceType>) {
    auto verletKK = dynamic_cast<VerletKokkos *>(update->integrate);
    if (verletKK) verletKK->request_force_clear_by_fix(this);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixOxdnaLRFKokkos<DeviceType>::setmask()
{
  int mask = 0;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaLRFKokkos<DeviceType>::min_setup_pre_force(int vflag)
{
  min_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaLRFKokkos<DeviceType>::min_pre_force(int /*vflag*/)
{
  compute_lrf_kokkos(0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaLRFKokkos<DeviceType>::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaLRFKokkos<DeviceType>::pre_force(int /*vflag*/)
{
  auto verletKK = dynamic_cast<VerletKokkos *>(update->integrate);
  compute_lrf_kokkos((verletKK && (update->whichflag == 1)) ? verletKK->force_clear_by_fix(this) : 0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaLRFKokkos<DeviceType>::compute_lrf_kokkos(int zero_forces_flag)
{
  zero_forces = zero_forces_flag;
  if (zero_forces) {
    f = atomKK->k_f.template view<DeviceType>();
    torque = atomKK->k_torque.template view<DeviceType>();
  }

  if (atom->nmax > static_cast<int>(k_nx.extent(0))) {
    MemKK::realloc_kokkos(k_nx, "FixOxdnaLRFKokkos:nx", atom->nmax);
    MemKK::realloc_kokkos(k_ny, "FixOxdnaLRFKokkos:ny", atom->nmax);
    MemKK::realloc_kokkos(k_nz, "FixOxdnaLRFKokkos:nz", atom->nmax);
    d_nx = k_nx.template view<DeviceType>();
    d_ny = k_ny.template view<DeviceType>();
    d_nz = k_nz.template view<DeviceType>();
  }
  if (atom->nmax > static_cast<int>(d_xn.extent(0)))
    d_xn = decltype(d_xn)(Kokkos::view_alloc(Kokkos::WithoutInitializing, "FixOxdnaLRFKokkos:xn"),
                          atom->nmax);

  atomKK->sync(execution_space, datamask_read);

  mask = atomKK->k_mask.template view<DeviceType>();
  x = atomKK->k_x.template view<DeviceType>();
  ellipsoid = atomKK->k_ellipsoid.template view<DeviceType>();
  bonus = avecEllipKK->k_bonus.template view<DeviceType>();

  copymode = 1;
  // Frames are needed for all owned + ghost atoms (the max index any neighbor
  // list or bond list can reference); the slots in [nlocal+nghost, nmax) are
  // never read, so iterate nall rather than the full allocated nmax.
  const int nall = atom->nlocal + atom->nghost;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixOxdnaLRFComputeQuatToXYZ>(0, nall), *this);
  copymode = 0;

  k_nx.template modify<DeviceType>();
  k_ny.template modify<DeviceType>();
  k_nz.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixOxdnaLRFKokkos<DeviceType>::operator()(TagFixOxdnaLRFComputeQuatToXYZ, const int &i) const
{
  // frame vectors nx, ny, nz (zero for atoms outside the group or without ellipsoid)
  KK_FLOAT n[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  const int ne = ellipsoid(i);
  if ((mask(i) & groupbit) && (ne >= 0)) {
    const KK_FLOAT q0 = static_cast<KK_FLOAT>(bonus(ne).quat[0]);
    const KK_FLOAT q1 = static_cast<KK_FLOAT>(bonus(ne).quat[1]);
    const KK_FLOAT q2 = static_cast<KK_FLOAT>(bonus(ne).quat[2]);
    const KK_FLOAT q3 = static_cast<KK_FLOAT>(bonus(ne).quat[3]);

    const KK_FLOAT two = 2.0;

    n[0] = Kokkos::fma(q0, q0, Kokkos::fma(q1, q1, -Kokkos::fma(q2, q2, q3 * q3)));
    n[1] = two * Kokkos::fma(q1, q2, q0 * q3);
    n[2] = two * Kokkos::fma(q1, q3, -q0 * q2);

    n[3] = two * Kokkos::fma(q1, q2, -q0 * q3);
    n[4] = Kokkos::fma(q0, q0, Kokkos::fma(q2, q2, -Kokkos::fma(q1, q1, q3 * q3)));
    n[5] = two * Kokkos::fma(q2, q3, q0 * q1);

    n[6] = two * Kokkos::fma(q1, q3, q0 * q2);
    n[7] = two * Kokkos::fma(q2, q3, -q0 * q1);
    n[8] = Kokkos::fma(q0, q0, q3 * q3 - Kokkos::fma(q1, q1, q2 * q2));
  }

  d_nx(i, 0) = n[0];
  d_nx(i, 1) = n[1];
  d_nx(i, 2) = n[2];
  d_ny(i, 0) = n[3];
  d_ny(i, 1) = n[4];
  d_ny(i, 2) = n[5];
  d_nz(i, 0) = n[6];
  d_nz(i, 1) = n[7];
  d_nz(i, 2) = n[8];

  // packed record: position and frame vectors in one row
  d_xn(i, 0) = x(i, 0);
  d_xn(i, 1) = x(i, 1);
  d_xn(i, 2) = x(i, 2);
  d_xn(i, 3) = 0.0;
  for (int k = 0; k < 9; k++) d_xn(i, 4 + k) = n[k];
  d_xn(i, 13) = 0.0;
  d_xn(i, 14) = 0.0;
  d_xn(i, 15) = 0.0;

  // in place of VerletKokkos::force_clear()
  if (zero_forces) {
    f(i, 0) = 0.0;
    f(i, 1) = 0.0;
    f(i, 2) = 0.0;
    torque(i, 0) = 0.0;
    torque(i, 1) = 0.0;
    torque(i, 2) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixOxdnaLRFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixOxdnaLRFKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS

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

#include "fix_drag_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDragKokkos<DeviceType>::FixDragKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDrag(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDragKokkos<DeviceType>::~FixDragKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDragKokkos<DeviceType>::init()
{
  FixDrag::init();

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Cannot (yet) use respa with fix drag/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDragKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
  force_flag = 0;

  // capture domain data for device kernels

  m_triclinic = domain->triclinic;
  m_xperiodic = domain->xperiodic;
  m_yperiodic = domain->yperiodic;
  m_zperiodic = domain->zperiodic;
  m_xprd = (KK_FLOAT) domain->xprd;
  m_yprd = (KK_FLOAT) domain->yprd;
  m_zprd = (KK_FLOAT) domain->zprd;
  m_xprd_half = (KK_FLOAT) domain->xprd_half;
  m_yprd_half = (KK_FLOAT) domain->yprd_half;
  m_zprd_half = (KK_FLOAT) domain->zprd_half;
  m_xy = (KK_FLOAT) domain->xy;
  m_xz = (KK_FLOAT) domain->xz;
  m_yz = (KK_FLOAT) domain->yz;

  double result[3] = {0.0, 0.0, 0.0};

  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixDrag>(0, nlocal), *this, result);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);

  ftotal[0] = result[0];
  ftotal[1] = result[1];
  ftotal[2] = result[2];
  force_flag = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixDragKokkos<DeviceType>::operator()(TagFixDrag, const int &i, value_type result) const
{
  if (mask[i] & groupbit) {
    KK_FLOAT dx = x(i,0) - (KK_FLOAT) xc;
    KK_FLOAT dy = x(i,1) - (KK_FLOAT) yc;
    KK_FLOAT dz = x(i,2) - (KK_FLOAT) zc;
    if (!xflag) dx = 0.0;
    if (!yflag) dy = 0.0;
    if (!zflag) dz = 0.0;
    minimum_image(dx, dy, dz);
    KK_FLOAT r = sqrt(dx*dx + dy*dy + dz*dz);
    if (r > (KK_FLOAT) delta) {
      KK_FLOAT prefactor = (KK_FLOAT) f_mag / r;
      KK_FLOAT fx = prefactor * dx;
      KK_FLOAT fy = prefactor * dy;
      KK_FLOAT fz = prefactor * dz;
      f(i,0) -= fx;
      f(i,1) -= fy;
      f(i,2) -= fz;
      result[0] -= fx;
      result[1] -= fy;
      result[2] -= fz;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixDragKokkos<DeviceType>::minimum_image(KK_FLOAT &dx, KK_FLOAT &dy, KK_FLOAT &dz) const
{
  if (!m_triclinic) {
    if (m_xperiodic) {
      if (dx > m_xprd_half) dx -= m_xprd;
      else if (dx < -m_xprd_half) dx += m_xprd;
    }
    if (m_yperiodic) {
      if (dy > m_yprd_half) dy -= m_yprd;
      else if (dy < -m_yprd_half) dy += m_yprd;
    }
    if (m_zperiodic) {
      if (dz > m_zprd_half) dz -= m_zprd;
      else if (dz < -m_zprd_half) dz += m_zprd;
    }
  } else {
    if (m_zperiodic) {
      if (dz > m_zprd_half) {
        dz -= m_zprd;
        dy -= m_yz;
        dx -= m_xz;
      } else if (dz < -m_zprd_half) {
        dz += m_zprd;
        dy += m_yz;
        dx += m_xz;
      }
    }
    if (m_yperiodic) {
      if (dy > m_yprd_half) {
        dy -= m_yprd;
        dx -= m_xy;
      } else if (dy < -m_yprd_half) {
        dy += m_yprd;
        dx += m_xy;
      }
    }
    if (m_xperiodic) {
      if (dx > m_xprd_half) dx -= m_xprd;
      else if (dx < -m_xprd_half) dx += m_xprd;
    }
  }
}

namespace LAMMPS_NS {
template class FixDragKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixDragKokkos<LMPHostType>;
#endif
}

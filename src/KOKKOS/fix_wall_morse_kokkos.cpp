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

#include "fix_wall_morse_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "memory_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallMorseKokkos<DeviceType>::FixWallMorseKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallMorse(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK;
  datamask_modify = F_MASK;

  memoryKK->create_kokkos(k_cutoff,  6, "wall_morse:cutoff");
  memoryKK->create_kokkos(k_coeff1,  6, "wall_morse:coeff1");
  memoryKK->create_kokkos(k_offset,  6, "wall_morse:offset");
  memoryKK->create_kokkos(k_alpha,   6, "wall_morse:alpha");
  memoryKK->create_kokkos(k_sigma,   6, "wall_morse:sigma");
  memoryKK->create_kokkos(k_epsilon, 6, "wall_morse:epsilon");

  d_cutoff  = k_cutoff.template view<DeviceType>();
  d_coeff1  = k_coeff1.template view<DeviceType>();
  d_offset  = k_offset.template view<DeviceType>();
  d_alpha   = k_alpha.template view<DeviceType>();
  d_sigma   = k_sigma.template view<DeviceType>();
  d_epsilon = k_epsilon.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallMorseKokkos<DeviceType>::~FixWallMorseKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_cutoff);
  memoryKK->destroy_kokkos(k_coeff1);
  memoryKK->destroy_kokkos(k_offset);
  memoryKK->destroy_kokkos(k_alpha);
  memoryKK->destroy_kokkos(k_sigma);
  memoryKK->destroy_kokkos(k_epsilon);
  memoryKK->destroy_kokkos(k_vatom, vatom);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallMorseKokkos<DeviceType>::precompute(int m_in)
{
  FixWallMorse::precompute(m_in);

  for (int i = 0; i < 6; i++) {
    k_cutoff.view_host()(i)  = cutoff[i];
    k_coeff1.view_host()(i)  = coeff1[i];
    k_offset.view_host()(i)  = offset[i];
    k_alpha.view_host()(i)   = alpha[i];
    k_sigma.view_host()(i)   = sigma[i];
    k_epsilon.view_host()(i) = epsilon[i];
  }

  k_cutoff.modify_host();
  k_coeff1.modify_host();
  k_offset.modify_host();
  k_alpha.modify_host();
  k_sigma.modify_host();
  k_epsilon.modify_host();

  k_cutoff.template sync<DeviceType>();
  k_coeff1.template sync<DeviceType>();
  k_offset.template sync<DeviceType>();
  k_alpha.template sync<DeviceType>();
  k_sigma.template sync<DeviceType>();
  k_epsilon.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallMorseKokkos<DeviceType>::post_force(int vflag)
{
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "wall_morse:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  FixWallMorse::post_force(vflag);

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

template <class DeviceType>
void FixWallMorseKokkos<DeviceType>::wall_particle(int m_in, int which, double coord_in)
{
  m = m_in;
  coord = coord_in;

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
void FixWallMorseKokkos<DeviceType>::operator()(const int &i, value_type result) const
{
  if (d_mask(i) & groupbit) {
    KK_FLOAT delta;
    if (side < 0) delta = d_x(i,dim) - coord;
    else delta = coord - d_x(i,dim);
    if (delta >= d_cutoff(m)) return;
    if (delta <= 0.0)
      Kokkos::abort("Particle on or inside fix wall surface");
    KK_FLOAT dr = delta - d_sigma(m);
    KK_FLOAT dexp = Kokkos::exp(-d_alpha(m) * dr);
    KK_FLOAT fwall = (KK_FLOAT) side * d_coeff1(m) * (dexp*dexp - dexp);
    d_f(i,dim) -= fwall;
    result[0] += d_epsilon(m) * (dexp*dexp - 2.0*dexp) - d_offset(m);
    result[m+1] += fwall;

    if (evflag) {
      KK_FLOAT vn;
      if (side < 0)
        vn = -fwall * delta;
      else
        vn = fwall * delta;
      v_tally(result, dim, i, vn);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallMorseKokkos<DeviceType>::v_tally(value_type result, int n, int i,
                                             KK_FLOAT vn) const
{
  if (vflag_global)
    result[n+7] += vn;

  if (vflag_atom)
    Kokkos::atomic_add(&(d_vatom(i,n)), vn);
}

namespace LAMMPS_NS {
template class FixWallMorseKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallMorseKokkos<LMPHostType>;
#endif
}

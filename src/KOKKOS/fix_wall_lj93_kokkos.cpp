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
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "fix_wall_lj93_kokkos.h"

#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"
#include "input.h"
#include "memory_kokkos.h"

#include <iostream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallLJ93Kokkos<DeviceType>::FixWallLJ93Kokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallLJ93(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | MASK_MASK;
  datamask_modify = F_MASK;

  memoryKK->create_kokkos(k_cutoff,6,"wall_lj93:cutoff");
  memoryKK->create_kokkos(k_coeff1,6,"wall_lj93:coeff1");
  memoryKK->create_kokkos(k_coeff2,6,"wall_lj93:coeff2");
  memoryKK->create_kokkos(k_coeff3,6,"wall_lj93:coeff3");
  memoryKK->create_kokkos(k_coeff4,6,"wall_lj93:coeff4");
  memoryKK->create_kokkos(k_offset,6,"wall_lj93:offset");

  d_cutoff = k_cutoff.template view<DeviceType>();
  d_coeff1 = k_coeff1.template view<DeviceType>();
  d_coeff2 = k_coeff2.template view<DeviceType>();
  d_coeff3 = k_coeff3.template view<DeviceType>();
  d_coeff4 = k_coeff4.template view<DeviceType>();
  d_offset = k_offset.template view<DeviceType>();
}

template<class DeviceType>
FixWallLJ93Kokkos<DeviceType>::~FixWallLJ93Kokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_cutoff);
  memoryKK->destroy_kokkos(k_coeff1);
  memoryKK->destroy_kokkos(k_coeff2);
  memoryKK->destroy_kokkos(k_coeff3);
  memoryKK->destroy_kokkos(k_coeff4);
  memoryKK->destroy_kokkos(k_offset);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallLJ93Kokkos<DeviceType>::precompute(int m)
{
  FixWallLJ93::precompute(m);

  for( int i=0 ; i<6 ; i++ ) {
    k_cutoff.h_view(i) = cutoff[i];
    k_coeff1.h_view(i) = coeff1[i];
    k_coeff2.h_view(i) = coeff2[i];
    k_coeff3.h_view(i) = coeff3[i];
    k_coeff4.h_view(i) = coeff4[i];
    k_offset.h_view(i) = offset[i];
  }

  k_cutoff.template modify<LMPHostType>();
  k_coeff1.template modify<LMPHostType>();
  k_coeff2.template modify<LMPHostType>();
  k_coeff3.template modify<LMPHostType>();
  k_coeff4.template modify<LMPHostType>();
  k_offset.template modify<LMPHostType>();

  k_cutoff.template sync<DeviceType>();
  k_coeff1.template sync<DeviceType>();
  k_coeff2.template sync<DeviceType>();
  k_coeff3.template sync<DeviceType>();
  k_coeff4.template sync<DeviceType>();
  k_offset.template sync<DeviceType>();

}


/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallLJ93Kokkos<DeviceType>::post_force(int vflag)
{

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"wall_lj93:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  FixWallLJ93::post_force(vflag);

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

template <class DeviceType>
void FixWallLJ93Kokkos<DeviceType>::wall_particle(int m_in, int which, double coord_in)
{
  m = m_in;
  coord = coord_in;

  atomKK->sync(execution_space,datamask_read);
  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  int nlocal = atomKK->nlocal;

  dim = which / 2;
  side = which % 2;
  if (side == 0) side = -1;

  double result[13] = {0.0};

  copymode = 1;
  Kokkos::parallel_reduce(nlocal,*this,result);
  copymode = 0;

  ewall[0] += result[0];
  ewall[m+1] += result[m+1];
  atomKK->modified(execution_space,datamask_modify);

  if (vflag_global) {
    virial[0] += result[7];
    virial[1] += result[8];
    virial[2] += result[9];
    virial[3] += result[10];
    virial[4] += result[11];
    virial[5] += result[12];
  }
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixWallLJ93Kokkos<DeviceType>::operator()(const int &i, value_type result) const {
  if (d_mask(i) & groupbit) {
    double delta;
    if (side < 0) delta = d_x(i,dim) - coord;
    else delta = coord - d_x(i,dim);
    if (delta >= d_cutoff(m)) return;
    if (delta <= 0.0)
      Kokkos::abort("Particle on or inside fix wall surface");
    double rinv = 1.0/delta;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;
    double r10inv = r4inv*r4inv*r2inv;
    double fwall = side * (d_coeff1(m)*r10inv - d_coeff2(m)*r4inv);
    d_f(i,dim) -= fwall;
    result[0] += d_coeff3(m)*r4inv*r4inv*rinv - d_coeff4(m)*r2inv*rinv - d_offset(m);
    result[m+1] += fwall;

    if (evflag) {
      double vn;
      if (side < 0)
        vn = -fwall * delta;
      else
        vn = fwall * delta;
      v_tally(result, dim, i, vn);
    }
  }
}

/* ----------------------------------------------------------------------
   tally virial component into global and per-atom accumulators
   n = index of virial component (0-5)
   i = local index of atom
   vn = nth component of virial for the interaction
   increment nth component of global virial by vn
   increment nth component of per-atom virial by vn
   this method can be used when fix computes forces in post_force()
   and the force depends on a distance to some external object
     e.g. fix wall/lj93: compute virial only on owned atoms
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixWallLJ93Kokkos<DeviceType>::v_tally(value_type result, int n, int i, double vn) const
{
  if (vflag_global)
    result[n+7] += vn;

  if (vflag_atom)
    Kokkos::atomic_add(&(d_vatom(i,n)),vn);
}

namespace LAMMPS_NS {
template class FixWallLJ93Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallLJ93Kokkos<LMPHostType>;
#endif
}

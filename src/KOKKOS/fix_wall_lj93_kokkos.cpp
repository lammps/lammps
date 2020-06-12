/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_wall_lj93_kokkos.h"
#include <cmath>
#include "atom_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
FixWallLJ93Kokkos<Space>::FixWallLJ93Kokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallLJ93(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  virial_flag = 0;
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixWallLJ93Kokkos<Space>::wall_particle(int m_in, int which, double coord_in)
{
  m = m_in;
  coord = coord_in;

  atomKK->sync(Space, X_MASK|F_MASK|MASK_MASK);
  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  DAT::tdual_int_scalar k_oneflag = DAT::tdual_int_scalar("fix:oneflag");
  d_oneflag = DualViewHelper<Space>::view(k_oneflag);

  int nlocal = atom->nlocal;

  dim = which / 2;
  side = which % 2;
  if (side == 0) side = -1;

  copymode = 1;
  FixWallLJ93KokkosFunctor<Space> wp_functor(this);
  Kokkos::parallel_reduce(nlocal,wp_functor,ewall);
  copymode = 0;

  atomKK->modified(Space, F_MASK);

  DualViewHelper<Space>::modify(k_oneflag);
  k_oneflag.sync_host();
  if (k_oneflag.h_view()) error->one(FLERR,"Particle on or inside fix wall surface");
}

template <ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixWallLJ93Kokkos<Space>::wall_particle_item(int i, value_type ewall) const {
  if (mask(i) & groupbit) {
    KK_FLOAT delta;
    if (side < 0) delta = x(i,dim) - coord;
    else delta = coord - x(i,dim);
    if (delta >= cutoff[m]) return;
    if (delta <= 0.0) {
      d_oneflag() = 1;
      return;
    }
    KK_FLOAT rinv = 1.0/delta;
    KK_FLOAT r2inv = rinv*rinv;
    KK_FLOAT r4inv = r2inv*r2inv;
    KK_FLOAT r10inv = r4inv*r4inv*r2inv;
    KK_FLOAT fwall = side * (coeff1[m]*r10inv - coeff2[m]*r4inv);
    f(i,dim) -= fwall;
    ewall[0] += coeff3[m]*r4inv*r4inv*rinv -
      coeff4[m]*r2inv*rinv - offset[m];
    ewall[m+1] += fwall;
  }
}

namespace LAMMPS_NS {
template class FixWallLJ93Kokkos<Device>;
template class FixWallLJ93Kokkos<Host>;
}

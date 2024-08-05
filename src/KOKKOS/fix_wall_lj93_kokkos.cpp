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
#include "math_special.h"
#include "memory_kokkos.h"
#include "modify_kokkos.h"
#include "variable.h"
#include "update.h"

using namespace LAMMPS_NS;
using MathSpecial::powint;

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

  memoryKK->create_kokkos(k_epsilon,6,"wall_lj93:epsilon");
  memoryKK->create_kokkos(k_sigma,6,"wall_lj93:sigma");
  memoryKK->create_kokkos(k_cutoff,6,"wall_lj93:cutoff");

  d_epsilon = k_epsilon.template view<DeviceType>();
  d_sigma = k_sigma.template view<DeviceType>();
  d_cutoff = k_cutoff.template view<DeviceType>();

  for( int i=0 ; i<6 ; i++ ) {
    k_epsilon.h_view(i) = epsilon[i];
    k_sigma.h_view(i) = sigma[i];
    k_cutoff.h_view(i) = cutoff[i];
  }

  k_epsilon.template modify<LMPHostType>();
  k_sigma.template modify<LMPHostType>();
  k_cutoff.template modify<LMPHostType>();

  k_epsilon.template sync<DeviceType>();
  k_sigma.template sync<DeviceType>();
  k_cutoff.template sync<DeviceType>();

  memoryKK->create_kokkos(d_coeff1,6,"wall_lj93:coeff1");
  memoryKK->create_kokkos(d_coeff2,6,"wall_lj93:coeff2");
  memoryKK->create_kokkos(d_coeff3,6,"wall_lj93:coeff3");
  memoryKK->create_kokkos(d_coeff4,6,"wall_lj93:coeff4");
  memoryKK->create_kokkos(d_offset,6,"wall_lj93:offset");

  memoryKK->create_kokkos(k_ewall,ewall,7,"wall_lj93:ewall");
  d_ewall = k_ewall.template view<DeviceType>();

}

template<class DeviceType>
FixWallLJ93Kokkos<DeviceType>::~FixWallLJ93Kokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_epsilon);
  memoryKK->destroy_kokkos(k_sigma);
  memoryKK->destroy_kokkos(k_cutoff);

  memoryKK->destroy_kokkos(d_coeff1);
  memoryKK->destroy_kokkos(d_coeff2);
  memoryKK->destroy_kokkos(d_coeff3);
  memoryKK->destroy_kokkos(d_coeff4);
  memoryKK->destroy_kokkos(d_offset);

  memoryKK->destroy_kokkos(k_vatom,vatom);
  memoryKK->destroy_kokkos(k_ewall,ewall);

}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallLJ93Kokkos<DeviceType>::precompute(int m)
{
  d_coeff1[m] = 6.0 / 5.0 * d_epsilon[m] * powint(d_sigma[m], 9);
  d_coeff2[m] = 3.0 * d_epsilon[m] * powint(d_sigma[m], 3);
  d_coeff3[m] = 2.0 / 15.0 * d_epsilon[m] * powint(d_sigma[m], 9);
  d_coeff4[m] = d_epsilon[m] * powint(d_sigma[m], 3);

  double rinv = 1.0 / d_cutoff[m];
  double r2inv = rinv * rinv;
  double r4inv = r2inv * r2inv;
  d_offset[m] = d_coeff3[m] * r4inv * r4inv * rinv - d_coeff4[m] * r2inv * rinv;
}


/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallLJ93Kokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  // virial setup

  v_init(vflag);

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"wall_lj93:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  // energy intialize.
  // eflag is used to track whether wall energies have been communicated.

  eflag = 0;
  for (int m = 0; m <= nwall; m++) d_ewall(m) = 0.0;

  // coord = current position of wall
  // evaluate variables if necessary, wrap with clear/add
  // for epsilon/sigma variables need to re-invoke precompute()

  if (varflag) modify->clearstep_compute();

  double coord;

  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(xindex[m]);
      if (wallwhich[m] < YLO)
        coord *= xscale;
      else if (wallwhich[m] < ZLO)
        coord *= yscale;
      else
        coord *= zscale;
    } else
      coord = coord0[m];
    if (wstyle[m] == VARIABLE) {
      if (estyle[m] == VARIABLE) {
        d_epsilon[m] = input->variable->compute_equal(eindex[m]);
        if (d_epsilon[m] < 0.0) error->all(FLERR, "Variable evaluation in fix wall gave bad value");
      }
      if (sstyle[m] == VARIABLE) {
        d_sigma[m] = input->variable->compute_equal(sindex[m]);
        if (d_sigma[m] < 0.0) error->all(FLERR, "Variable evaluation in fix wall gave bad value");
      }
      precompute(m);
    }

    wall_particle(m, wallwhich[m], coord);

  }

  k_ewall.template modify<DeviceType>();
  k_ewall.template sync<LMPHostType>();

  if (varflag) modify->addstep_compute(update->ntimestep + 1);

  atomKK->modified(execution_space,F_MASK);

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
KOKKOS_INLINE_FUNCTION
void FixWallLJ93Kokkos<DeviceType>::wall_particle(int m_in, int which, double coord_in)
{
  m = m_in;
  coord = coord_in;
  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  int nlocal = atomKK->nlocal;

  dim = which / 2;
  side = which % 2;
  if (side == 0) side = -1;

  double result[13];

  copymode = 1;
  FixWallLJ93KokkosFunctor<DeviceType> functor(this);
  Kokkos::parallel_reduce(nlocal,functor,result);
  copymode = 0;

  Kokkos::atomic_add(&(d_ewall[0]),result[0]);
  Kokkos::atomic_add(&(d_ewall[m+1]),result[m+1]);

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
void FixWallLJ93Kokkos<DeviceType>::wall_particle_item(int i, value_type result) const {
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

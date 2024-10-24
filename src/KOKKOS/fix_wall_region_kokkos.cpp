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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_wall_region_kokkos.h"

#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"
#include "kokkos_base.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "region.h"
#include "region_block_kokkos.h"
#include "region_sphere_kokkos.h"

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;

enum { LJ93, LJ126, LJ1043, COLLOID, HARMONIC, MORSE };

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallRegionKokkos<DeviceType>::FixWallRegionKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallRegion(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | MASK_MASK;
  datamask_modify = F_MASK;
}

template<class DeviceType>
FixWallRegionKokkos<DeviceType>::~FixWallRegionKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallRegionKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  // virial setup

  v_init(vflag);

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"wall_region:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  if (style == COLLOID) d_radius = atomKK->k_radius.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  int nlocal = atomKK->nlocal;

  region->prematch();

  // virial setup

  v_init(vflag);

  // region->match() ensures particle is in region or on surface, else error
  // if returned contact dist r = 0, is on surface, also an error
  // in COLLOID case, r <= radius is an error
  // initilize ewall after region->prematch(),
  //   so a dynamic region can access last timestep values

  // energy intialize
  // eflag is used to track whether wall energies have been communicated

  eflag = 0;
  double result[10];
  copymode = 1;

  if(auto *regionKK = dynamic_cast<RegBlockKokkos<DeviceType>*>(region)) {
    FixWallRegionKokkosFunctor<DeviceType,class RegBlockKokkos<DeviceType>> functor(this,regionKK);
    Kokkos::parallel_reduce(nlocal,functor,result);
  } else if (auto *regionKK = dynamic_cast<RegSphereKokkos<DeviceType>*>(region)){
    FixWallRegionKokkosFunctor<DeviceType,class RegSphereKokkos<DeviceType>> functor(this,regionKK);
    Kokkos::parallel_reduce(nlocal,functor,result);
  }

  copymode = 0;
  for( int i=0 ; i<4 ; i++ ) ewall[i] = result[i];

  if (vflag_global) {
    virial[0] += result[4];
    virial[1] += result[5];
    virial[2] += result[6];
    virial[3] += result[7];
    virial[4] += result[8];
    virial[5] += result[9];
  }

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

template<class DeviceType>
template<class T>
KOKKOS_INLINE_FUNCTION
void FixWallRegionKokkos<DeviceType>::wall_particle(T regionKK, const int i, value_type result) const {
  if (d_mask(i) & groupbit) {

    if (!regionKK->match_kokkos(d_x(i,0), d_x(i,1), d_x(i,2))) Kokkos::abort("Particle outside surface of region used in fix wall/region");

    double rinv, tooclose;

    if (style == COLLOID)
      tooclose = d_radius(i);
    else
      tooclose = 0.0;

    int n = regionKK->surface_kokkos(d_x(i,0), d_x(i,1), d_x(i,2), cutoff);

    for ( int m = 0; m < n; m++) {

      double r = regionKK->d_contact[m].r;
      double delx = regionKK->d_contact[m].delx;
      double dely = regionKK->d_contact[m].dely;
      double delz = regionKK->d_contact[m].delz;

      if (r <= tooclose)
        Kokkos::abort("Particle outside surface of region used in fix wall/region");
      else
        rinv = 1.0 / r;

      double fwallKK, engKK;

      if (style == LJ93) engKK = lj93(r,fwallKK);
      else if (style == LJ126) engKK = lj126(r,fwallKK);
      else if (style == LJ1043) engKK = lj1043(r,fwallKK);
      else if (style == MORSE) engKK = morse(r,fwallKK);
      else if (style == COLLOID) engKK = colloid(r,d_radius(i),fwallKK);
      else engKK = harmonic(r,fwallKK);

      double fx = fwallKK * delx * rinv;
      double fy = fwallKK * dely * rinv;
      double fz = fwallKK * delz * rinv;
      d_f(i,0) += fx;
      d_f(i,1) += fy;
      d_f(i,2) += fz;
      result[1] -= fx;
      result[2] -= fy;
      result[3] -= fz;
      result[0] += engKK;
      if (evflag) {
        double v[6] = {
          fx * delx,
          fy * dely,
          fz * delz,
          fx * dely,
          fx * delz,
          fy * delz
        };
        v_tally(result,i,v);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   LJ 9/3 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixWallRegionKokkos<DeviceType>::lj93(double r, double& fwallKK) const
{
  double rinv = 1.0 / r;
  double r2inv = rinv * rinv;
  double r4inv = r2inv * r2inv;
  double r10inv = r4inv * r4inv * r2inv;
  fwallKK = coeff1 * r10inv - coeff2 * r4inv;
  return coeff3 * r4inv * r4inv * rinv - coeff4 * r2inv * rinv - offset;
}

/* ----------------------------------------------------------------------
   LJ 12/6 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixWallRegionKokkos<DeviceType>::lj126(double r, double& fwallKK) const
{
  double rinv = 1.0 / r;
  double r2inv = rinv * rinv;
  double r6inv = r2inv * r2inv * r2inv;
  fwallKK = r6inv * (coeff1 * r6inv - coeff2) * rinv;
  return r6inv * (coeff3 * r6inv - coeff4) - offset;
}

/* ----------------------------------------------------------------------
   LJ 10/4/3 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixWallRegionKokkos<DeviceType>::lj1043(double r, double& fwallKK) const
{
  double rinv = 1.0 / r;
  double r2inv = rinv * rinv;
  double r4inv = r2inv * r2inv;
  double r10inv = r4inv * r4inv * r2inv;
  fwallKK = coeff5 * r10inv * rinv - coeff6 * r4inv * rinv - coeff7 * powint(r + coeff4, -4);
  return coeff1 * r10inv - coeff2 * r4inv - coeff3 * powint(r + coeff4, -3) - offset;
}

/* ----------------------------------------------------------------------
   Morse interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixWallRegionKokkos<DeviceType>::morse(double r, double& fwallKK) const
{
  double dr = r - sigma;
  double dexp = exp(-alpha * dr);
  fwallKK = coeff1 * (dexp * dexp - dexp);
  return epsilon * (dexp * dexp - 2.0 * dexp) - offset;
}

/* ----------------------------------------------------------------------
   colloid interaction for finite-size particle of rad with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixWallRegionKokkos<DeviceType>::colloid(double r, double rad, double& fwallKK) const
{
  double new_coeff2 = coeff2 * rad * rad * rad;
  double diam = 2.0 * rad;

  double rad2 = rad * rad;
  double rad4 = rad2 * rad2;
  double rad8 = rad4 * rad4;
  double delta2 = rad2 - r * r;
  double rinv = 1.0 / delta2;
  double r2inv = rinv * rinv;
  double r4inv = r2inv * r2inv;
  double r8inv = r4inv * r4inv;
  fwallKK = coeff1 *
          (rad8 * rad + 27.0 * rad4 * rad2 * rad * r * r + 63.0 * rad4 * rad * powint(r, 4) +
           21.0 * rad2 * rad * powint(r, 6)) *
          r8inv -
      new_coeff2 * r2inv;

  double r2 = 0.5 * diam - r;
  double rinv2 = 1.0 / r2;
  double r2inv2 = rinv2 * rinv2;
  double r4inv2 = r2inv2 * r2inv2;
  double r3 = r + 0.5 * diam;
  double rinv3 = 1.0 / r3;
  double r2inv3 = rinv3 * rinv3;
  double r4inv3 = r2inv3 * r2inv3;
  return coeff3 *
          ((-3.5 * diam + r) * r4inv2 * r2inv2 * rinv2 +
           (3.5 * diam + r) * r4inv3 * r2inv3 * rinv3) -
      coeff4 * ((-diam * r + r2 * r3 * (log(-r2) - log(r3))) * (-rinv2) * rinv3) - offset;
}

/* ----------------------------------------------------------------------
   harmonic interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixWallRegionKokkos<DeviceType>::harmonic(double r, double& fwallKK) const
{
  double dr = cutoff - r;
  fwallKK = 2.0 * epsilon * dr;
  return epsilon * dr * dr;
}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   i = local index of atom
   v = total virial for the interaction
   increment global virial by v
   increment per-atom virial by v
   this method can be used when fix computes forces in post_force()
   and the force depends on a distance to some external object
     e.g. fix wall/lj93: compute virial only on owned atoms
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixWallRegionKokkos<DeviceType>::v_tally(value_type result, int i, double *v) const
{
  if (vflag_global) {
    result[4] += v[0];
    result[5] += v[1];
    result[6] += v[2];
    result[7] += v[3];
    result[8] += v[4];
    result[9] += v[5];
  }

  if (vflag_atom) {
    Kokkos::atomic_add(&(d_vatom(i,0)),v[0]);
    Kokkos::atomic_add(&(d_vatom(i,1)),v[1]);
    Kokkos::atomic_add(&(d_vatom(i,2)),v[2]);
    Kokkos::atomic_add(&(d_vatom(i,3)),v[3]);
    Kokkos::atomic_add(&(d_vatom(i,4)),v[4]);
    Kokkos::atomic_add(&(d_vatom(i,5)),v[5]);
  }
}

namespace LAMMPS_NS {
template class FixWallRegionKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallRegionKokkos<LMPHostType>;
#endif
}

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

#include <type_traits>

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
  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK;
  datamask_modify = F_MASK;
}

template<class DeviceType>
FixWallRegionKokkos<DeviceType>::~FixWallRegionKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallRegionKokkos<DeviceType>::init()
{
  FixWallRegion::init();

  // without this check a region w/o KOKKOS support was silently ignored:
  // no wall forces and undefined energy/virial contributions

  if (!dynamic_cast<RegBlockKokkos<DeviceType>*>(region) &&
      !dynamic_cast<RegSphereKokkos<DeviceType>*>(region))
    error->all(FLERR,"Fix wall/region/kk requires region style block/kk or sphere/kk");
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallRegionKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  // virial setup

  // the per-atom virial is accumulated into a dual view, so the plain
  // base-class vatom array must not be allocated here (alloc = 0)

  v_init(vflag,0);

  // reallocate the per-atom virial dual view if necessary

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
    k_vatom.sync_host();
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallRegionKokkos<DeviceType>::wall_particle(T regionKK, const int i, value_type result) const {
  if (d_mask(i) & groupbit) {

    if (!regionKK->match_kokkos(static_cast<double>(d_x(i,0)), static_cast<double>(d_x(i,1)), static_cast<double>(d_x(i,2)))) Kokkos::abort("Particle outside surface of region used in fix wall/region");

    KK_FLOAT rinv, tooclose;

    if (style == COLLOID)
      tooclose = d_radius(i);
    else
      tooclose = 0.0;

    // the contact list lives on the stack, so that concurrently running
    // threads do not overwrite each other's contacts

    Region::Contact contact[std::remove_pointer_t<T>::MAXCONTACT];

    int n = regionKK->surface_kokkos(static_cast<double>(d_x(i,0)),
                                     static_cast<double>(d_x(i,1)),
                                     static_cast<double>(d_x(i,2)), cutoff, contact);

    for ( int m = 0; m < n; m++) {

      KK_FLOAT r = static_cast<KK_FLOAT>(contact[m].r);
      KK_FLOAT delx = static_cast<KK_FLOAT>(contact[m].delx);
      KK_FLOAT dely = static_cast<KK_FLOAT>(contact[m].dely);
      KK_FLOAT delz = static_cast<KK_FLOAT>(contact[m].delz);

      if (r <= tooclose)
        Kokkos::abort("Particle outside surface of region used in fix wall/region");
      else
        rinv = static_cast<KK_FLOAT>(1.0) / r;

      KK_FLOAT fwallKK, engKK;

      if (style == LJ93) engKK = lj93(r,fwallKK);
      else if (style == LJ126) engKK = lj126(r,fwallKK);
      else if (style == LJ1043) engKK = lj1043(r,fwallKK);
      else if (style == MORSE) engKK = morse(r,fwallKK);
      else if (style == COLLOID) engKK = colloid(r,d_radius(i),fwallKK);
      else engKK = harmonic(r,fwallKK);

      KK_FLOAT fx = fwallKK * delx * rinv;
      KK_FLOAT fy = fwallKK * dely * rinv;
      KK_FLOAT fz = fwallKK * delz * rinv;
      d_f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
      d_f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
      d_f(i,2) += static_cast<KK_ACC_FLOAT>(fz);
      result[1] -= static_cast<double>(fx);
      result[2] -= static_cast<double>(fy);
      result[3] -= static_cast<double>(fz);
      result[0] += static_cast<double>(engKK);
      if (evflag) {
        KK_FLOAT v[6] = {
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT FixWallRegionKokkos<DeviceType>::lj93(KK_FLOAT r, KK_FLOAT& fwallKK) const
{
  const KK_FLOAT coeff1_kk = static_cast<KK_FLOAT>(coeff1);
  const KK_FLOAT coeff2_kk = static_cast<KK_FLOAT>(coeff2);
  const KK_FLOAT coeff3_kk = static_cast<KK_FLOAT>(coeff3);
  const KK_FLOAT coeff4_kk = static_cast<KK_FLOAT>(coeff4);
  const KK_FLOAT offset_kk = static_cast<KK_FLOAT>(offset);
  KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
  KK_FLOAT r2inv = rinv * rinv;
  KK_FLOAT r4inv = r2inv * r2inv;
  KK_FLOAT r10inv = r4inv * r4inv * r2inv;
  fwallKK = coeff1_kk * r10inv - coeff2_kk * r4inv;
  return coeff3_kk * r4inv * r4inv * rinv - coeff4_kk * r2inv * rinv - offset_kk;
}

/* ----------------------------------------------------------------------
   LJ 12/6 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT FixWallRegionKokkos<DeviceType>::lj126(KK_FLOAT r, KK_FLOAT& fwallKK) const
{
  const KK_FLOAT coeff1_kk = static_cast<KK_FLOAT>(coeff1);
  const KK_FLOAT coeff2_kk = static_cast<KK_FLOAT>(coeff2);
  const KK_FLOAT coeff3_kk = static_cast<KK_FLOAT>(coeff3);
  const KK_FLOAT coeff4_kk = static_cast<KK_FLOAT>(coeff4);
  const KK_FLOAT offset_kk = static_cast<KK_FLOAT>(offset);
  KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
  KK_FLOAT r2inv = rinv * rinv;
  KK_FLOAT r6inv = r2inv * r2inv * r2inv;
  fwallKK = r6inv * (coeff1_kk * r6inv - coeff2_kk) * rinv;
  return r6inv * (coeff3_kk * r6inv - coeff4_kk) - offset_kk;
}

/* ----------------------------------------------------------------------
   LJ 10/4/3 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT FixWallRegionKokkos<DeviceType>::lj1043(KK_FLOAT r, KK_FLOAT& fwallKK) const
{
  const KK_FLOAT coeff1_kk = static_cast<KK_FLOAT>(coeff1);
  const KK_FLOAT coeff2_kk = static_cast<KK_FLOAT>(coeff2);
  const KK_FLOAT coeff3_kk = static_cast<KK_FLOAT>(coeff3);
  const KK_FLOAT coeff4_kk = static_cast<KK_FLOAT>(coeff4);
  const KK_FLOAT coeff5_kk = static_cast<KK_FLOAT>(coeff5);
  const KK_FLOAT coeff6_kk = static_cast<KK_FLOAT>(coeff6);
  const KK_FLOAT coeff7_kk = static_cast<KK_FLOAT>(coeff7);
  const KK_FLOAT offset_kk = static_cast<KK_FLOAT>(offset);
  KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
  KK_FLOAT r2inv = rinv * rinv;
  KK_FLOAT r4inv = r2inv * r2inv;
  KK_FLOAT r10inv = r4inv * r4inv * r2inv;
  fwallKK = coeff5_kk * r10inv * rinv - coeff6_kk * r4inv * rinv - coeff7_kk * powint(r + coeff4_kk, -4);
  return coeff1_kk * r10inv - coeff2_kk * r4inv - coeff3_kk * powint(r + coeff4_kk, -3) - offset_kk;
}

/* ----------------------------------------------------------------------
   Morse interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT FixWallRegionKokkos<DeviceType>::morse(KK_FLOAT r, KK_FLOAT& fwallKK) const
{
  const KK_FLOAT sigma_kk = static_cast<KK_FLOAT>(sigma);
  const KK_FLOAT alpha_kk = static_cast<KK_FLOAT>(alpha);
  const KK_FLOAT coeff1_kk = static_cast<KK_FLOAT>(coeff1);
  const KK_FLOAT epsilon_kk = static_cast<KK_FLOAT>(epsilon);
  const KK_FLOAT offset_kk = static_cast<KK_FLOAT>(offset);
  KK_FLOAT dr = r - sigma_kk;
  KK_FLOAT dexp = Kokkos::exp(-alpha_kk * dr);
  fwallKK = coeff1_kk * (dexp * dexp - dexp);
  return epsilon_kk * (dexp * dexp - static_cast<KK_FLOAT>(2.0) * dexp) - offset_kk;
}

/* ----------------------------------------------------------------------
   colloid interaction for finite-size particle of rad with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT FixWallRegionKokkos<DeviceType>::colloid(KK_FLOAT r, KK_FLOAT rad, KK_FLOAT& fwallKK) const
{
  const KK_FLOAT coeff1_kk = static_cast<KK_FLOAT>(coeff1);
  const KK_FLOAT coeff2_kk = static_cast<KK_FLOAT>(coeff2);
  const KK_FLOAT coeff3_kk = static_cast<KK_FLOAT>(coeff3);
  const KK_FLOAT coeff4_kk = static_cast<KK_FLOAT>(coeff4);
  const KK_FLOAT offset_kk = static_cast<KK_FLOAT>(offset);
  KK_FLOAT new_coeff2 = coeff2_kk * rad * rad * rad;
  KK_FLOAT diam = static_cast<KK_FLOAT>(2.0) * rad;

  KK_FLOAT rad2 = rad * rad;
  KK_FLOAT rad4 = rad2 * rad2;
  KK_FLOAT rad8 = rad4 * rad4;
  KK_FLOAT delta2 = rad2 - r * r;
  KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / delta2;
  KK_FLOAT r2inv = rinv * rinv;
  KK_FLOAT r4inv = r2inv * r2inv;
  KK_FLOAT r8inv = r4inv * r4inv;
  fwallKK = coeff1_kk *
          (rad8 * rad + static_cast<KK_FLOAT>(27.0) * rad4 * rad2 * rad * r * r + static_cast<KK_FLOAT>(63.0) * rad4 * rad * powint(r, 4) +
           static_cast<KK_FLOAT>(21.0) * rad2 * rad * powint(r, 6)) *
          r8inv -
      new_coeff2 * r2inv;

  KK_FLOAT r2 = static_cast<KK_FLOAT>(0.5) * diam - r;
  KK_FLOAT rinv2 = static_cast<KK_FLOAT>(1.0) / r2;
  KK_FLOAT r2inv2 = rinv2 * rinv2;
  KK_FLOAT r4inv2 = r2inv2 * r2inv2;
  KK_FLOAT r3 = r + static_cast<KK_FLOAT>(0.5) * diam;
  KK_FLOAT rinv3 = static_cast<KK_FLOAT>(1.0) / r3;
  KK_FLOAT r2inv3 = rinv3 * rinv3;
  KK_FLOAT r4inv3 = r2inv3 * r2inv3;
  return coeff3_kk *
          ((static_cast<KK_FLOAT>(-3.5) * diam + r) * r4inv2 * r2inv2 * rinv2 +
           (static_cast<KK_FLOAT>(3.5) * diam + r) * r4inv3 * r2inv3 * rinv3) -
      coeff4_kk * ((-diam * r + r2 * r3 * (Kokkos::log(-r2) - Kokkos::log(r3))) * (-rinv2) * rinv3) - offset_kk;
}

/* ----------------------------------------------------------------------
   harmonic interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT FixWallRegionKokkos<DeviceType>::harmonic(KK_FLOAT r, KK_FLOAT& fwallKK) const
{
  const KK_FLOAT cutoff_kk = static_cast<KK_FLOAT>(cutoff);
  const KK_FLOAT epsilon_kk = static_cast<KK_FLOAT>(epsilon);
  KK_FLOAT dr = cutoff_kk - r;
  fwallKK = static_cast<KK_FLOAT>(2.0) * epsilon_kk * dr;
  return epsilon_kk * dr * dr;
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallRegionKokkos<DeviceType>::v_tally(value_type result, int i, KK_FLOAT *v) const
{
  if (vflag_global) {
    result[4] += static_cast<double>(v[0]);
    result[5] += static_cast<double>(v[1]);
    result[6] += static_cast<double>(v[2]);
    result[7] += static_cast<double>(v[3]);
    result[8] += static_cast<double>(v[4]);
    result[9] += static_cast<double>(v[5]);
  }

  if (vflag_atom) {
    Kokkos::atomic_add(&(d_vatom(i,0)),static_cast<KK_ACC_FLOAT>(v[0]));
    Kokkos::atomic_add(&(d_vatom(i,1)),static_cast<KK_ACC_FLOAT>(v[1]));
    Kokkos::atomic_add(&(d_vatom(i,2)),static_cast<KK_ACC_FLOAT>(v[2]));
    Kokkos::atomic_add(&(d_vatom(i,3)),static_cast<KK_ACC_FLOAT>(v[3]));
    Kokkos::atomic_add(&(d_vatom(i,4)),static_cast<KK_ACC_FLOAT>(v[4]));
    Kokkos::atomic_add(&(d_vatom(i,5)),static_cast<KK_ACC_FLOAT>(v[5]));
  }
}

namespace LAMMPS_NS {
template class FixWallRegionKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallRegionKokkos<LMPHostType>;
#endif
}

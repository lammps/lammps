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

#include "fix_wall_piston_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallPistonKokkos<DeviceType>::FixWallPistonKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallPiston(lmp, narg, arg),
#ifndef LMP_KOKKOS_DEBUG_RNG
  rand_pool(tseed + comm->me)
#else
  rand_pool(tseed + comm->me, lmp)
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK;
  datamask_modify = X_MASK | V_MASK | F_MASK;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(randomt,tseed + comm->me);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallPistonKokkos<DeviceType>::~FixWallPistonKokkos()
{
  if (copymode) return;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallPistonKokkos<DeviceType>::init()
{
  FixWallPiston::init();

  // only the zlo piston has a device kernel; the base style already rejects
  // the other faces for the NL ramps, and the plain ramps only move zlo

  if (!zloflag)
    error->all(FLERR,"Fix wall/piston/kk is only implemented for the zlo wall");

  if (tempflag) {
    const int ntypes = atom->ntypes;
    k_gfactor1 = DAT::tdual_kkfloat_1d("wall/piston:gfactor1",ntypes+1);
    k_gfactor2 = DAT::tdual_kkfloat_1d("wall/piston:gfactor2",ntypes+1);
    d_gfactor1 = k_gfactor1.template view<DeviceType>();
    d_gfactor2 = k_gfactor2.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallPistonKokkos<DeviceType>::post_integrate()
{
  // all of the piston trajectory is scalar host work

  double zlo = z0;

  const double t = (update->ntimestep - update->beginstep) * update->dt;
  const double tott = (update->endstep - update->beginstep) * update->dt;
  const double tt = t * t;
  const double ttt = tt * t;
  const double tttt = tt * tt;
  const double t0p5 = sqrt(t/tott);
  const double t1p5 = t0p5*t0p5*t0p5;
  const double t2p5 = t1p5*t0p5*t0p5;

  if (rampflag) {
    paccelx = maxvx / tott;
    paccely = maxvy / tott;
    paccelz = maxvz / tott;
    zlo = z0 + 0.5 * paccelz * tt; vz = paccelz * t;
  } else if (rampNL1flag) {
    paccelz = maxvz / tott;
    angfreq = MY_2PI / (0.5 * tott);
    zlo = z0 + paccelz * (0.5*tt + 1.0/(angfreq*angfreq) -
                          1.0/(angfreq*angfreq)*cos(angfreq*t));
    vz = paccelz * (t + 1.0/angfreq*sin(angfreq*t));
  } else if (rampNL2flag) {
    paccelz = maxvz / tott;
    angfreq = 3.0*MY_2PI / tott;
    zlo = z0 + paccelz * (0.5*tt + 4.0/(3.0*angfreq*angfreq)*
                          (1.0-cos(angfreq*t)) +
                          1.0/(6.0*angfreq*angfreq)*(1.0-cos(2.0*angfreq*t)));
    vz = paccelz * (t + 4.0/(3.0*angfreq)*sin(angfreq*t) +
                    1.0/(3.0*angfreq)*sin(2.0*angfreq*t));
  } else if (rampNL3flag) {
    paccelz = maxvz / tott;
    zlo = z0 + paccelz*tott*tott/2.5 * t2p5;
    vz = paccelz * tott * t1p5;
  } else if (rampNL4flag) {
    paccelz = maxvz / tott;
    zlo = z0 + paccelz/tott/3.0 * ttt;
    vz = paccelz / tott * tt;
  } else if (rampNL5flag) {
    paccelz = maxvz / tott;
    zlo = z0 + paccelz/tott/tott/4.0 * tttt;
    vz = paccelz / tott / tott * ttt;
  } else {
    zlo = z0 + vz * t;
  }

  if ((update->ntimestep % 1000 == 0) && (comm->me == 0))
    utils::logmesg(lmp,"SHOCK: step {} t {} zpos {} vz {} az {} zlo {}\n",
                   update->ntimestep, t, zlo, vz, paccelz, domain->boxlo[2]);

  atomKK->sync(execution_space,datamask_read);

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  l_zlo = static_cast<KK_FLOAT>(zlo);
  l_vz = static_cast<KK_FLOAT>(vz);
  l_roughdist = static_cast<KK_FLOAT>(roughdist);
  for (int k = 0; k < 3; k++) {
    l_boxlo[k] = static_cast<KK_FLOAT>(domain->boxlo[k]);
    l_boxhi[k] = static_cast<KK_FLOAT>(domain->boxhi[k]);
  }

  const int nlocal = atom->nlocal;

  copymode = 1;
  if (roughflag)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallPistonReflect<1>>(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallPistonReflect<0>>(0,nlocal),*this);
  copymode = 0;

  // optional Langevin region ahead of the piston

  if (tempflag) {
    l_tsqrt = static_cast<KK_FLOAT>(sqrt(t_target));
    l_zcut = static_cast<KK_FLOAT>(domain->boxlo[2] + t_extent);

    if (atom->mass) {
      for (int i = 1; i <= atom->ntypes; i++) {
        gfactor1[i] = -atom->mass[i] / t_period / force->ftm2v;
        gfactor2[i] = sqrt(atom->mass[i]) *
          sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) / force->ftm2v;
        k_gfactor1.view_host()(i) = static_cast<KK_FLOAT>(gfactor1[i]);
        k_gfactor2.view_host()(i) = static_cast<KK_FLOAT>(gfactor2[i]);
      }
      k_gfactor1.template modify<LMPHostType>();
      k_gfactor2.template modify<LMPHostType>();
      k_gfactor1.template sync<DeviceType>();
      k_gfactor2.template sync<DeviceType>();
    } else {
      l_gamma1_pref = static_cast<KK_FLOAT>(-1.0 / t_period / force->ftm2v);
      l_gamma2_pref = static_cast<KK_FLOAT>(
        sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) / force->ftm2v);
    }

    copymode = 1;
    if (atom->mass)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallPistonTemp<0>>(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallPistonTemp<1>>(0,nlocal),*this);
    copymode = 0;
  }

  atomKK->modified(execution_space,datamask_modify);
}

/* ----------------------------------------------------------------------
   reflect atoms that have fallen behind the moving piston face
------------------------------------------------------------------------- */

template<class DeviceType>
template<int ROUGH>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallPistonKokkos<DeviceType>::operator()(TagFixWallPistonReflect<ROUGH>,
                                                 const int &i) const
{
  if (d_mask[i] & groupbit) {

    // roughoff is a local here: the base style keeps it in a member, which
    // would be a race between threads

    KK_FLOAT roughoff = 0.0;
    if (ROUGH) {
      roughoff += l_roughdist*Kokkos::fabs((d_x(i,0) - l_boxlo[0])/
                                           (l_boxhi[0]-l_boxlo[0])-static_cast<KK_FLOAT>(0.5));
      roughoff += l_roughdist*Kokkos::fabs((d_x(i,1) - l_boxlo[1])/
                                           (l_boxhi[1]-l_boxlo[1])-static_cast<KK_FLOAT>(0.5));
    }
    if (d_x(i,2) < l_zlo - roughoff) {
      d_x(i,2) = static_cast<KK_FLOAT>(2.0) * (l_zlo - roughoff) - d_x(i,2);
      d_v(i,2) = static_cast<KK_FLOAT>(2.0) * l_vz - d_v(i,2);
    }
  }
}

/* ----------------------------------------------------------------------
   Langevin thermostat applied to the slab just ahead of the piston
------------------------------------------------------------------------- */

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallPistonKokkos<DeviceType>::operator()(TagFixWallPistonTemp<RMASS>, const int &i) const
{
  if (d_x(i,2) > l_zcut) return;

  KK_FLOAT gamma1, gamma2;
  if (RMASS) {
    gamma1 = l_gamma1_pref * d_rmass(i);
    gamma2 = l_gamma2_pref * Kokkos::sqrt(d_rmass(i)) * l_tsqrt;
  } else {
    gamma1 = d_gfactor1(d_type(i));
    gamma2 = d_gfactor2(d_type(i)) * l_tsqrt;
  }

  rand_type rand_gen = rand_pool.get_state();
  const KK_FLOAT r0 = static_cast<KK_FLOAT>(rand_gen.drand() - 0.5);
  const KK_FLOAT r1 = static_cast<KK_FLOAT>(rand_gen.drand() - 0.5);
  const KK_FLOAT r2 = static_cast<KK_FLOAT>(rand_gen.drand() - 0.5);
  rand_pool.free_state(rand_gen);

  // the per-type branch of the base style damps v_z relative to the piston,
  // the per-atom-mass branch damps the lab-frame v_z

  d_f(i,0) += static_cast<KK_ACC_FLOAT>(gamma1*d_v(i,0) + gamma2*r0);
  d_f(i,1) += static_cast<KK_ACC_FLOAT>(gamma1*d_v(i,1) + gamma2*r1);
  if (RMASS) d_f(i,2) += static_cast<KK_ACC_FLOAT>(gamma1*d_v(i,2) + gamma2*r2);
  else d_f(i,2) += static_cast<KK_ACC_FLOAT>(gamma1*(d_v(i,2)-l_vz) + gamma2*r2);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallPistonKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallPistonKokkos<LMPHostType>;
#endif
}

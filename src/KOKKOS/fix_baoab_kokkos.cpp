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

#include "fix_baoab_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixBAOABKokkos<DeviceType>::FixBAOABKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixBAOAB(lmp, narg, arg),
#ifndef LMP_KOKKOS_DEBUG_RNG
  rand_pool(seed + comm->me)
#else
  rand_pool(seed + comm->me, lmp)
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK;
  datamask_modify = X_MASK | V_MASK;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(random,seed + comm->me);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixBAOABKokkos<DeviceType>::~FixBAOABKokkos()
{
  if (copymode) return;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBAOABKokkos<DeviceType>::init()
{
  FixBAOAB::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix baoab/kk");
}

/* ----------------------------------------------------------------------
   B-A-O-A: half force kick, half drift, exact Ornstein-Uhlenbeck thermostat,
   half drift, with an optional correction that removes the net random
   momentum the O step injected
------------------------------------------------------------------------- */

template<class DeviceType>
void FixBAOABKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space,datamask_read);

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  atomKK->k_mass.template sync<DeviceType>();

  // update target temperature (supports ramping via run start/stop)

  compute_target();

  // kT in velocity-squared * mass units

  const double kT = force->boltz * t_target / force->mvv2e;

  l_dtf = static_cast<KK_FLOAT>(dtf);
  l_dtby2 = static_cast<KK_FLOAT>(dtby2);
  l_c1 = static_cast<KK_FLOAT>(c1);
  l_kT = static_cast<KK_FLOAT>(kT);
  l_one_minus_c1sq = static_cast<KK_FLOAT>(1.0 - c1*c1);
  l_mvv2e = static_cast<KK_FLOAT>(force->mvv2e);

  const int rmass_flag = (atomKK->rmass) ? 1 : 0;

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // result = {fran_total[0..2], mass_total, energy_onestep}

  double result[5] = {0.0,0.0,0.0,0.0,0.0};

  copymode = 1;
  if (rmass_flag) {
    if (zeroflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixBAOABInitial<1,1>>(0,nlocal),*this,result);
    else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixBAOABInitial<1,0>>(0,nlocal),*this,result);
  } else {
    if (zeroflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixBAOABInitial<0,1>>(0,nlocal),*this,result);
    else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixBAOABInitial<0,0>>(0,nlocal),*this,result);
  }
  copymode = 0;

  energy_onestep = result[4];

  // zero net random momentum across all MPI ranks

  if (zeroflag) {
    double buf[4] = {result[0],result[1],result[2],result[3]};
    double bufall[4];
    MPI_Allreduce(buf,bufall,4,MPI_DOUBLE,MPI_SUM,world);

    if (bufall[3] > 0.0) {
      const double inv_mtot = 1.0 / bufall[3];
      for (int k = 0; k < 3; k++)
        l_vcorr[k] = static_cast<KK_FLOAT>(bufall[k] * inv_mtot);

      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixBAOABZeroMomentum>(0,nlocal),*this);
      copymode = 0;
    }
  }

  atomKK->modified(execution_space,datamask_modify);

  energy += energy_onestep;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBAOABKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space,V_MASK | F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  atomKK->k_mass.template sync<DeviceType>();

  l_dtf = static_cast<KK_FLOAT>(dtf);

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixBAOABFinal<1>>(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixBAOABFinal<0>>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int RMASS, int ZERO>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixBAOABKokkos<DeviceType>::operator()(TagFixBAOABInitial<RMASS,ZERO>, const int &i,
                                            value_type result) const
{
  if (!(d_mask[i] & groupbit)) return;

  const KK_FLOAT mi = RMASS ? d_rmass(i) : d_mass(d_type(i));
  const KK_FLOAT invmass = static_cast<KK_FLOAT>(1.0) / mi;

  // B step: half velocity kick from conservative forces

  const KK_FLOAT dtfm = l_dtf * invmass;
  d_v(i,0) += dtfm * static_cast<KK_FLOAT>(d_f(i,0));
  d_v(i,1) += dtfm * static_cast<KK_FLOAT>(d_f(i,1));
  d_v(i,2) += dtfm * static_cast<KK_FLOAT>(d_f(i,2));

  // A step: half position drift

  d_x(i,0) += l_dtby2 * d_v(i,0);
  d_x(i,1) += l_dtby2 * d_v(i,1);
  d_x(i,2) += l_dtby2 * d_v(i,2);

  // O step: exact Ornstein-Uhlenbeck thermostat over the full timestep

  const KK_FLOAT c2 = Kokkos::sqrt(l_kT * invmass * l_one_minus_c1sq);

  const KK_FLOAT ke_before = static_cast<KK_FLOAT>(0.5) * l_mvv2e * mi *
    (d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2));

  rand_type rand_gen = rand_pool.get_state();
  const KK_FLOAT r0 = static_cast<KK_FLOAT>(rand_gen.normal());
  const KK_FLOAT r1 = static_cast<KK_FLOAT>(rand_gen.normal());
  const KK_FLOAT r2 = static_cast<KK_FLOAT>(rand_gen.normal());
  rand_pool.free_state(rand_gen);

  d_v(i,0) = l_c1 * d_v(i,0) + c2 * r0;
  d_v(i,1) = l_c1 * d_v(i,1) + c2 * r1;
  d_v(i,2) = l_c1 * d_v(i,2) + c2 * r2;

  if (ZERO) {
    result[0] += static_cast<double>(mi * c2 * r0);
    result[1] += static_cast<double>(mi * c2 * r1);
    result[2] += static_cast<double>(mi * c2 * r2);
    result[3] += static_cast<double>(mi);
  }

  const KK_FLOAT ke_after = static_cast<KK_FLOAT>(0.5) * l_mvv2e * mi *
    (d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2));
  result[4] += static_cast<double>(ke_before - ke_after);

  // A step: second half position drift

  d_x(i,0) += l_dtby2 * d_v(i,0);
  d_x(i,1) += l_dtby2 * d_v(i,1);
  d_x(i,2) += l_dtby2 * d_v(i,2);
}

/* ----------------------------------------------------------------------
   remove the net random momentum, and undo its effect on the second A step
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixBAOABKokkos<DeviceType>::operator()(TagFixBAOABZeroMomentum, const int &i) const
{
  if (!(d_mask[i] & groupbit)) return;

  d_v(i,0) -= l_vcorr[0];
  d_v(i,1) -= l_vcorr[1];
  d_v(i,2) -= l_vcorr[2];

  d_x(i,0) -= l_dtby2 * l_vcorr[0];
  d_x(i,1) -= l_dtby2 * l_vcorr[1];
  d_x(i,2) -= l_dtby2 * l_vcorr[2];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixBAOABKokkos<DeviceType>::operator()(TagFixBAOABFinal<RMASS>, const int &i) const
{
  if (!(d_mask[i] & groupbit)) return;

  // B step: half velocity kick from the new forces

  const KK_FLOAT mi = RMASS ? d_rmass(i) : d_mass(d_type(i));
  const KK_FLOAT dtfm = l_dtf / mi;
  d_v(i,0) += dtfm * static_cast<KK_FLOAT>(d_f(i,0));
  d_v(i,1) += dtfm * static_cast<KK_FLOAT>(d_f(i,1));
  d_v(i,2) += dtfm * static_cast<KK_FLOAT>(d_f(i,2));
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixBAOABKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixBAOABKokkos<LMPHostType>;
#endif
}

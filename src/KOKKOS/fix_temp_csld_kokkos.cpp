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

#include "fix_temp_csld_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "kokkos.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixTempCSLDKokkos<DeviceType>::FixTempCSLDKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixTempCSLD(lmp, narg, arg),
#ifndef LMP_KOKKOS_DEBUG_RNG
  rand_pool(utils::inumeric(FLERR, arg[6], false, lmp) + comm->me)
#else
  rand_pool(utils::inumeric(FLERR, arg[6], false, lmp) + comm->me, lmp)
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(random, utils::inumeric(FLERR, arg[6], false, lmp) + comm->me);
#endif

  // the base class created the internal temperature compute relying on the
  // -sf kk suffix.  If this style is requested with an explicit /kk suffix but
  // without -sf kk, that compute is not a KOKKOS style and would force a
  // host/device sync every step.  Recreate it as temp/kk in that case.

  if (tflag) {
    Compute *c = modify->get_compute_by_id(id_temp);
    if (c && !c->kokkosable) {
      modify->delete_compute(id_temp);
      modify->add_compute(fmt::format("{} {} temp/kk", id_temp, group->names[igroup]));
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixTempCSLDKokkos<DeviceType>::~FixTempCSLDKokkos()
{
  if (copymode) return;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempCSLDKokkos<DeviceType>::init()
{
  FixTempCSLD::init();

  // the CSLD bias path mixes the bias of the *saved* velocities with freshly
  // drawn random velocities, which does not map onto remove_bias_all_kk
  // (that operates on atom->v).  Only the unbiased case is supported on device.

  if (which == BIAS)
    error->all(FLERR, "Fix temp/csld/kk does not support a biased temperature compute");

  KokkosLMP::warn_nonkokkos_compute(lmp, style, temperature, "temperature");
}

/* ----------------------------------------------------------------------
   compute the temperature scalar on the device when possible, otherwise via
   a host sync round-trip (correct but slow; e.g. fix_modify temp <non-kk>)
------------------------------------------------------------------------- */

template<class DeviceType>
double FixTempCSLDKokkos<DeviceType>::temp_scalar()
{
  if (temperature->kokkosable) return temperature->compute_scalar();

  atomKK->sync(temperature->execution_space, temperature->datamask_read);
  double t = temperature->compute_scalar();
  atomKK->modified(temperature->execution_space, temperature->datamask_modify);
  atomKK->sync(execution_space, temperature->datamask_modify);
  return t;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempCSLDKokkos<DeviceType>::end_of_step()
{
  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR, "Fix {} variable returned negative temperature", style);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  double t_current = temp_scalar();
  double ekin_old = t_current * 0.5 * temperature->dof * force->boltz;

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  // (re)allocate device holding space for the existing velocities

  if (nmax < atom->nlocal) {
    nmax = atom->nlocal + 1;
    d_vhold = typename AT::t_kkfloat_1d_3("csld/kk:vhold", nmax);
  }

  l_rmass_flag = (atom->rmass) ? 1 : 0;

  atomKK->sync(execution_space, V_MASK|MASK_MASK|TYPE_MASK|RMASS_MASK);
  atomKK->k_mass.sync<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  if (atom->rmass) rmass = atomKK->k_rmass.view<DeviceType>();
  else mass = atomKK->k_mass.view<DeviceType>();
  int nlocal = atom->nlocal;

  // The CSLD thermostat is a linear combination of old and new velocities,
  // where the new ones are randomly chosen from a gaussian distribution.
  // see Bussi and Parrinello, Phys. Rev. E (2007).
  // step 1: save current velocities, replace with random ones

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixTempCSLDInitial>(0,nlocal),*this);
  copymode = 0;
  atomKK->modified(execution_space, V_MASK);

  // mixing factors (c2 uses the temperature of the freshly drawn velocities)

  const double c1 = exp(-update->dt/t_period);
  const double c2 = sqrt((1.0-c1*c1)*t_target/temp_scalar());
  l_c1 = static_cast<KK_FLOAT>(c1);
  l_c2 = static_cast<KK_FLOAT>(c2);

  // step 2: mix the saved and random velocities

  v = atomKK->k_v.view<DeviceType>();
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixTempCSLDFinal>(0,nlocal),*this);
  copymode = 0;
  atomKK->modified(execution_space, V_MASK);

  // tally the kinetic energy transferred between heat bath and system

  t_current = temp_scalar();
  energy += ekin_old - t_current * 0.5 * temperature->dof * force->boltz;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixTempCSLDKokkos<DeviceType>::operator()(TagFixTempCSLDInitial, const int &i) const {
  if (mask[i] & groupbit) {
    KK_FLOAT m = l_rmass_flag ? rmass[i] : mass[type[i]];
    const KK_FLOAT factor = 1.0/Kokkos::sqrt(m);
    rand_type rand_gen = rand_pool.get_state();
    d_vhold(i,0) = v(i,0); v(i,0) = rand_gen.normal()*factor;
    d_vhold(i,1) = v(i,1); v(i,1) = rand_gen.normal()*factor;
    d_vhold(i,2) = v(i,2); v(i,2) = rand_gen.normal()*factor;
    rand_pool.free_state(rand_gen);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixTempCSLDKokkos<DeviceType>::operator()(TagFixTempCSLDFinal, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) = d_vhold(i,0)*l_c1 + v(i,0)*l_c2;
    v(i,1) = d_vhold(i,1)*l_c1 + v(i,1)*l_c2;
    v(i,2) = d_vhold(i,2)*l_c1 + v(i,2)*l_c2;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixTempCSLDKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixTempCSLDKokkos<LMPHostType>;
#endif
}

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

#include "fix_temp_csvr_kokkos.h"

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
#include "update.h"
#include "variable.h"

#include <cmath>
#include <type_traits>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixTempCSVRKokkos<DeviceType>::FixTempCSVRKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixTempCSVR(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  // the base class created the internal temperature compute relying on the
  // -sf kk suffix.  If this style is requested with an explicit /kk suffix but
  // without -sf kk, that compute is not a KOKKOS style and would force a
  // host/device sync every step.  Recreate it as temp/kk in that case.

  if (tflag) {
    Compute *c = modify->get_compute_by_id(id_temp);
    if (c && !c->kokkosable) {
      modify->delete_compute(id_temp);
      // match this fix's host/device flavor so a /kk/host fix gets a
      // host-space helper (the two are identical on CPU-only builds)
      const char *tempstyle =
        std::is_same_v<DeviceType,LMPDeviceType> ? "temp/kk" : "temp/kk/host";
      modify->add_compute(fmt::format("{} {} {}", id_temp, group->names[igroup], tempstyle));
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempCSVRKokkos<DeviceType>::init()
{
  FixTempCSVR::init();
  KokkosLMP::warn_nonkokkos_compute(lmp, style, temperature, "temperature");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempCSVRKokkos<DeviceType>::end_of_step()
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
      error->one(FLERR, "Fix {} variable {} returned negative temperature", style, tstr);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // current temperature: temp/kk reduces on the device and returns a scalar;
  // a non-kokkosable temperature needs an explicit sync (correct but slow)

  double t_current;
  if (temperature->kokkosable)
    t_current = temperature->compute_scalar();
  else {
    atomKK->sync(temperature->execution_space, temperature->datamask_read);
    t_current = temperature->compute_scalar();
    atomKK->modified(temperature->execution_space, temperature->datamask_modify);
    atomKK->sync(execution_space, temperature->datamask_modify);
  }

  const double efactor = 0.5 * temperature->dof * force->boltz;
  const double ekin_old = t_current * efactor;
  const double ekin_new = t_target * efactor;

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  // compute velocity scaling factor on root node and broadcast
  // (a handful of host RNG draws -- not per-atom)

  double lamda;
  if (comm->me == 0) lamda = resamplekin(ekin_old, ekin_new);
  MPI_Bcast(&lamda,1,MPI_DOUBLE,0,world);

  // remove velocity bias (group-wide) on the device if needed

  if (which == BIAS) {
    if (temperature->kokkosable) temperature->remove_bias_all_kk();
    else {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      temperature->remove_bias_all();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      atomKK->sync(execution_space, temperature->datamask_modify);
    }
  }

  // rescale velocities by lamda on the device

  atomKK->sync(execution_space, V_MASK|MASK_MASK);
  auto v = atomKK->k_v.view<DeviceType>();
  auto mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;
  auto groupbit = this->groupbit;

  const KK_FLOAT lamda_kk = static_cast<KK_FLOAT>(lamda);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal), LAMMPS_LAMBDA(int i) {
    if (mask[i] & groupbit) {
      v(i,0) *= lamda_kk;
      v(i,1) *= lamda_kk;
      v(i,2) *= lamda_kk;
    }
  });

  atomKK->modified(execution_space, V_MASK);

  if (which == BIAS) {
    if (temperature->kokkosable) temperature->restore_bias_all_kk();
    else {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      temperature->restore_bias_all();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      atomKK->sync(execution_space, temperature->datamask_modify);
    }
  }

  // tally the kinetic energy transferred between heat bath and system

  energy += ekin_old * (1.0 - lamda*lamda);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixTempCSVRKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixTempCSVRKokkos<LMPHostType>;
#endif
}

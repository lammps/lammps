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

#include "fix_temp_rescale_kokkos.h"

#include "atom_kokkos.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"
#include "atom_masks.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixTempRescaleKokkos<DeviceType>::FixTempRescaleKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixTempRescale(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempRescaleKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(temperature->execution_space,temperature->datamask_read);
  double t_current = temperature->compute_scalar();
  atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  atomKK->sync(execution_space,temperature->datamask_modify);

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  // protect against division by zero

  if (t_current == 0.0)
    error->all(FLERR,"Computed temperature for fix temp/rescale cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR, "Fix temp/rescale variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocity of appropriate atoms if outside window
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since factor is multiplied by v

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    energy += (t_current-t_target) * efactor;

    auto v = atomKK->k_v.view<DeviceType>();
    auto mask = atomKK->k_mask.view<DeviceType>();
    int nlocal = atom->nlocal;
    auto groupbit = this->groupbit;

    if (which == NOBIAS) {
      atomKK->sync(temperature->execution_space,temperature->datamask_read);
      temperature->remove_bias_all();
      atomKK->modified(temperature->execution_space,temperature->datamask_modify);
      atomKK->sync(execution_space,temperature->datamask_modify);
    }

    atomKK->sync(execution_space,V_MASK|MASK_MASK);

    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal), LAMMPS_LAMBDA(int i) {
      if (mask[i] & groupbit) {
        v(i,0) *= factor;
        v(i,1) *= factor;
        v(i,2) *= factor;
      }
    });

    atomKK->modified(execution_space,V_MASK);

    if (which == NOBIAS) {
      atomKK->sync(temperature->execution_space,temperature->datamask_read);
      temperature->restore_bias_all();
      atomKK->modified(temperature->execution_space,temperature->datamask_modify);
      atomKK->sync(execution_space,temperature->datamask_modify);

    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixTempRescaleKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixTempRescaleKokkos<LMPHostType>;
#endif
}

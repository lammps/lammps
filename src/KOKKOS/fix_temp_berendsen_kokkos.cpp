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

#include "fix_temp_berendsen_kokkos.h"

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
FixTempBerendsenKokkos<DeviceType>::FixTempBerendsenKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixTempBerendsen(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempBerendsenKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(temperature->execution_space,temperature->datamask_read);
  double t_current = temperature->compute_scalar();
  atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  atomKK->sync(execution_space,temperature->datamask_modify);

  double tdof = temperature->dof;

  // there is nothing to do, if there are no degrees of freedom

  if (tdof < 1) return;

  if (t_current == 0.0)
    error->all(FLERR, "Computed current temperature for fix temp/berendsen must not be 0.0");

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
      error->one(FLERR, "Fix temp/berendsen variable {} returned negative temperature",
                 input->variable->names[tvar]);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocities by lamda
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since lamda is multiplied by v

  double lamda = sqrt(1.0 + update->dt/t_period*(t_target/t_current - 1.0));
  double efactor = 0.5 * force->boltz * tdof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

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
      v(i,0) *= lamda;
      v(i,1) *= lamda;
      v(i,2) *= lamda;
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
/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixTempBerendsenKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixTempBerendsenKokkos<LMPHostType>;
#endif
}

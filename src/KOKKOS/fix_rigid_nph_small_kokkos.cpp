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
   Contributing author: Trung Dac Nguyen (ORNL) for the host fix rigid/nph/small
------------------------------------------------------------------------- */

#include "fix_rigid_nph_small_kokkos.h"

#include "error.h"
#include "modify.h"

#include <string>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNPHSmallKokkos<DeviceType>::FixRigidNPHSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNHSmallKokkos<DeviceType>(lmp, narg, arg)
{
  // other setting are made by parent

  this->scalar_flag = 1;
  this->restart_global = 1;
  this->extscalar = 1;

  // error checks

  if (this->pstat_flag == 0)
    this->error->all(FLERR,"Pressure control must be used with fix nph/small");
  if (this->tstat_flag == 1)
    this->error->all(FLERR,"Temperature control must not be used with fix nph/small");
  if (this->p_start[0] < 0.0 || this->p_start[1] < 0.0 || this->p_start[2] < 0.0 ||
      this->p_stop[0] < 0.0 || this->p_stop[1] < 0.0 || this->p_stop[2] < 0.0)
    this->error->all(FLERR,"Target pressure for fix rigid/nph cannot be < 0.0");

  // convert input periods to frequency

  this->p_freq[0] = this->p_freq[1] = this->p_freq[2] = 0.0;

  if (this->p_flag[0]) this->p_freq[0] = 1.0 / this->p_period[0];
  if (this->p_flag[1]) this->p_freq[1] = 1.0 / this->p_period[1];
  if (this->p_flag[2]) this->p_freq[2] = 1.0 / this->p_period[2];

  // create a new compute temp style (group all)

  // create the temperature compute as temp/kk: Modify::add_compute only applies
  // the /kk suffix when -sf kk is active, so requesting this style with an
  // explicit suffix would otherwise get a host compute and force a full x/v
  // host<->device transfer every step (see the KOKKOS package instructions,
  // "Internal helper computes/fixes must be KOKKOS too").  compute pressure has
  // no /kk variant, so only the temperature is promoted.
  this->id_temp = utils::strdup(std::string(this->id)+"_temp");
  this->modify->add_compute(fmt::format("{} all temp/kk",this->id_temp));
  this->tcomputeflag = 1;

  // create a new compute pressure style (group all), pass id_temp as 4th arg

  this->id_press = utils::strdup(std::string(this->id)+"_press");
  this->modify->add_compute(fmt::format("{} all pressure {}",this->id_press,this->id_temp));
  this->pcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidNPHSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidNPHSmallKokkos<LMPHostType>;
#endif
}

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
   Contributing author: Trung Dac Nguyen (ORNL) for the host fix rigid/npt/small
------------------------------------------------------------------------- */

#include "fix_rigid_npt_small_kokkos.h"

#include "error.h"
#include "modify.h"

#include <string>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNPTSmallKokkos<DeviceType>::FixRigidNPTSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNHSmallKokkos<DeviceType>(lmp, narg, arg)
{
  // other setting are made by parent

  this->scalar_flag = 1;
  this->restart_global = 1;
  this->extscalar = 1;

  // error checks

  if (this->tstat_flag == 0 || this->pstat_flag == 0)
    this->error->all(FLERR,"Did not set temp or press for fix rigid/npt/small");
  if (this->t_start <= 0.0 || this->t_stop <= 0.0)
    this->error->all(FLERR,"Target temperature for fix rigid/npt/small cannot be 0.0");
  if (this->p_start[0] < 0.0 || this->p_start[1] < 0.0 || this->p_start[2] < 0.0 ||
      this->p_stop[0] < 0.0 || this->p_stop[1] < 0.0 || this->p_stop[2] < 0.0)
    this->error->all(FLERR,"Target pressure for fix rigid/npt/small cannot be < 0.0");

  if (this->t_period <= 0.0) this->error->all(FLERR,"Fix rigid/npt/small period must be > 0.0");

  // thermostat chain parameters

  if (this->t_chain < 1) this->error->all(FLERR,"Fix rigid npt/small t_chain should not be less than 1");
  if (this->t_iter < 1) this->error->all(FLERR,"Fix rigid npt/small t_chain should not be less than 1");
  if (this->t_order != 3 && this->t_order != 5)
    this->error->all(FLERR,"Fix rigid npt/small t_order must be 3 or 5");

  // convert input periods to frequency

  this->t_freq = 0.0;
  this->p_freq[0] = this->p_freq[1] = this->p_freq[2] = 0.0;

  this->t_freq = 1.0 / this->t_period;
  if (this->p_flag[0]) this->p_freq[0] = 1.0 / this->p_period[0];
  if (this->p_flag[1]) this->p_freq[1] = 1.0 / this->p_period[1];
  if (this->p_flag[2]) this->p_freq[2] = 1.0 / this->p_period[2];

  // create a new compute temp style (group all)

  this->id_temp = utils::strdup(std::string(this->id)+"_temp");
  this->modify->add_compute(fmt::format("{} all temp",this->id_temp));
  this->tcomputeflag = 1;

  // create a new compute pressure style (group all), pass id_temp as 4th arg

  this->id_press = utils::strdup(std::string(this->id)+"_press");
  this->modify->add_compute(fmt::format("{} all pressure {}",this->id_press,this->id_temp));
  this->pcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidNPTSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidNPTSmallKokkos<LMPHostType>;
#endif
}

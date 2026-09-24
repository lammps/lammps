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
   Contributing author: Trung Dac Nguyen (ORNL) for the host fix rigid/nvt/small
------------------------------------------------------------------------- */

#include "fix_rigid_nvt_small_kokkos.h"

#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNVTSmallKokkos<DeviceType>::FixRigidNVTSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNHSmallKokkos<DeviceType>(lmp, narg, arg)
{
  // other settings are made by parent

  this->scalar_flag = 1;
  this->restart_global = 1;
  this->extscalar = 1;

  // error checking; convert input period to frequency

  if (this->tstat_flag == 0)
    this->error->all(FLERR,"Did not set temp for fix rigid/nvt/small");
  if (this->t_start < 0.0 || this->t_stop <= 0.0)
    this->error->all(FLERR,"Target temperature for fix rigid/nvt/small cannot be 0.0");
  if (this->t_period <= 0.0) this->error->all(FLERR,"Fix rigid/nvt/small period must be > 0.0");
  this->t_freq = 1.0 / this->t_period;

  if (this->t_chain < 1) this->error->all(FLERR,"Fix rigid nvt/small t_chain should not be less than 1");
  if (this->t_iter < 1) this->error->all(FLERR,"Fix rigid nvt/small t_iter should not be less than 1");
  if (this->t_order != 3 && this->t_order != 5)
    this->error->all(FLERR,"Fix rigid nvt/small t_order must be 3 or 5");
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidNVTSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidNVTSmallKokkos<LMPHostType>;
#endif
}

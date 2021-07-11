// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_npt_kokkos.h"

#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNPTKokkos<DeviceType>::FixNPTKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNHKokkos<DeviceType>(lmp, narg, arg)
{
  this->kokkosable = 1;
  if (!this->tstat_flag)
    this->error->all(FLERR,"Temperature control must be used with fix npt/kk");
  if (!this->pstat_flag)
    this->error->all(FLERR,"Pressure control must be used with fix npt/kk");

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  // and thus its KE/temperature contribution should use group all

  this->id_temp = utils::strdup(std::string(this->id) + "_temp");
  this->modify->add_compute(fmt::format("{} all temp/kk",this->id_temp));
  this->tcomputeflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  this->id_press = utils::strdup(std::string(this->id) + "_press");
  this->modify->add_compute(fmt::format("{} all pressure {}",
                                        this->id_press, this->id_temp));
  this->pcomputeflag = 1;
}

namespace LAMMPS_NS {
template class FixNPTKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNPTKokkos<LMPHostType>;
#endif
}


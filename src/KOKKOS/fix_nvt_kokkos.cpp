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

#include "fix_nvt_kokkos.h"

#include "error.h"
#include "group.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVTKokkos<DeviceType>::FixNVTKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNHKokkos<DeviceType>(lmp, narg, arg)
{
  this->kokkosable = 1;
  if (!this->tstat_flag)
    this->error->all(FLERR,"Temperature control must be used with fix nvt/kk");
  if (this->pstat_flag)
    this->error->all(FLERR,"Pressure control can not be used with fix nvt/kk");

  // create a new compute temp style
  // id = fix-ID + temp

  this->id_temp = utils::strdup(std::string(this->id)+"_temp");
  this->modify->add_compute(fmt::format("{} {} temp/kk",this->id_temp,this->group->names[this->igroup]));
  this->tcomputeflag = 1;
}

namespace LAMMPS_NS {
template class FixNVTKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVTKokkos<LMPHostType>;
#endif
}


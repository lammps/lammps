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

#include "fix_nvt_sphere_kokkos.h"

#include "error.h"
#include "group.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixNVTSphereKokkos<DeviceType>::FixNVTSphereKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixNHSphereKokkos<DeviceType>(lmp, narg, arg)
{
  if (!this->tstat_flag)
    this->error->all(FLERR, "Temperature control must be used with fix nvt/sphere/kk");
  if (this->pstat_flag)
    this->error->all(FLERR, "Pressure control can not be used with fix nvt/sphere/kk");

  // create a new compute temp/sphere style
  // id = fix-ID + temp

  this->id_temp = utils::strdup(std::string(this->id) + "_temp");
  this->modify->add_compute(
      fmt::format("{} {} temp/sphere/kk", this->id_temp,
                  this->group->names[this->igroup]));
  this->tcomputeflag = 1;
}

namespace LAMMPS_NS {
template class FixNVTSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVTSphereKokkos<LMPHostType>;
#endif
}

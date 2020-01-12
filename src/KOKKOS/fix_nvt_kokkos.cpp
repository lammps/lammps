/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nvt_kokkos.h"
#include <cstring>
#include "group.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVTKokkos<DeviceType>::FixNVTKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNHKokkos<DeviceType>(lmp, narg, arg)
{
  this->kokkosable = 1;
  if (!this->tstat_flag)
    this->error->all(FLERR,"Temperature control must be used with fix nvt");
  if (this->pstat_flag)
    this->error->all(FLERR,"Pressure control can not be used with fix nvt");

  // create a new compute temp style
  // id = fix-ID + temp

  int n = strlen(this->id) + 6;
  this->id_temp = new char[n];
  strcpy(this->id_temp,this->id);
  strcat(this->id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = this->id_temp;
  newarg[1] = this->group->names[this->igroup];
  newarg[2] = (char *) "temp/kk";

  this->modify->add_compute(3,newarg);
  delete [] newarg;
  this->tcomputeflag = 1;
}

namespace LAMMPS_NS {
template class FixNVTKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixNVTKokkos<LMPHostType>;
#endif
}


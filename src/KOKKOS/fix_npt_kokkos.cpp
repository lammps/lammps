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

#include <cstring>
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
    this->error->all(FLERR,"Temperature control must be used with fix npt");
  if (!this->pstat_flag)
    this->error->all(FLERR,"Pressure control must be used with fix npt");

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  // and thus its KE/temperature contribution should use group all

  int n = strlen(this->id) + 6;
  this->id_temp = new char[n];
  strcpy(this->id_temp,this->id);
  strcat(this->id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = this->id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp/kk";

  this->modify->add_compute(3,newarg);
  delete [] newarg;
  this->tcomputeflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(this->id) + 7;
  this->id_press = new char[n];
  strcpy(this->id_press,this->id);
  strcat(this->id_press,"_press");

  newarg = new char*[4];
  newarg[0] = this->id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = this->id_temp;
  this->modify->add_compute(4,newarg);
  delete [] newarg;
  this->pcomputeflag = 1;
}

namespace LAMMPS_NS {
template class FixNPTKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixNPTKokkos<LMPHostType>;
#endif
}


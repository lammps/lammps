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

#include "fix_nvt_gpu.h"

#include "error.h"
#include "group.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVTGPU::FixNVTGPU(LAMMPS *lmp, int narg, char **arg) :
  FixNHGPU(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/gpu");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt/gpu");

  // create a new compute temp style
  // id = fix-ID + temp

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} temp",id_temp,group->names[igroup]));
  tcomputeflag = 1;
}

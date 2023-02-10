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

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#include "fix_nvt_uef.h"
#include "error.h"

using namespace LAMMPS_NS;

FixNVTUef::FixNVTUef(LAMMPS *lmp, int narg, char **arg) :
  FixNHUef(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/uef");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can't be used with fix nvt/uef");
}



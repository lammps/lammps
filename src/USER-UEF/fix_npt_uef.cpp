/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#include "fix_npt_uef.h"
#include "error.h"

using namespace LAMMPS_NS;

FixNPTUef::FixNPTUef(LAMMPS *lmp, int narg, char **arg) :
  FixNHUef(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix npt/uef");
  if (!pstat_flag)
    error->all(FLERR,"Pressure control must be used with fix npt/uef");
}

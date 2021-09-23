/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Taylor Barnes (MolSSI)
   MolSSI Driver Interface (MDI) support for LAMMPS
------------------------------------------------------------------------- */

#include "error.h"
#include "fix_mdi_engine2.h"
#include "mdi_engine2.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMDIEngine2::FixMDIEngine2(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal fix mdi/engine command");
}

/* ---------------------------------------------------------------------- */

int FixMDIEngine2::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine2::post_integrate()
{
  mdi_engine->engine_node("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine2::min_pre_force(int vflag)
{
  mdi_engine->engine_node("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine2::post_force(int vflag)
{
  mdi_engine->engine_node("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine2::min_post_force(int vflag)
{
  mdi_engine->engine_node("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine2::end_of_step()
{
  mdi_engine->engine_node("@ENDSTEP");
}

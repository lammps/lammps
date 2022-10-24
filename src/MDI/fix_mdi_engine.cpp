/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "fix_mdi_engine.h"

#include "error.h"
#include "update.h"

#include "mdi_engine.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMDIEngine::FixMDIEngine(LAMMPS *_lmp, int narg, char **arg) : Fix(_lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal fix mdi/engine command");
}

/* ---------------------------------------------------------------------- */

int FixMDIEngine::setmask()
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

void FixMDIEngine::setup(int /*vflag*/)
{
  // engine is now at FORCES node

  mdi_engine->engine_node("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_integrate()
{
  // engine is now at COORDS node for MD

  mdi_engine->engine_node("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_pre_force(int /*vflag*/)
{
  // engine is now at COORDS node for minimizer

  mdi_engine->engine_node("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_force(int /*vflag*/)
{
  // engine is now at FORCES node for MD

  mdi_engine->engine_node("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_post_force(int /*vflag*/)
{
  // engine is now at FORCES node for minimizer

  mdi_engine->engine_node("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::end_of_step()
{
  // engine is now at ENDSTEP node for MD

  mdi_engine->engine_node("@ENDSTEP");
}

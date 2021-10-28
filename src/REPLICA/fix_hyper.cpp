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

#include "fix_hyper.h"
#include "update.h"
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixHyper::FixHyper(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  hyperflag = 0;
  ntimestep_initial = update->ntimestep;
}

/* ----------------------------------------------------------------------
   extract hyper flag setting for all Fixes that perform hyperdynamics
------------------------------------------------------------------------- */

void *FixHyper::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"hyperflag") == 0) {
    return &hyperflag;
  }
  return nullptr;
}


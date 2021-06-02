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

#include "compute_deprecated.h"

#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDeprecated::ComputeDeprecated(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  std::string my_style = style;

  if (my_style == "DEPRECATED") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\nCompute style 'DEPRECATED' is a dummy style\n\n");
    return;
  }
  error->all(FLERR,"This compute style is no longer available");
}

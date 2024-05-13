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

#include "compute_deprecated.h"

#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDeprecated::ComputeDeprecated(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  std::string my_style = style;

  if (my_style == "DEPRECATED") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp, "\nCompute style 'DEPRECATED' is a dummy style\n\n");
    return;
  } else if (my_style == "mesont") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,
                     "\nCompute style 'mesont' and the associated pair style have been "
                     "removed. Please use pair style 'mesocnt' instead.\n\n");
  }
  error->all(FLERR, "This compute style is no longer available");
}

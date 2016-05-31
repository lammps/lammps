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

#include "fix_ave_spatial.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAveSpatial::FixAveSpatial(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  const char *message = "\n"
    "NOTE: The fix ave/spatial command has been replaced by the\n"
    "more general fix ave/chunk and compute chunk/atom commands.\n"
    "All ave/spatial functionality is available in the new commands.\n"
    "These ave/spatial keywords & options are part of fix ave/chunk:\n"
    "  Nevery, Nrepeat, Nfreq, input values, norm, ave, file, overwrite, title123\n"
    "These ave/spatial keywords & options for binning are part of compute chunk/atom:\n"
    "  dim, origin, delta, region, bound, discard, units\n\n";

  if (comm->me == 0) {
    if (screen) fprintf(screen,message);
    if (logfile) fprintf(logfile,message);
  }

  error->all(FLERR,"The fix ave/spatial command has been removed from LAMMPS");
}

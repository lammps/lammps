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
#include "fix_deprecated.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

static void writemsg(LAMMPS *lmp, const char *msg, int abend=1)
{
  if (lmp->comm->me == 0) {
    if (lmp->screen) fputs(msg,lmp->screen);
    if (lmp->logfile) fputs(msg,lmp->logfile);
  }
  if (abend)
    lmp->error->all(FLERR,"This fix style is no longer available");
}

/* ---------------------------------------------------------------------- */

FixDeprecated::FixDeprecated(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"DEPRECATED") == 0) {
    writemsg(lmp,"\nFix style 'DEPRECATED' is a dummy style\n\n",0);

  } else if (strncmp(style,"ave/spatial",11) == 0) {
    writemsg(lmp,"\nFix styles 'ave/spatial' and 'ave/spatial/sphere' have "
             "been replaced\nby the more general fix ave/chunk and compute "
             "chunk/atom commands.\nAll ave/spatial and ave/spatial/sphere "
             "functionality is available in these\nnew commands. These "
             "ave/spatial keywords & options are part of fix ave/chunk:\n"
             "  Nevery, Nrepeat, Nfreq, input values, norm, ave, file, "
             "overwrite, title123\nThese ave/spatial keywords & options for "
             "binning are part of compute chunk/atom:\n  dim, origin, delta,"
             " region, bound, discard, units\n\n");
  }
}

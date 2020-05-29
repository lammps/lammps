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

#include "dump_deprecated.h"
#include <cstring>
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
    lmp->error->all(FLERR,"This dump style is no longer available");
}

/* ---------------------------------------------------------------------- */

DumpDeprecated::DumpDeprecated(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (strcmp(style,"DEPRECATED") == 0) {
    writemsg(lmp,"\nDump style 'DEPRECATED' is a dummy style\n\n",0);

  }
}

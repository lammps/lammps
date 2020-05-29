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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "bond_deprecated.h"
#include <cstring>
#include "bond_hybrid.h"
#include "comm.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

static void writemsg(LAMMPS *lmp, const char *msg, int abend=1)
{
  if (lmp->comm->me == 0) {
    if (lmp->screen) fputs(msg,lmp->screen);
    if (lmp->logfile) fputs(msg,lmp->logfile);
  }
  if (abend)
    lmp->error->all(FLERR,"This bond style is no longer available");
}

/* ---------------------------------------------------------------------- */

void BondDeprecated::settings(int, char **)
{
  const char *my_style = force->bond_style;

  // hybrid substyles are created in BondHybrid::settings(), so when this is
  // called, our style was just added at the end of the list of substyles

  if (strncmp(my_style,"hybrid",6) == 0) {
    BondHybrid *hybrid = (BondHybrid *)force->bond;
    my_style = hybrid->keywords[hybrid->nstyles];
  }

  if (strcmp(my_style,"DEPRECATED") == 0) {
    writemsg(lmp,"\nBond style 'DEPRECATED' is a dummy style\n\n",0);

  }
}



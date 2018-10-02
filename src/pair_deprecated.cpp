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

#include "pair_deprecated.h"
#include "pair_hybrid.h"
#include "comm.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairDeprecated::PairDeprecated(LAMMPS *lmp) : Pair(lmp)
{
  if (strcmp(force->pair_style,"deprecated") == 0) {
    const char *message = "\n"
    "NOTE: The pair style 'deprecated' is a dummy fix style that was added to\n"
    "LAMMPS in order to print suitable error messages for deleted features.\n\n";

    if (comm->me == 0) {
      if (screen) fputs(message,screen);
      if (logfile) fputs(message,logfile);
    }
  }
  if (strncmp(force->pair_style,"reax",11) == 0) {
    const char *message = "\n"
    "NOTE: The pair style 'reax' has been removed from LAMMPS after the\n"
    "## November 2018 stable release. Its functionality has long before\n"
    "been superseded by pair styles 'reax/c' and 'reax/c/kk'\n\n";

    if (comm->me == 0) {
      if (screen) fputs(message,screen);
      if (logfile) fputs(message,logfile);
    }
  }
  error->all(FLERR,"This pair_style command has been removed from LAMMPS");
}


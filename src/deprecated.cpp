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
   Contributing authors:  Axel Kohlmeyer (Temple U),
------------------------------------------------------------------------- */

#include "deprecated.h"
#include <string>
#include "comm.h"
#include "error.h"
#include "input.h"
#include "utils.h"

using namespace LAMMPS_NS;

static void writemsg(LAMMPS *lmp, const std::string &msg, int abend=1)
{
  if (lmp->comm->me == 0) utils::logmesg(lmp,msg);
  if (abend)
    lmp->error->all(FLERR,"This command is no longer available");
}

/* ---------------------------------------------------------------------- */

void Deprecated::command(int /* narg */, char ** /* arg */)
{
  const std::string cmd = input->command;

  if (cmd == "DEPRECATED") {
    writemsg(lmp,"\nCommand 'DEPRECATED' is a dummy command\n\n",0);
  } else if (cmd == "reset_ids") {
    writemsg(lmp,"\n'reset_ids' has been renamed to 'reset_atom_ids'\n\n");
  }
}

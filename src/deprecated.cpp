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

/* ---------------------------------------------------------------------- */

void Deprecated::command(int /* narg */, char ** /* arg */)
{
  const std::string cmd = input->command;

  if (cmd == "DEPRECATED") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\nCommand 'DEPRECATED' is a dummy command\n\n");
    return;
  } else if (cmd == "reset_ids") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\n'reset_ids' has been renamed to 'reset_atom_ids'\n\n");
  }
  error->all(FLERR,"This command is no longer available");
}

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

/* ----------------------------------------------------------------------
   Contributing authors:  Axel Kohlmeyer (Temple U),
------------------------------------------------------------------------- */

#include "deprecated.h"

#include "comm.h"
#include "error.h"
#include "input.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void Deprecated::command(int narg, char **arg)
{
  const std::string cmd = input->command;

  if (cmd == "DEPRECATED") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\nCommand 'DEPRECATED' is a dummy command\n\n");
    return;
  } else if (cmd == "reset_ids") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\n'reset_ids' has been renamed to 'reset_atom_ids'\n\n");
  } else if (utils::strmatch(cmd,"^kim_")) {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\nWARNING: 'kim_<command>' has been renamed to "
                      "'kim <command>'. Please update your input.\n\n");
    std::string newcmd("kim");
    newcmd += " " + cmd.substr(4);
    for (int i=0; i < narg; ++i) {
       newcmd.append(1,' ');
       newcmd.append(arg[i]);
    }
    input->one(newcmd);
    return;
  }
  error->all(FLERR,"This command is no longer available");
}

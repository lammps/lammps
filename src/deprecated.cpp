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
    if (lmp->comm->me == 0) utils::logmesg(lmp, "\nCommand 'DEPRECATED' is a dummy command\n\n");
    return;
  } else if (cmd == "box") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp, "\nThe 'box' command has been removed and will be ignored\n\n");
    return;
  } else if (utils::strmatch(cmd, "^kim_")) {
    std::string newcmd("kim");
    newcmd += " " + cmd.substr(4);
    if (lmp->comm->me == 0)
      utils::logmesg(lmp, "\nWARNING: '{}' has been renamed to '{}'. Please update your input.\n\n",
                     cmd, newcmd);
    for (int i = 0; i < narg; ++i) {
      newcmd.append(1, ' ');
      newcmd.append(arg[i]);
    }
    input->one(newcmd);
    return;
  } else if (utils::strmatch(cmd, "^reset_")) {
    std::string newcmd("reset_atoms");
    if ((cmd == "reset_ids") || (cmd == "reset_atom_ids")) newcmd += " id";
    if (cmd == "reset_mol_ids") newcmd += " mol";
    if (lmp->comm->me == 0)
      utils::logmesg(lmp, "\nWARNING: '{}' has been renamed to '{}'. Please update your input.\n\n",
                     cmd, newcmd);
    for (int i = 0; i < narg; ++i) {
      newcmd.append(1, ' ');
      newcmd.append(arg[i]);
    }
    input->one(newcmd);
    return;
  } else if ((cmd == "message") || (cmd == "server")) {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp, "\nThe MESSAGE package has been replaced by the MDI package.\n\n");
  }
  error->all(FLERR, "This command is no longer available");
}

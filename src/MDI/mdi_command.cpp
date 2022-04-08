/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mdi_command.h"

#include "error.h"
#include "mdi_engine.h"
#include "mdi_plugin.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   mdi command: engine or plugin
---------------------------------------------------------------------- */

void MDICommand::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal mdi command");

  if (strcmp(arg[0], "engine") == 0) {
    MDIEngine(lmp, narg - 1, &arg[1]);
  } else if (strcmp(arg[0], "plugin") == 0) {
    MDIPlugin(lmp, narg - 1, &arg[1]);
  } else
    error->all(FLERR, "Illegal mdi command");
}

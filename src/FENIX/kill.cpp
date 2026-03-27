/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "kill.h"
#include "universe.h"
#include "comm.h"
#include "error.h"
#include "utils.h"

#include <signal.h>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Kill::Kill(LAMMPS* lmp_in) : Command(lmp), lmp(lmp_in) {};

/* ---------------------------------------------------------------------- */

void Kill::command(int narg, char** arg) {
  int kill_rank = -1;
  if(narg < 1) utils::missing_cmd_args(FLERR, "kill", lmp->error);

  for(int i = 0; i < narg; i++){
    if(strcmp(arg[i], "rank") == 0){
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "kill rank", lmp->error);
      kill_rank = utils::inumeric(FLERR, arg[++i], false, lmp);
    } else {
      lmp->error->all(FLERR, "Invalid argument kill {}", arg[i]);
    }
  }

  if(lmp->comm->me == kill_rank) kill();
}

/* ---------------------------------------------------------------------- */

[[noreturn]] /*static*/ void Kill::kill(){
  raise(SIGTERM);
  std::abort();
}

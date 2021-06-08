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

#include "server.h"
#include <cstring>
#include "error.h"

// customize by adding a new server protocol include and enum

#include "server_md.h"
#include "server_mc.h"

using namespace LAMMPS_NS;

enum{MD,MC};

/* ---------------------------------------------------------------------- */

void Server::command(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal server command");

  if (lmp->clientserver != 2)
    error->all(FLERR,"Message command not used to setup LAMMPS as a server");

  // customize by adding a new server protocol

  int protocol;
  if (strcmp(arg[0],"md") == 0) protocol = MD;
  else if (strcmp(arg[0],"mc") == 0) protocol = MC;
  else error->all(FLERR,"Unknown message protocol");

  if (protocol == MD) {
    ServerMD *server = new ServerMD(lmp);
    server->loop();
  } else if (protocol == MC) {
    ServerMC *server = new ServerMC(lmp);
    server->loop();
  }
}

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

#include "message.h"
#include <cstring>
#include "error.h"

// CSlib interface

#include "cslib.h"

using namespace LAMMPS_NS;
using namespace CSLIB_NS;

/* ---------------------------------------------------------------------- */

void Message::command(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal message command");

  int clientserver;
  if (strcmp(arg[0],"client") == 0) clientserver = 1;
  else if (strcmp(arg[0],"server") == 0) clientserver = 2;
  else error->all(FLERR,"Illegal message command");
  lmp->clientserver = clientserver;

  // validate supported protocols

  if ((strcmp(arg[1],"md") != 0) && (strcmp(arg[1],"mc") != 0))
    error->all(FLERR,"Unknown message protocol");

  // instantiate CSlib with chosen communication mode

  if (strcmp(arg[2],"file") == 0 || strcmp(arg[2],"zmq") == 0 ||
      strcmp(arg[2],"mpi/two") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal message command");
    lmp->cslib = new CSlib(clientserver-1,arg[2],arg[3],&world);

  } else if (strcmp(arg[2],"mpi/one") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal message command");
    if (!lmp->cscomm)
      error->all(FLERR,"Message mpi/one mode, but -mpi cmdline arg not used");
    lmp->cslib = new CSlib(clientserver-1,arg[2],&lmp->cscomm,&world);

  } else error->all(FLERR,"Illegal message command");

  // perform initial handshake between client and server
  // other code being coupled to must perform similar operation
  // client sends protocol with msgID = 0
  // server matches it and replies

  CSlib *cs = (CSlib *) lmp->cslib;

  if (clientserver == 1) {
    cs->send(0,1);
    cs->pack_string(1,arg[1]);

    int nfield;
    int *fieldID,*fieldtype,*fieldlen;
    int msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    if (msgID != 0) error->one(FLERR,"Bad initial client/server handshake");

  } else {
    int nfield;
    int *fieldID,*fieldtype,*fieldlen;
    int msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    if (msgID != 0) error->one(FLERR,"Bad initial client/server handshake");
    char *pstr = cs->unpack_string(1);
    if (strcmp(pstr,arg[1]) != 0)
      error->one(FLERR,"Mismatch in client/server protocol");

    cs->send(0,0);
  }
}

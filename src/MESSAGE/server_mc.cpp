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

#include "server_mc.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "integrate.h"
#include "input.h"
#include "output.h"
#include "thermo.h"
#include "error.h"

// CSlib interface

#include "cslib.h"

using namespace LAMMPS_NS;
using namespace CSLIB_NS;

enum{NATOMS=1,EINIT,DISPLACE,ACCEPT,RUN};

/* ---------------------------------------------------------------------- */

ServerMC::ServerMC(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ServerMC::loop()
{
  int m;
  double xold[3];
  tagint atomid;

  CSlib *cs = (CSlib *) lmp->cslib;

  if (domain->box_exist == 0)
    error->all(FLERR,"Server command before simulation box is defined");

  if (!atom->map_style) error->all(FLERR,"Server mc requires atom map");
  if (atom->tag_enable == 0) error->all(FLERR,"Server mc requires atom IDs");
  if (sizeof(tagint) != 4) error->all(FLERR,"Server mc requires 32-bit atom IDs");

  // initialize LAMMPS for dynamics

  input->one("run 0");

  // loop on messages
  // receive a message, process it, send return message if necessary

  int msgID,nfield;
  int *fieldID,*fieldtype,*fieldlen;

  while (1) {
    msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    if (msgID < 0) break;

    if (msgID == NATOMS) {

      cs->send(msgID,1);
      cs->pack_int(1,atom->natoms);

    } else if (msgID == EINIT) {

      double dval;
      output->thermo->evaluate_keyword((char *) "pe",&dval);

      cs->send(msgID,2);
      cs->pack_double(1,dval);
      double *coords = NULL;
      if (atom->nlocal) coords = &atom->x[0][0];
      cs->pack_parallel(2,4,atom->nlocal,atom->tag,3,coords);

    } else if (msgID == DISPLACE) {

      atomid = cs->unpack_int(1);
      double *xnew = (double *) cs->unpack(2);
      double **x = atom->x;

      m = atom->map(atomid);
      if (m >= 0 && m < atom->nlocal) {
        xold[0] = x[m][0];
        xold[1] = x[m][1];
        xold[2] = x[m][2];
        x[m][0] = xnew[0];
        x[m][1] = xnew[1];
        x[m][2] = xnew[2];
      }

      input->one("run 0");
      double dval;
      output->thermo->evaluate_keyword((char *) "pe",&dval);

      cs->send(msgID,1);
      cs->pack_double(1,dval);

    } else if (msgID == ACCEPT) {

      int accept = cs->unpack_int(1);
      double **x = atom->x;

      if (!accept) {
        m = atom->map(atomid);
        if (m >= 0 && m < atom->nlocal) {
          x[m][0] = xold[0];
          x[m][1] = xold[1];
          x[m][2] = xold[2];
        }
      }

      cs->send(msgID,0);

    } else if (msgID == RUN) {

      int nsteps = cs->unpack_int(1);

      update->nsteps = nsteps;
      update->firststep = update->ntimestep;
      update->laststep = update->ntimestep + nsteps;

      update->integrate->setup(1);
      update->integrate->run(nsteps);

      cs->send(msgID,0);

    } else error->all(FLERR,"Server received unrecognized message");
  }

  // final reply to client

  cs->send(0,0);

  // clean up

  delete cs;
  lmp->cslib = NULL;
}

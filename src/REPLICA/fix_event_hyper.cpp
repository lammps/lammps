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

#include "fix_event_hyper.h"

#include "comm.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEventHyper::FixEventHyper(LAMMPS *lmp, int narg, char **arg) :
  FixEvent(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix event command");

  restart_global = 1;

  event_number = 0;
  event_timestep = update->ntimestep;
  clock = 0;
}

/* ----------------------------------------------------------------------
   save current atom coords as an event (via call to base class)
   called when an event occurs in some replica
   set event_timestep = when event occurred in a particular replica
   update clock = elapsed time since last event, across all replicas
------------------------------------------------------------------------- */

void FixEventHyper::store_event_hyper(bigint ntimestep, int delta_clock)
{
  store_event();
  event_timestep = ntimestep;
  clock += delta_clock;
  event_number++;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixEventHyper::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = ubuf(event_number).d;
  list[n++] = ubuf(event_timestep).d;
  list[n++] = ubuf(clock).d;
  list[n++] = ubuf(replica_number).d;
  list[n++] = ubuf(correlated_event).d;
  list[n++] = ubuf(ncoincident).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixEventHyper::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;

  event_number = (int) ubuf(list[n++]).i;
  event_timestep = (bigint) ubuf(list[n++]).i;
  clock = (bigint) ubuf(list[n++]).i;
  replica_number = (int) ubuf(list[n++]).i;
  correlated_event = (int) ubuf(list[n++]).i;
  ncoincident = (int) ubuf(list[n++]).i;
}

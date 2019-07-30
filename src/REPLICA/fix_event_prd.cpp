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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "fix_event_prd.h"
#include "comm.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEventPRD::FixEventPRD(LAMMPS *lmp, int narg, char **arg) :
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

void FixEventPRD::store_event_prd(bigint ntimestep, int delta_clock)
{
  store_event();
  event_timestep = ntimestep;
  clock += delta_clock;
  event_number++;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixEventPRD::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = event_number;
  list[n++] = event_timestep;
  list[n++] = clock;
  list[n++] = replica_number;
  list[n++] = correlated_event;
  list[n++] = ncoincident;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixEventPRD::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  event_number = static_cast<int> (list[n++]);
  event_timestep = static_cast<bigint> (list[n++]);
  clock = static_cast<bigint> (list[n++]);
  replica_number = static_cast<int> (list[n++]);
  correlated_event = static_cast<int> (list[n++]);
  ncoincident = static_cast<int> (list[n++]);
}

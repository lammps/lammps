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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "fix_event_tad.h"
#include "comm.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEventTAD::FixEventTAD(LAMMPS *lmp, int narg, char **arg) :
  FixEvent(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix event command");

  restart_global = 1;

  event_number = 0;
  event_timestep = update->ntimestep;
  tlo = 0.0;
  ebarrier = 0.0;
}

/* ----------------------------------------------------------------------
   save current atom coords as an event (via call to base class)
   called when an event occurs in some replica
   set event_timestep = when event occurred
------------------------------------------------------------------------- */

void FixEventTAD::store_event_tad(bigint ntimestep)
{
  store_event();
  event_timestep = ntimestep;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixEventTAD::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = event_number;
  list[n++] = event_timestep;
  list[n++] = tlo;
  list[n++] = ebarrier;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixEventTAD::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;

  event_number = static_cast<int> (list[n++]);
  event_timestep = static_cast<int> (list[n++]);
  tlo = list[n++];
  ebarrier = list[n++];
}

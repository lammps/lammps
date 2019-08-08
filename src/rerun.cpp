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

#include "rerun.h"
#include <cstring>
#include "read_dump.h"
#include "domain.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "output.h"
#include "finish.h"
#include "timer.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Rerun::Rerun(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Rerun::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Rerun command before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal rerun command");

  // list of dump files = args until a keyword

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) break;
    if (strcmp(arg[iarg],"last") == 0) break;
    if (strcmp(arg[iarg],"every") == 0) break;
    if (strcmp(arg[iarg],"skip") == 0) break;
    if (strcmp(arg[iarg],"start") == 0) break;
    if (strcmp(arg[iarg],"stop") == 0) break;
    if (strcmp(arg[iarg],"dump") == 0) break;
    iarg++;
  }
  int nfile = iarg;
  if (nfile == 0 || nfile == narg) error->all(FLERR,"Illegal rerun command");

  // parse optional args up until "dump"
  // use MAXBIGINT -1 so Output can add 1 to it and still be a big int

  bigint first = 0;
  bigint last = MAXBIGINT - 1;
  int nevery = 0;
  int nskip = 1;
  int startflag = 0;
  int stopflag = 0;
  bigint start = -1;
  bigint stop = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      first = force->bnumeric(FLERR,arg[iarg+1]);
      if (first < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"last") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      last = force->bnumeric(FLERR,arg[iarg+1]);
      if (last < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nevery = force->inumeric(FLERR,arg[iarg+1]);
      if (nevery < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"skip") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nskip = force->inumeric(FLERR,arg[iarg+1]);
      if (nskip <= 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      startflag = 1;
      start = force->bnumeric(FLERR,arg[iarg+1]);
      if (start < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"stop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      stopflag = 1;
      stop = force->bnumeric(FLERR,arg[iarg+1]);
      if (stop < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dump") == 0) {
      break;
    } else error->all(FLERR,"Illegal rerun command");
  }

  int nremain = narg - iarg - 1;
  if (nremain <= 0) error->all(FLERR,"Illegal rerun command");
  if (first > last) error->all(FLERR,"Illegal rerun command");
  if (startflag && stopflag && start > stop)
    error->all(FLERR,"Illegal rerun command");

  // pass list of filenames to ReadDump
  // along with post-"dump" args and post-"format" args

  ReadDump *rd = new ReadDump(lmp);

  rd->store_files(nfile,arg);
  if (nremain)
    nremain = rd->fields_and_keywords(nremain,&arg[narg-nremain]);
  else nremain = rd->fields_and_keywords(0,NULL);
  if (nremain) rd->setup_reader(nremain,&arg[narg-nremain]);
  else rd->setup_reader(0,NULL);

  // perform the pseudo run
  // invoke lmp->init() only once
  // read all relevant snapshots
  // use setup_minimal() since atoms are already owned by correct procs
  // addstep_compute_all() insures energy/virial computed on every snapshot

  update->whichflag = 1;

  if (startflag) update->beginstep = update->firststep = start;
  else update->beginstep = update->firststep = first;
  if (stopflag) update->endstep = update->laststep = stop;
  else update->endstep = update->laststep = last;

  int firstflag = 1;
  int ndump = 0;

  lmp->init();

  timer->init();
  timer->barrier_start();

  bigint ntimestep = rd->seek(first,0);
  if (ntimestep < 0)
    error->all(FLERR,"Rerun dump file does not contain requested snapshot");

  while (1) {
    ndump++;
    rd->header(firstflag);
    update->reset_timestep(ntimestep);
    rd->atoms();

    modify->init();
    update->integrate->setup_minimal(1);
    modify->end_of_step();
    if (firstflag) output->setup();
    else if (output->next) output->write(ntimestep);

    firstflag = 0;
    ntimestep = rd->next(ntimestep,last,nevery,nskip);
    if (stopflag && ntimestep > stop)
      error->all(FLERR,"Read rerun dump file timestep > specified stop");
    if (ntimestep < 0) break;
  }

  // insure thermo output on last dump timestep

  output->next_thermo = update->ntimestep;
  output->write(update->ntimestep);

  timer->barrier_stop();

  update->integrate->cleanup();

  // set update->nsteps to ndump for Finish stats to print

  update->nsteps = ndump;

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;

  // clean-up

  delete rd;
}

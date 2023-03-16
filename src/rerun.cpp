// clang-format off
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

#include "rerun.h"

#include "domain.h"
#include "error.h"
#include "finish.h"
#include "input.h"
#include "integrate.h"
#include "modify.h"
#include "output.h"
#include "read_dump.h"
#include "timer.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

#define EPSDT 1.0e-6
/* ---------------------------------------------------------------------- */

Rerun::Rerun(LAMMPS *lmp) : Command(lmp) {}

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
    if (strcmp(arg[iarg],"post") == 0) break;
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
  int postflag = 0;
  bigint start = -1;
  bigint stop = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      first = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (first < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"last") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      last = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (last < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"skip") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nskip = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nskip <= 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      startflag = 1;
      start = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (start < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"stop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      stopflag = 1;
      stop = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (stop < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      postflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
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

  auto rd = new ReadDump(lmp);

  rd->store_files(nfile,arg);
  if (nremain)
    nremain = rd->fields_and_keywords(nremain,&arg[narg-nremain]);
  else nremain = rd->fields_and_keywords(0,nullptr);
  if (nremain) rd->setup_reader(nremain,&arg[narg-nremain]);
  else rd->setup_reader(0,nullptr);

  // perform the pseudo run
  // invoke lmp->init() only once
  // read all relevant snapshots
  // use setup_minimal() since atoms are already owned by correct procs
  // addstep_compute_all() ensures energy/virial computed on every snapshot

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

  while (true) {
    ndump++;
    rd->header(firstflag);
    update->reset_timestep(ntimestep, false);
    rd->atoms();

    modify->init();
    update->integrate->setup_minimal(1);
    modify->end_of_step();

    // fix up the "next_dump" settings for dumps if the sequence differs from what is read in

    for (int idump = 0; idump < output->ndump; ++idump) {
      // dumps triggered by timestep
      if (output->mode_dump[idump] == 0) {
        // rerun has advanced the timestep faster than the dump expected so we need to catch it up
        if (output->next_dump[idump] < ntimestep) {
          // equidistant dumps
          if (output->every_dump[idump]) {
            // if current step compatible with dump frequency, adjust next_dump setting for dump
            if (ntimestep % output->every_dump[idump] == 0) output->next_dump[idump] = ntimestep;
          } else {
            // next dump is determined by variable.
            // advance next timestep computation from variable until it is equal or larger
            // than the current dump timestep; trigger dump only if equal.
            bigint savedstep = update->ntimestep;
            update->ntimestep = output->next_dump[idump];
            bigint nextdump;
            do {
              nextdump = (bigint) input->variable->compute_equal(output->ivar_dump[idump]);
              update->ntimestep = nextdump;
            } while (nextdump < ntimestep);
            output->next_dump[idump] = nextdump;
            update->ntimestep = savedstep;
          }
        }
      } else {
        // dumps triggered by time
        double tcurrent = update->atime + (ntimestep - update->atimestep) * update->dt;
        // rerun time has moved beyond expected time for dump so we need to catch it up
        if (output->next_time_dump[idump] < tcurrent) {
          // equidistant dumps in time
          if (output->every_time_dump[idump] > 0.0) {
            // trigger dump if current time is within +/- half a timestep of the every interval
            double every = output->every_time_dump[idump];
            double rest = fabs(tcurrent/every - round(tcurrent/every)) * every/update->dt;
            if (rest < 0.5) {
              output->next_dump[idump] = ntimestep;
              output->next_time_dump[idump] = tcurrent;
            } else {
              double nexttime = (floor(tcurrent/every) + 1.0) * every;
              output->next_dump[idump] = update->ntimestep
                + (bigint) ((nexttime - (update->atime + (update->ntimestep - update->atimestep) *
                                         update->dt) - EPSDT*update->dt) / update->dt) + 1;
              output->next_time_dump[idump] = nexttime;
            }
          } else {
            // next dump time is determined by variable.
            // advance next time computation from variable until is is equal or larger
            // than the current time/timestep
            bigint savedstep = update->ntimestep;
            update->ntimestep = output->next_dump[idump];
            double nexttime = output->next_time_dump[idump];
            bigint nextstep = output->next_dump[idump];
            while (nextstep < ntimestep) {
              nexttime = input->variable->compute_equal(output->ivar_dump[idump]);
              nextstep = update->ntimestep
                + (bigint) ((nexttime - (update->atime + (update->ntimestep - update->atimestep) *
                                         update->dt) - EPSDT*update->dt) / update->dt) + 1;
              update->ntimestep = nextstep;
            };
            if (ntimestep > 0) {
              output->next_time_dump[idump] = nexttime;
              output->next_dump[idump] = nextstep;
            }
            update->ntimestep = savedstep;
          }
        }
      }
    }

    output->next_dump_any = ntimestep;
    if (firstflag) output->setup();
    else if (output->next) output->write(ntimestep);

    firstflag = 0;
    ntimestep = rd->next(ntimestep,last,nevery,nskip);
    if (stopflag && ntimestep > stop)
      error->all(FLERR,"Read rerun dump file timestep {} > specified stop {}", ntimestep, stop);
    if (ntimestep < 0) break;
  }

  // ensure thermo output on last dump timestep

  output->next_thermo = update->ntimestep;
  output->write(update->ntimestep);

  timer->barrier_stop();

  update->integrate->cleanup();

  // set update->nsteps to ndump for Finish stats to print

  update->nsteps = ndump;

  Finish finish(lmp);
  finish.end(postflag);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;

  // clean-up

  delete rd;
}

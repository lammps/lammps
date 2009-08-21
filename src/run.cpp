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

#include "stdlib.h"
#include "string.h"
#include "run.h"
#include "domain.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "output.h"
#include "finish.h"
#include "input.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 2048

/* ---------------------------------------------------------------------- */

Run::Run(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Run::command(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal run command");

  if (domain->box_exist == 0)
    error->all("Run command before simulation box is defined");

  int nsteps = atoi(arg[0]);

  // parse optional args

  int uptoflag = 0;
  int startflag = 0;
  int stopflag = 0;
  int start,stop;
  int preflag = 1;
  int postflag = 1;
  int nevery = 0;
  char *commandstr = NULL;
    
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"upto") == 0) {
      if (iarg+1 > narg) error->all("Illegal run command");
      uptoflag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      startflag = 1;
      start = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"stop") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      stopflag = 1;
      stop = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pre") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) preflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) preflag = 1;
      else error->all("Illegal run command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) postflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) postflag = 1;
      else error->all("Illegal run command");
      iarg += 2;

    // generate commandstr if last arg is not NULL
    // commandstr = concatenation of all remaining args
    // if an arg has spaces, enclose in quotes since input parser removed them

    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+3 > narg) error->all("Illegal run command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all("Illegal run command");
      if (strcmp(arg[iarg+2],"NULL") != 0) {
	commandstr = new char[MAXLINE];
	commandstr[0] = '\0';
	for (int jarg = iarg+2; jarg < narg; jarg++) {
	  if (strchr(arg[jarg],' ')) strcat(commandstr,"\"");
	  strcat(commandstr,arg[jarg]);
	  if (strchr(arg[jarg],' ')) strcat(commandstr,"\"");
	  strcat(commandstr," ");
	}
      }
      iarg = narg;
    } else error->all("Illegal run command");
  }

  // adjust nsteps if upto was specified

  if (uptoflag) nsteps -= update->ntimestep;

  // error check

  if (uptoflag && nsteps < 0)
    error->all("Run command upto value is before current timestep");
  if (startflag && start > update->ntimestep)
    error->all("Run command start value is after start of run");
  if (stopflag && stop < update->ntimestep + nsteps)
    error->all("Run command stop value is before end of run");

  // perform a single run
  // use start/stop to set begin/end step
  // if pre or 1st run, do System init/setup,
  //   else just init timer and setup output
  // if post, do full Finish, else just print time

  update->whichflag = 1;

  if (nevery == 0) {
    update->nsteps = nsteps;
    update->firststep = update->ntimestep;
    update->laststep = update->ntimestep + nsteps;

    if (startflag) update->beginstep = start;
    else update->beginstep = update->firststep;
    if (stopflag) update->endstep = stop;
    else update->endstep = update->laststep;

    if (preflag || update->first_update == 0) {
      lmp->init();
      update->integrate->setup();
    } else {
      timer->init();
      output->setup(0);
    }

    timer->barrier_start(TIME_LOOP);
    update->integrate->run(nsteps);
    timer->barrier_stop(TIME_LOOP);

    update->integrate->cleanup();

    Finish finish(lmp);
    finish.end(postflag);

  // perform multiple runs optionally interleaved with invocation of a command
  // use start/stop to set begin/end step
  // if pre or 1st iteration of multiple runs, do System init/setup,
  //   else just init timer and setup output
  // if post or last iteration, do full Finish, else just print time

  } else {
    int iter = 0;
    int nleft = nsteps;
    while (nleft > 0 || iter == 0) {
      nsteps = MIN(nleft,nevery);

      update->nsteps = nsteps;
      update->firststep = update->ntimestep;
      update->laststep = update->ntimestep + nsteps;

      if (startflag) update->beginstep = start;
      else update->beginstep = update->firststep;
      if (stopflag) update->endstep = stop;
      else update->endstep = update->laststep;

      if (preflag || iter == 0) {
	lmp->init();
	update->integrate->setup();
      } else {
	timer->init();
	output->setup(0);
      }

      timer->barrier_start(TIME_LOOP);
      update->integrate->run(nsteps);
      timer->barrier_stop(TIME_LOOP);

      update->integrate->cleanup();

      Finish finish(lmp);
      if (postflag || nleft <= nsteps) finish.end(1);
      else finish.end(0);

      // command string may invoke computes via command with variable
      // so wrap with clearstep/addstep

      if (commandstr) {
	modify->clearstep_compute();
	char *command = input->one(commandstr);
	modify->addstep_compute(update->ntimestep + nevery);
      }

      nleft -= nsteps;
      iter++;
    }
  }

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
  delete [] commandstr;
}

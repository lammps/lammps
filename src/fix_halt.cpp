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

#include "fix_halt.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "update.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "atom.h"
#include "neighbor.h"
#include "modify.h"
#include "comm.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{BONDMAX,TLIMIT,DISKFREE,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ,XOR};
enum{HARD,SOFT,CONTINUE};
enum{NOMSG,YESMSG};

/* ---------------------------------------------------------------------- */

FixHalt::FixHalt(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idvar(NULL), dlimit_path(NULL)
{
  if (narg < 7) error->all(FLERR,"Illegal fix halt command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix halt command");

  // comparison args

  idvar = NULL;
  int iarg = 4;

  if (strcmp(arg[iarg],"tlimit") == 0) {
    attribute = TLIMIT;
  } else if (strcmp(arg[iarg],"diskfree") == 0) {
    attribute = DISKFREE;
    dlimit_path = new char[2];
    strcpy(dlimit_path,".");
  } else if (strcmp(arg[iarg],"bondmax") == 0) {
    attribute = BONDMAX;
  } else if (strncmp(arg[iarg],"v_",2) == 0) {
    attribute = VARIABLE;
    int n = strlen(arg[iarg]);
    idvar = new char[n];
    strcpy(idvar,&arg[iarg][2]);
    ivar = input->variable->find(idvar);
    if (ivar < 0) error->all(FLERR,"Could not find fix halt variable name");
    if (input->variable->equalstyle(ivar) == 0)
      error->all(FLERR,"Fix halt variable is not equal-style variable");
  } else error->all(FLERR,"Invalid fix halt attribute");

  ++iarg;
  if (strcmp(arg[iarg],"<") == 0) operation = LT;
  else if (strcmp(arg[iarg],"<=") == 0) operation = LE;
  else if (strcmp(arg[iarg],">") == 0) operation = GT;
  else if (strcmp(arg[iarg],">=") == 0) operation = GE;
  else if (strcmp(arg[iarg],"==") == 0) operation = EQ;
  else if (strcmp(arg[iarg],"!=") == 0) operation = NEQ;
  else if (strcmp(arg[iarg],"|^") == 0) operation = XOR;
  else error->all(FLERR,"Invalid fix halt operator");

  ++iarg;
  value = force->numeric(FLERR,arg[iarg]);

  // parse optional args

  eflag = SOFT;
  msgflag = YESMSG;
  ++iarg;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"error") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix halt command");
      if (strcmp(arg[iarg+1],"hard") == 0) eflag = HARD;
      else if (strcmp(arg[iarg+1],"soft") == 0) eflag = SOFT;
      else if (strcmp(arg[iarg+1],"continue") == 0) eflag = CONTINUE;
      else error->all(FLERR,"Illegal fix halt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"message") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix halt command");
      if (strcmp(arg[iarg+1],"no") == 0) msgflag = NOMSG;
      else if (strcmp(arg[iarg+1],"yes") == 0) msgflag = YESMSG;
      else error->all(FLERR,"Illegal fix halt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"path") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix halt command");
      ++iarg;
      int len = strlen(arg[iarg])+1;
      delete[] dlimit_path;
      dlimit_path = new char[len];
      // strip off quotes, if present
      if ( ((arg[iarg][0] == '"') || (arg[iarg][0] == '\''))
           && (arg[iarg][0] == arg[iarg][len-2]) ) {
        strcpy(dlimit_path,&arg[iarg][1]);
        dlimit_path[len-3] = '\0';
      } else strcpy(dlimit_path,arg[iarg]);
      ++iarg;
    } else error->all(FLERR,"Illegal fix halt command");
  }

  // add nfirst to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  if (attribute == VARIABLE) {
    const bigint nfirst = (update->ntimestep/nevery)*nevery + nevery;
    modify->addstep_compute_all(nfirst);
  }
}

/* ---------------------------------------------------------------------- */

FixHalt::~FixHalt()
{
  delete [] idvar;
  delete [] dlimit_path;
}

/* ---------------------------------------------------------------------- */

int FixHalt::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHalt::init()
{
  // set ivar from current variable list

  if (attribute == VARIABLE) {
    ivar = input->variable->find(idvar);
    if (ivar < 0) error->all(FLERR,"Could not find fix halt variable name");
    if (input->variable->equalstyle(ivar) == 0)
      error->all(FLERR,"Fix halt variable is not equal-style variable");
  }

  // settings used by TLIMIT

  nextstep = (update->ntimestep/nevery)*nevery + nevery;
  thisstep = -1;
  tratio = 0.5;

  // check if disk limit is supported

  if (attribute == DISKFREE) {
    if (diskfree() < 0.0)
      error->all(FLERR,"Disk limit not supported by OS or illegal path");
  }
}

/* ---------------------------------------------------------------------- */

void FixHalt::min_post_force(int /* vflag */)
{
  if (update->ntimestep == thisstep) return;
  if ((update->ntimestep % nevery) == 0) end_of_step();
  thisstep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixHalt::end_of_step()
{
  // variable evaluation may invoke computes so wrap with clear/add

  double attvalue;

  if (attribute == TLIMIT) {
    if (update->ntimestep != nextstep) return;
    attvalue = tlimit();
  } else if (attribute == DISKFREE) {
    attvalue = diskfree();
  } else if (attribute == BONDMAX) {
    attvalue = bondmax();
  } else {
    modify->clearstep_compute();
    attvalue = input->variable->compute_equal(ivar);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // check if halt is triggered, else just return

  if (operation == LT) {
    if (attvalue >= value) return;
  } else if (operation == LE) {
    if (attvalue > value) return;
  } else if (operation == GT) {
    if (attvalue <= value) return;
  } else if (operation == GE) {
    if (attvalue < value) return;
  } else if (operation == EQ) {
    if (attvalue != value) return;
  } else if (operation == NEQ) {
    if (attvalue == value) return;
  } else if (operation == XOR) {
    if ((attvalue == 0.0 && value == 0.0) ||
        (attvalue != 0.0 && value != 0.0)) return;
  }

  // hard halt -> exit LAMMPS
  // soft/continue halt -> trigger timer to break from run loop
  // print message with ID of fix halt in case multiple instances

  char str[128];
  sprintf(str,"Fix halt condition for fix-id %s met on step "
          BIGINT_FORMAT " with value %g",
          id, update->ntimestep, attvalue);

  if (eflag == HARD) {
    error->all(FLERR,str);
  } else if (eflag == SOFT || eflag == CONTINUE) {
    if (comm->me == 0 && msgflag == YESMSG) error->message(FLERR,str);
    timer->force_timeout();
  }
}

/* ----------------------------------------------------------------------
   reset expired timer setting to original value, if requested
------------------------------------------------------------------------- */

void FixHalt::post_run()
{
  // continue halt -> subsequent runs are allowed

  if (eflag == CONTINUE) timer->reset_timeout();
}

/* ----------------------------------------------------------------------
   compute max length of any bond using Neighbor bondlist for each proc
------------------------------------------------------------------------- */

double FixHalt::bondmax()
{
  double **x = atom->x;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int i1,i2;
  double delx,dely,delz,rsq;
  double maxone = 0.0;

  for (int n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    maxone = MAX(rsq,maxone);
  }

  double maxall;
  MPI_Allreduce(&maxone,&maxall,1,MPI_DOUBLE,MPI_MAX,world);

  return sqrt(maxall);
}

/* ----------------------------------------------------------------------
   compute synced elapsed time
   reset nextstep = estimate of timestep when run will end
   first project to 1/2 the run time, thereafter to end of run
------------------------------------------------------------------------- */

double FixHalt::tlimit()
{
  double cpu = timer->elapsed(Timer::TOTAL);
  MPI_Bcast(&cpu,1,MPI_DOUBLE,0,world);

  if (cpu < value) {
    bigint elapsed = update->ntimestep - update->firststep;
    bigint final = update->firststep +
      static_cast<bigint> (tratio*value/cpu * elapsed);
    nextstep = (final/nevery)*nevery + nevery;
    if (nextstep == update->ntimestep) nextstep += nevery;
    tratio = 1.0;
  }

  return cpu;
}

/* ----------------------------------------------------------------------
   determine available disk space, if supported. Return -1 if not.
------------------------------------------------------------------------- */
#if defined(__linux) || defined(__APPLE__)
#include <sys/statvfs.h>
#endif
double FixHalt::diskfree()
{
#if defined(__linux) || defined(__APPLE__)
  struct statvfs fs;
  double disk_free = -1.0;

  if (dlimit_path) {
    disk_free = 1.0e100;
    int rv = statvfs(dlimit_path,&fs);
    if (rv == 0) {
#if defined(__linux)
      disk_free = fs.f_bavail*fs.f_bsize/1048576.0;
#elif defined(__APPLE__)
      disk_free = fs.f_bavail*fs.f_frsize/1048576.0;
#endif
    } else
      disk_free = -1.0;

    MPI_Bcast(&disk_free,1,MPI_DOUBLE,0,world);
  }
  return disk_free;
#else
  return -1.0;
#endif
}

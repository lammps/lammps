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

#include <mpi.h>
#include <cmath>
#include <cstring>
#include "fix_halt.h"
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

enum{BONDMAX,TLIMIT,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ,XOR};
enum{HARD,SOFT,CONTINUE};
enum{NOMSG,YESMSG};

/* ---------------------------------------------------------------------- */

FixHalt::FixHalt(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idvar(NULL)
{
  if (narg < 7) error->all(FLERR,"Illegal fix halt command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix halt command");

  // comparison args

  idvar = NULL;

  if (strcmp(arg[4],"tlimit") == 0) attribute = TLIMIT;
  else if (strcmp(arg[4],"bondmax") == 0) attribute = BONDMAX;
  else if (strncmp(arg[4],"v_",2) == 0) {
    attribute = VARIABLE;
    int n = strlen(arg[4]);
    idvar = new char[n];
    strcpy(idvar,&arg[4][2]);
    ivar = input->variable->find(idvar);
    if (ivar < 0) error->all(FLERR,"Could not find fix halt variable name");
    if (input->variable->equalstyle(ivar) == 0)
      error->all(FLERR,"Fix halt variable is not equal-style variable");
  } else error->all(FLERR,"Invalid fix halt attribute");

  if (strcmp(arg[5],"<") == 0) operation = LT;
  else if (strcmp(arg[5],"<=") == 0) operation = LE;
  else if (strcmp(arg[5],">") == 0) operation = GT;
  else if (strcmp(arg[5],">=") == 0) operation = GE;
  else if (strcmp(arg[5],"==") == 0) operation = EQ;
  else if (strcmp(arg[5],"!=") == 0) operation = NEQ;
  else if (strcmp(arg[5],"|^") == 0) operation = XOR;
  else error->all(FLERR,"Invalid fix halt operator");

  value = force->numeric(FLERR,arg[6]);

  // parse optional args

  eflag = SOFT;
  msgflag = YESMSG;

  int iarg = 7;
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
  sprintf(str,"Fix halt %s condition met on step " BIGINT_FORMAT " with value %g",
          id,update->ntimestep,attvalue);

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

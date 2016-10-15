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

#include <stdlib.h>
#include <string.h>
#include "fix_halt.h"
#include "update.h"
#include "force.h"
#include "input.h"
#include "modify.h"
#include "comm.h"
#include "timer.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{BONDMAX,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ,XOR};
enum{SOFT,HARD};

/* ---------------------------------------------------------------------- */

FixHalt::FixHalt(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idvar(NULL)
{
  if (narg < 7) error->all(FLERR,"Illegal fix halt command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix halt command");

  // comparison args

  idvar = NULL;

  if (strcmp(arg[4],"bondmax") == 0) attribute = BONDMAX;
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

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"error") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix halt command");
      if (strcmp(arg[iarg+1],"soft") == 0) eflag = SOFT;
      else if (strcmp(arg[iarg+1],"hard") == 0) eflag = HARD;
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
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHalt::init()
{
  // set ivar from current variable list

  ivar = input->variable->find(idvar);
  if (ivar < 0) error->all(FLERR,"Could not find fix halt variable name");
  if (input->variable->equalstyle(ivar) == 0)
    error->all(FLERR,"Fix halt variable is not equal-style variable");
}

/* ---------------------------------------------------------------------- */

void FixHalt::end_of_step()
{
  // make a copy of string to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

  double attvalue;

  if (attribute == BONDMAX) attvalue = bondmax();
  else {
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
  // also print ID of fix halt in case multiple ?

  if (eflag == HARD) error->all(FLERR,"Fix halt condition triggered");
  else { // eflag == SOFT
    if (comm->me == 0) error->message(FLERR,"Fix halt condition triggered");
    timer->force_timeout();
  }
}

/* ---------------------------------------------------------------------- */

double FixHalt::bondmax()
{
  // could add custom keyword methods
  // e.g. one to compute max bond length, with no sqrt()

  return 0.0;
}

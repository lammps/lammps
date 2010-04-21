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
   Contributing author: Christina Payne (Vanderbilt U)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_efield.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixEfield::FixEfield(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix efield command");

  efactor = force->qe2f;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else ex = efactor * atof(arg[3]);

  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[4][2]);
  } else ey = efactor * atof(arg[4]);

  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[5][2]);
  } else ez = efactor * atof(arg[5]);

  maxatom = 0;
  efield = NULL;
}

/* ---------------------------------------------------------------------- */

FixEfield::~FixEfield()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  memory->destroy_2d_double_array(efield);
}

/* ---------------------------------------------------------------------- */

int FixEfield::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEfield::init()
{
  // require an atom style with charge defined

  if (atom->q == NULL)
    error->all("Must use charged atom style with fix efield");

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all("Variable name for fix efield does not exist");
    if (input->variable->equalstyle(xvar)) xvarstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xvarstyle = ATOM;
    else error->all("Variable for fix efield is invalid style");
  } else xvarstyle = NONE;
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all("Variable name for fix efield does not exist");
    if (input->variable->equalstyle(yvar)) yvarstyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) yvarstyle = ATOM;
    else error->all("Variable for fix efield is invalid style");
  } else yvarstyle = NONE;
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all("Variable name for fix efield does not exist");
    if (input->variable->equalstyle(zvar)) zvarstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zvarstyle = ATOM;
    else error->all("Variable for fix efield is invalid style");
  } else zvarstyle = NONE;

  if (xvarstyle == ATOM || yvarstyle == ATOM || zvarstyle == ATOM)
    varflag = ATOM;
  else if (xvarstyle == EQUAL || yvarstyle == EQUAL || zvarstyle == EQUAL)
    varflag = EQUAL;
  else varflag = NONE;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixEfield::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ----------------------------------------------------------------------
   apply F = qE
------------------------------------------------------------------------- */

void FixEfield::post_force(int vflag)
{
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // reallocate efield array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy_2d_double_array(efield);
    efield = memory->create_2d_double_array(maxatom,3,"efield:efield");
  }

  if (varflag == NONE) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	f[i][0] += q[i]*ex;
	f[i][1] += q[i]*ey;
	f[i][2] += q[i]*ez;
      }

  } else {
    if (xstr) {
      if (xvarstyle == EQUAL)
	ex = efactor * input->variable->compute_equal(xvar);
      else 
	input->variable->compute_atom(xvar,igroup,&efield[0][0],3,0);
    }
    if (ystr) {
      if (yvarstyle == EQUAL)
	ey = efactor * input->variable->compute_equal(yvar);
      else 
	input->variable->compute_atom(yvar,igroup,&efield[0][1],3,0);
    }
    if (zstr) {
      if (zvarstyle == EQUAL)
	ez = efactor * input->variable->compute_equal(zvar);
      else 
	input->variable->compute_atom(zvar,igroup,&efield[0][2],3,0);
    }

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (xvarstyle == ATOM) f[i][0] += q[i]*efield[i][0];
	else f[i][0] += q[i]*ex;
	if (yvarstyle == ATOM) f[i][1] += q[i]*efield[i][1];
	else f[i][1] += q[i]*ey;
	if (zvarstyle == ATOM) f[i][2] += q[i]*efield[i][2];
	else f[i][2] += q[i]*ez;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixEfield::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixEfield::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

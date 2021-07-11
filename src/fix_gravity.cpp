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

#include "fix_gravity.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_const.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{CHUTE,SPHERICAL,VECTOR};
enum{CONSTANT,EQUAL};          // same as FixPour

/* ---------------------------------------------------------------------- */

FixGravity::FixGravity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  mstr(nullptr), vstr(nullptr), pstr(nullptr), tstr(nullptr),
  xstr(nullptr), ystr(nullptr), zstr(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix gravity command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  mstr = vstr = pstr = tstr = xstr = ystr = zstr = nullptr;
  mstyle = vstyle = pstyle = tstyle = xstyle = ystyle = zstyle = CONSTANT;

  if (utils::strmatch(arg[3],"^v_")) {
    mstr = utils::strdup(arg[3]+2);
    mstyle = EQUAL;
  } else {
    magnitude = utils::numeric(FLERR,arg[3],false,lmp);
    mstyle = CONSTANT;
  }

  int iarg = 4;

  if (strcmp(arg[4],"chute") == 0) {
    if (narg < 6) error->all(FLERR,"Illegal fix gravity command");
    style = CHUTE;
    if (utils::strmatch(arg[5],"^v_")) {
      vstr = utils::strdup(arg[5]+2);
      vstyle = EQUAL;
    } else {
      vert = utils::numeric(FLERR,arg[5],false,lmp);
      vstyle = CONSTANT;
    }
    iarg = 6;

  } else if (strcmp(arg[4],"spherical") == 0) {
    if (narg < 7) error->all(FLERR,"Illegal fix gravity command");
    style = SPHERICAL;
    if (utils::strmatch(arg[5],"^v_")) {
      pstr = utils::strdup(arg[5]+2);
      pstyle = EQUAL;
    } else {
      phi = utils::numeric(FLERR,arg[5],false,lmp);
      pstyle = CONSTANT;
    }
    if (utils::strmatch(arg[6],"^v_")) {
      tstr = utils::strdup(arg[6]+2);
      tstyle = EQUAL;
    } else {
      theta = utils::numeric(FLERR,arg[6],false,lmp);
      tstyle = CONSTANT;
    }
    iarg = 7;

  } else if (strcmp(arg[4],"vector") == 0) {
    if (narg < 8) error->all(FLERR,"Illegal fix gravity command");
    style = VECTOR;
    if (utils::strmatch(arg[5],"^v_")) {
      xstr = utils::strdup(arg[5]+2);
      xstyle = EQUAL;
    } else {
      xdir = utils::numeric(FLERR,arg[5],false,lmp);
      xstyle = CONSTANT;
    }
    if (utils::strmatch(arg[6],"^v_")) {
      ystr = utils::strdup(arg[6]+2);
      ystyle = EQUAL;
    } else {
      ydir = utils::numeric(FLERR,arg[6],false,lmp);
      ystyle = CONSTANT;
    }
    if (utils::strmatch(arg[7],"^v_")) {
      zstr = utils::strdup(arg[7]+2);
      zstyle = EQUAL;
    } else {
      zdir = utils::numeric(FLERR,arg[7],false,lmp);
      zstyle = CONSTANT;
    }
    iarg = 8;

  } else error->all(FLERR,"Illegal fix gravity command");

  // optional keywords

  disable = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"disable") == 0) {
      disable = 1;
      iarg++;
    } else error->all(FLERR,"Illegal fix gravity command");
  }

  // initializations

  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;

  eflag = 0;
  egrav = 0.0;

  // set gravity components once and for all if CONSTANT

  varflag = CONSTANT;
  if (mstyle != CONSTANT || vstyle != CONSTANT || pstyle != CONSTANT ||
      tstyle != CONSTANT || xstyle != CONSTANT || ystyle != CONSTANT ||
      zstyle != CONSTANT) varflag = EQUAL;

  if (varflag == CONSTANT) set_acceleration();

}

/* ---------------------------------------------------------------------- */

FixGravity::~FixGravity()
{
  if (copymode) return;

  delete [] mstr;
  delete [] vstr;
  delete [] pstr;
  delete [] tstr;
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
}

/* ---------------------------------------------------------------------- */

int FixGravity::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGravity::init()
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  // check variables

  if (mstr) {
    mvar = input->variable->find(mstr);
    if (mvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(mvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (vstr) {
    vvar = input->variable->find(vstr);
    if (vvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(tvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void FixGravity::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixGravity::post_force(int /*vflag*/)
{
  // update gravity due to variables

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
    if (vstyle == EQUAL) vert = input->variable->compute_equal(vvar);
    if (pstyle == EQUAL) phi = input->variable->compute_equal(pvar);
    if (tstyle == EQUAL) theta = input->variable->compute_equal(tvar);
    if (xstyle == EQUAL) xdir = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) ydir = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zdir = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);

    set_acceleration();
  }

  // just exit if application of force is disabled

  if (disable) return;

  // apply gravity force to each particle

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;

  eflag = 0;
  egrav = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] += massone*xacc;
        f[i][1] += massone*yacc;
        f[i][2] += massone*zacc;
        egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] += massone*xacc;
        f[i][1] += massone*yacc;
        f[i][2] += massone*zacc;
        egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixGravity::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGravity::set_acceleration()
{
  if (style == CHUTE || style == SPHERICAL) {
    if (style == CHUTE) {
      phi = 0.0;
      theta = 180.0 - vert;
    }
    if (domain->dimension == 3) {
      xgrav = sin(degree2rad * theta) * cos(degree2rad * phi);
      ygrav = sin(degree2rad * theta) * sin(degree2rad * phi);
      zgrav = cos(degree2rad * theta);
    } else {
      xgrav = sin(degree2rad * theta);
      ygrav = cos(degree2rad * theta);
      zgrav = 0.0;
    }
  } else if (style == VECTOR) {
    if (domain->dimension == 3) {
      double length = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);
      xgrav = xdir/length;
      ygrav = ydir/length;
      zgrav = zdir/length;
    } else {
      double length = sqrt(xdir*xdir + ydir*ydir);
      xgrav = xdir/length;
      ygrav = ydir/length;
      zgrav = 0.0;
    }
  }

  gvec[0] = xacc = magnitude*xgrav;
  gvec[1] = yacc = magnitude*ygrav;
  gvec[2] = zacc = magnitude*zgrav;
}

/* ----------------------------------------------------------------------
   potential energy in gravity field
------------------------------------------------------------------------- */

double FixGravity::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&egrav,&egrav_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return egrav_all;
}

/* ----------------------------------------------------------------------
   extract current gravity direction vector
------------------------------------------------------------------------- */

void *FixGravity::extract(const char *name, int &dim)
{
  if (strcmp(name,"gvec") == 0) {
    dim = 1;
    return (void *) gvec;
  }
  return nullptr;
}

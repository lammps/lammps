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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_gravity.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{CHUTE,SPHERICAL,GRADIENT,VECTOR};

/* ---------------------------------------------------------------------- */

FixGravity::FixGravity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix gravity command");

  time_depend = 1;

  magnitude = atof(arg[3]);

  if (strcmp(arg[4],"chute") == 0) {
    if (narg != 6) error->all("Illegal fix gravity command");
    style = CHUTE;
    phi = 0.0;
    theta = 180.0 - atof(arg[5]);
  } else if (strcmp(arg[4],"spherical") == 0) {
    if (narg != 7) error->all("Illegal fix gravity command");
    style = SPHERICAL;
    phi = atof(arg[5]);
    theta = atof(arg[6]);
  } else if (strcmp(arg[4],"gradient") == 0) {
    if (narg != 9) error->all("Illegal fix gravity command");
    style = GRADIENT;
    phi = atof(arg[5]);
    theta = atof(arg[6]);
    phigrad = atof(arg[7]);
    thetagrad = atof(arg[8]);
  } else if (strcmp(arg[4],"vector") == 0) {
    if (narg != 8) error->all("Illegal fix gravity command");
    style = VECTOR;
    xdir = atof(arg[5]);
    ydir = atof(arg[6]);
    zdir = atof(arg[7]);
  } else error->all("Illegal fix gravity command");

  double PI = 2.0 * asin(1.0);
  degree2rad = 2.0*PI / 360.0;

  if (style == CHUTE || style == SPHERICAL || style == GRADIENT) {
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

  time_initial = update->ntimestep;
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
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  dt = update->dt;

  xacc = magnitude*xgrav;
  yacc = magnitude*ygrav;
  zacc = magnitude*zgrav;
}

/* ---------------------------------------------------------------------- */

void FixGravity::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixGravity::post_force(int vflag)
{
  // update direction of gravity vector if gradient style

  if (style == GRADIENT) {
    if (domain->dimension == 3) {
      double phi_current = degree2rad * 
	(phi + (update->ntimestep-time_initial)*dt*phigrad*360.0);
      double theta_current = degree2rad * 
	(theta + (update->ntimestep-time_initial)*dt*thetagrad*360.0);
      xgrav = sin(theta_current) * cos(phi_current);
      ygrav = sin(theta_current) * sin(phi_current);
      zgrav = cos(theta_current);
    } else {
      double theta_current = degree2rad * 
	(theta + (update->ntimestep-time_initial)*dt*thetagrad*360.0);
      xgrav = sin(theta_current);
      ygrav = cos(theta_current);
    }
    xacc = magnitude*xgrav;
    yacc = magnitude*ygrav;
    zacc = magnitude*zgrav;
  }

  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = rmass[i];
	f[i][0] += massone*xacc;
	f[i][1] += massone*yacc;
	f[i][2] += massone*zacc;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = mass[type[i]];
	f[i][0] += massone*xacc;
	f[i][1] += massone*yacc;
	f[i][2] += massone*zacc;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixGravity::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

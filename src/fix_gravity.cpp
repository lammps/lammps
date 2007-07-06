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
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixGravity::FixGravity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix gravity command");

  if (strcmp(arg[3],"chute") == 0) {
    if (narg != 5) error->all("Illegal fix gravity command");
    dynamic = 0;
    granular = 1;
    phi = 0.0;
    theta = 180.0 - atof(arg[4]);
  } else if (strcmp(arg[3],"spherical") == 0) {
    if (narg != 6) error->all("Illegal fix gravity command");
    dynamic = 0;
    granular = 1;
    phi = atof(arg[4]);
    theta = atof(arg[5]);
  } else if (strcmp(arg[3],"gradient") == 0) {
    if (narg != 8) error->all("Illegal fix gravity command");
    dynamic = 1;
    granular = 1;
    phi = atof(arg[4]);
    theta = atof(arg[5]);
    phigrad = atof(arg[6]);
    thetagrad = atof(arg[7]);
  } else if (strcmp(arg[3],"vector") == 0) {
    if (narg != 8) error->all("Illegal fix gravity command");
    dynamic = 0;
    granular = 0;
    magnitude = atof(arg[4]);
    xdir = atof(arg[5]);
    ydir = atof(arg[6]);
    zdir = atof(arg[7]);
  } else error->all("Illegal fix gravity command");

  double PI = 2.0 * asin(1.0);
  degree2rad = 2.0*PI / 360.0;
  time_initial = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixGravity::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGravity::init()
{
  dt = update->dt;

  if (granular) {
    if (domain->dimension == 3) {
      xgrav = sin(degree2rad * theta) * cos(degree2rad * phi);
      ygrav = sin(degree2rad * theta) * sin(degree2rad * phi);
      zgrav = cos(degree2rad * theta);
    } else {
      xgrav = sin(degree2rad * theta);
      ygrav = cos(degree2rad * theta);
      zgrav = 0.0;
    }
  } else {
    if (domain->dimension == 3) {
      double length = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);
      xgrav = magnitude * xdir/length;
      ygrav = magnitude * ydir/length;
      zgrav = magnitude * zdir/length;
    } else {
      double length = sqrt(xdir*xdir + ydir*ydir);
      xgrav = magnitude * xdir/length;
      ygrav = magnitude * ydir/length;
      zgrav = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGravity::setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixGravity::post_force(int vflag)
{
  // update direction of gravity vector if dynamic

  if (dynamic) {
    double phi_current = degree2rad * 
      (phi + (update->ntimestep-time_initial)*dt*phigrad*360.0);
    double theta_current = degree2rad * 
      (theta + (update->ntimestep-time_initial)*dt*thetagrad*360.0);
    xgrav = sin(theta_current) * cos(phi_current);
    ygrav = sin(theta_current) * sin(phi_current);
    zgrav = cos(theta_current);
  }

  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = mass[type[i]];
	f[i][0] += massone*xgrav;
	f[i][1] += massone*ygrav;
	f[i][2] += massone*zgrav;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	f[i][0] += rmass[i]*xgrav;
	f[i][1] += rmass[i]*ygrav;
	f[i][2] += rmass[i]*zgrav;
      }
  }
}

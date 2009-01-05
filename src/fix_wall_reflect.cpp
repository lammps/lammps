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

#include "string.h"
#include "fix_wall_reflect.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallReflect::FixWallReflect(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix wall/reflect command");

  xloflag = xhiflag = yloflag = yhiflag = zloflag = zhiflag = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"xlo") == 0) xloflag = 1;
    else if (strcmp(arg[iarg],"xhi") == 0) xhiflag = 1;
    else if (strcmp(arg[iarg],"ylo") == 0) yloflag = 1;
    else if (strcmp(arg[iarg],"yhi") == 0) yhiflag = 1;
    else if (strcmp(arg[iarg],"zlo") == 0) zloflag = 1;
    else if (strcmp(arg[iarg],"zhi") == 0) zhiflag = 1;
    else error->all("Illegal fix wall/reflect command");
  }

  if ((xloflag || xhiflag) && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if ((yloflag || yhiflag) && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if ((zloflag || zhiflag) && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");
}

/* ---------------------------------------------------------------------- */

int FixWallReflect::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::post_integrate()
{
  double xlo = domain->boxlo[0];
  double xhi = domain->boxhi[0];
  double ylo = domain->boxlo[1];
  double yhi = domain->boxhi[1];
  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (xloflag && x[i][0] < xlo) {
	x[i][0] = xlo + (xlo - x[i][0]);
	v[i][0] = -v[i][0];
      }
      if (xhiflag && x[i][0] > xhi) {
	x[i][0] = xhi - (x[i][0] - xhi);
	v[i][0] = -v[i][0];
      }
      if (yloflag && x[i][1] < ylo) {
	x[i][1] = ylo + (ylo - x[i][1]);
	v[i][1] = -v[i][1];
      }
      if (yhiflag && x[i][1] > yhi) {
	x[i][1] = yhi - (x[i][1] - yhi);
	v[i][1] = -v[i][1];
      }
      if (zloflag && x[i][2] < zlo) {
	x[i][2] = zlo + (zlo - x[i][2]);
	v[i][2] = -v[i][2];
      }
      if (zhiflag && x[i][2] > zhi) {
	x[i][2] = zhi - (x[i][2] - zhi);
	v[i][2] = -v[i][2];
      }
    }
  }
}

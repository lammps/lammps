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
#include "stdlib.h"
#include "fix_set_force.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixSetForce::FixSetForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix setforce command");

  flagx = flagy = flagz = 1;
  if (strcmp(arg[3],"NULL") == 0) flagx = 0;
  else xvalue = atof(arg[3]);
  if (strcmp(arg[4],"NULL") == 0) flagy = 0;
  else yvalue = atof(arg[4]);
  if (strcmp(arg[5],"NULL") == 0) flagz = 0;
  else zvalue = atof(arg[5]);
}

/* ---------------------------------------------------------------------- */

int FixSetForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetForce::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixSetForce::setup()
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(1,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixSetForce::min_setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixSetForce::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      foriginal[0] += f[i][0];
      foriginal[1] += f[i][1];
      foriginal[2] += f[i][2];
      if (flagx) f[i][0] = xvalue;
      if (flagy) f[i][1] = yvalue;
      if (flagz) f[i][2] = zvalue;
    }
}

/* ---------------------------------------------------------------------- */

void FixSetForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  // set force to desired value on outermost level, 0.0 on other levels

  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	foriginal[0] += f[i][0];
	foriginal[1] += f[i][1];
	foriginal[2] += f[i][2];
	if (flagx) f[i][0] = 0.0;
	if (flagy) f[i][1] = 0.0;
	if (flagz) f[i][2] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   allow setforce values to be printed with thermo output
   n = 1,2,3 = components of total force on fix group before reset
------------------------------------------------------------------------- */

double FixSetForce::thermo(int n)
{
  if (n >= 1 && n <= 3) {
    double ftotal;
    MPI_Allreduce(&foriginal[n-1],&ftotal,1,MPI_DOUBLE,MPI_SUM,world);
    return ftotal;
  } else return 0.0;
}

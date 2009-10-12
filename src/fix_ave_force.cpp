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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ave_force.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAveForce::FixAveForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all("Illegal fix aveforce command");

  vector_flag = 1;
  size_vector = 3;
  scalar_vector_freq = 1;
  extvector = 1;

  xflag = yflag = zflag = 1;
  if (strcmp(arg[3],"NULL") == 0) xflag = 0;
  else xvalue = atof(arg[3]);
  if (strcmp(arg[4],"NULL") == 0) yflag = 0;
  else yvalue = atof(arg[4]);
  if (strcmp(arg[5],"NULL") == 0) zflag = 0;
  else zvalue = atof(arg[5]);

  // optional args

  iregion = -1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix aveforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all("Fix aveforce region ID does not exist");
      iarg += 2;
    } else error->all("Illegal fix aveforce command");

  }

  foriginal_all[0] = foriginal_all[1] = 
    foriginal_all[2] = foriginal_all[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixAveForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveForce::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAveForce::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixAveForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAveForce::post_force(int vflag)
{
  // sum forces on participating atoms

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double foriginal[4];
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (iregion >= 0 && 
          !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	continue;

      foriginal[0] += f[i][0];
      foriginal[1] += f[i][1];
      foriginal[2] += f[i][2];
      foriginal[3] += 1.0;
    }

  // average the force on participating atoms
  // add in requested amount

  MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);

  int ncount = static_cast<int> (foriginal_all[3]);
  if (ncount == 0) return;

  double fave[3];
  fave[0] = foriginal_all[0]/ncount + xvalue;
  fave[1] = foriginal_all[1]/ncount + yvalue;
  fave[2] = foriginal_all[2]/ncount + zvalue;

  // set force of all participating atoms to same value
  // only for active dimensions

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (iregion >= 0 && 
          !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	continue;

      if (xflag) f[i][0] = fave[0];
      if (yflag) f[i][1] = fave[1];
      if (zflag) f[i][2] = fave[2];
    }
}

/* ---------------------------------------------------------------------- */

void FixAveForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  // ave + extra force on outermost level
  // just ave on inner levels

  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double foriginal[4];
    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (iregion >= 0 && 
	    !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	  continue;

	foriginal[0] += f[i][0];
	foriginal[1] += f[i][1];
	foriginal[2] += f[i][2];
	foriginal[3] += 1.0;
      }

    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);

    int ncount = static_cast<int> (foriginal_all[3]);
    if (ncount == 0) return;

    double fave[3];
    fave[0] = foriginal_all[0]/ncount;
    fave[1] = foriginal_all[1]/ncount;
    fave[2] = foriginal_all[2]/ncount;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (iregion >= 0 && 
	    !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	  continue;

	if (xflag) f[i][0] = fave[0];
	if (yflag) f[i][1] = fave[1];
	if (zflag) f[i][2] = fave[2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixAveForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAveForce::compute_vector(int n)
{
  return foriginal_all[n];
}

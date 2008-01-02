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
#include "stdlib.h"
#include "string.h"
#include "fix_dt_reset.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "domain.h"
#include "lattice.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "fix.h"
#include "output.h"
#include "dump.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixDtReset::FixDtReset(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix dt/reset command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 1;
  scalar_vector_freq = 1;
  extscalar = 0;
  extvector = 0;

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix dt/reset command");

  minbound = maxbound = 1;
  tmin = tmax = 0.0;
  if (strcmp(arg[4],"INF") == 0) minbound = 0;
  else tmin = atof(arg[4]);
  if (strcmp(arg[5],"INF") == 0) maxbound = 0;
  else tmax = atof(arg[5]);
  xmax = atof(arg[6]);

  if (minbound && tmin < 0.0) error->all("Illegal fix dt/reset command");
  if (maxbound && tmax < 0.0) error->all("Illegal fix dt/reset command");
  if (minbound && maxbound && tmin > tmax)
    error->all("Illegal fix dt/reset command");
  if (xmax <= 0.0) error->all("Illegal fix dt/reset command");

  int scaleflag = 1;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix dt/reset command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal fix dt/reset command");
      iarg += 2;
    } else error->all("Illegal fix dt/reset command");
  }

  // setup scaling, based on xlattice parameter

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of fix dt/reset with undefined lattice");
  if (scaleflag) xmax *= domain->lattice->xlattice;

  // initializations

  t_elapsed = 0.0;
  laststep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixDtReset::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDtReset::init()
{
  // set rRESPA flag

  respaflag = 0;
  if (strcmp(update->integrate_style,"respa") == 0) respaflag = 1;

  // check for DCD or XTC dumps

  for (int i = 0; i < output->ndump; i++)
    if ((strcmp(output->dump[i]->style,"dcd") == 0 ||
	strcmp(output->dump[i]->style,"xtc") == 0) && comm->me == 0)
      error->warning("Dump dcd/xtc timestamp may be wrong with fix dt/reset");
}

/* ---------------------------------------------------------------------- */

void FixDtReset::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixDtReset::end_of_step()
{
  // accumulate total time based on previous timestep

  t_elapsed += (update->ntimestep - laststep) * update->dt;

  // compute vmax and amax of any atom in group

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vsq,fsq,asq,ms;
  double bound[2],bound_all[2];
  bound[0] = bound[1] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      fsq = f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
      if (mass) ms = mass[type[i]];
      else ms = rmass[i];
      asq = fsq/ms/ms;
      bound[0] = MAX(bound[0],vsq);
      bound[1] = MAX(bound[1],asq);
    }

  MPI_Allreduce(bound,bound_all,2,MPI_DOUBLE,MPI_MAX,world);

  // set new timestep

  double dt,dt1,dt2;

  if (bound_all[0] > 0.0) dt1 = xmax/sqrt(bound_all[0]);
  else dt1 = BIG;
  if (bound_all[1] > 0.0) dt2 = sqrt(2.0 * xmax/sqrt(bound_all[1]));
  else dt2 = BIG;

  dt = MIN(dt1,dt2);
  if (minbound) dt = MAX(dt,tmin);
  if (maxbound) dt = MIN(dt,tmax);
  
  // reset update->dt and other classes that depend on it
  // rRESPA, pair style, fixes

  laststep = update->ntimestep;
  if (dt == update->dt) return;

  update->dt = dt;
  if (respaflag) update->integrate->reset_dt();
  if (force->pair) force->pair->reset_dt();
  for (int i = 0; i < modify->nfix; i++) modify->fix[i]->reset_dt();
}

/* ---------------------------------------------------------------------- */

double FixDtReset::compute_scalar()
{
  return update->dt;
}

/* ---------------------------------------------------------------------- */

double FixDtReset::compute_vector(int n)
{
  return t_elapsed;
}

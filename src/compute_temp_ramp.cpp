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
#include <cstdlib>
#include <cstring>
#include "compute_temp_ramp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "fix.h"
#include "domain.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempRamp::ComputeTempRamp(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  // parse optional args

  scaleflag = 1;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/ramp command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal compute temp/ramp command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/ramp command");
  }

  // setup scaling

  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // read standard args and apply scaling

  if (strcmp(arg[3],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[3],"vy") == 0) v_dim = 1;
  else if (strcmp(arg[3],"vz") == 0) v_dim = 2;
  else error->all(FLERR,"Illegal compute temp/ramp command");

  if (v_dim == 0) {
    v_lo = xscale*force->numeric(FLERR,arg[4]);
    v_hi = xscale*force->numeric(FLERR,arg[5]);
  } else if (v_dim == 1) {
    v_lo = yscale*force->numeric(FLERR,arg[4]);
    v_hi = yscale*force->numeric(FLERR,arg[5]);
  } else if (v_dim == 2) {
    v_lo = zscale*force->numeric(FLERR,arg[4]);
    v_hi = zscale*force->numeric(FLERR,arg[5]);
  }

  if (strcmp(arg[6],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[6],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[6],"z") == 0) coord_dim = 2;
  else error->all(FLERR,"Illegal compute temp/ramp command");

  if (coord_dim == 0) {
    coord_lo = xscale*force->numeric(FLERR,arg[7]);
    coord_hi = xscale*force->numeric(FLERR,arg[8]);
  } else if (coord_dim == 1) {
    coord_lo = yscale*force->numeric(FLERR,arg[7]);
    coord_hi = yscale*force->numeric(FLERR,arg[8]);
  } else if (coord_dim == 2) {
    coord_lo = zscale*force->numeric(FLERR,arg[7]);
    coord_hi = zscale*force->numeric(FLERR,arg[8]);
  }

  maxbias = 0;
  vbiasall = NULL;
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempRamp::~ComputeTempRamp()
{
  memory->destroy(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRamp::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempRamp::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempRamp::compute_scalar()
{
  double fraction,vramp,vthermal[3];

  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      vthermal[0] = v[i][0];
      vthermal[1] = v[i][1];
      vthermal[2] = v[i][2];
      vthermal[v_dim] -= vramp;
      if (rmass)
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2]) * rmass[i];
      else
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2]) * mass[type[i]];
    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRamp::compute_vector()
{
  int i;
  double fraction,vramp,vthermal[3];

  invoked_vector = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      vthermal[0] = v[i][0];
      vthermal[1] = v[i][1];
      vthermal[2] = v[i][2];
      vthermal[v_dim] -= vramp;

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * vthermal[0]*vthermal[0];
      t[1] += massone * vthermal[1]*vthermal[1];
      t[2] += massone * vthermal[2]*vthermal[2];
      t[3] += massone * vthermal[0]*vthermal[1];
      t[4] += massone * vthermal[0]*vthermal[2];
      t[5] += massone * vthermal[1]*vthermal[2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRamp::remove_bias(int i, double *v)
{
  double fraction = (atom->x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
  fraction = MAX(fraction,0.0);
  fraction = MIN(fraction,1.0);
  vbias[v_dim] = v_lo + fraction*(v_hi - v_lo);
  v[v_dim] -= vbias[v_dim];
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRamp::remove_bias_thr(int i, double *v, double *b)
{
  double fraction = (atom->x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
  fraction = MAX(fraction,0.0);
  fraction = MIN(fraction,1.0);
  b[v_dim] = v_lo + fraction*(v_hi - v_lo);
  v[v_dim] -= b[v_dim];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRamp::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/ramp:vbiasall");
  }

  double fraction;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (atom->x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vbiasall[i][v_dim] = v_lo + fraction*(v_hi - v_lo);
      v[i][v_dim] -= vbiasall[i][v_dim];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempRamp::restore_bias(int /*i*/, double *v)
{
  v[v_dim] += vbias[v_dim];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeTempRamp::restore_bias_thr(int /*i*/, double *v, double *b)
{
  v[v_dim] += b[v_dim];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempRamp::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      v[i][v_dim] += vbiasall[i][v_dim];
}

/* ---------------------------------------------------------------------- */

double ComputeTempRamp::memory_usage()
{
  double bytes = 3*maxbias * sizeof(double);
  return bytes;
}

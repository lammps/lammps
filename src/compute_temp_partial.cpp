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
#include "compute_temp_partial.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempPartial::ComputeTempPartial(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal compute temp/partial command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  xflag = force->inumeric(FLERR,arg[3]);
  yflag = force->inumeric(FLERR,arg[4]);
  zflag = force->inumeric(FLERR,arg[5]);
  if ((xflag != 0 && xflag != 1) || (yflag != 0 && yflag != 1)
      || (zflag != 0 && zflag != 1))
    error->all(FLERR,"Illegal compute temp/partial command");
  if (zflag && domain->dimension == 2)
    error->all(FLERR,"Compute temp/partial cannot use vz for 2d systemx");

  maxbias = 0;
  vbiasall = NULL;
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempPartial::~ComputeTempPartial()
{
  memory->destroy(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartial::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ----------------------------------------------------------------------
   DOF for a body of N atoms with S constraints (e.g. from SHAKE)
   DOF = nper/dim (dim*N - S), where dim = dimensionality = 2 or 3
------------------------------------------------------------------------- */

void ComputeTempPartial::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  int nper = xflag+yflag+zflag;
  dof = nper * natoms_temp;

  // distribute extra dofs evenly across active dimensions

  dof -= (1.0*nper/domain->dimension)*(fix_dof + extra_dof);
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputeTempPartial::dof_remove(int /*i*/)
{
  int nper = xflag+yflag+zflag;
  return (domain->dimension - nper);
}

/* ---------------------------------------------------------------------- */

double ComputeTempPartial::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] +
              zflag*v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] +
              zflag*v[i][2]*v[i][2]) * mass[type[i]];
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartial::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

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
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * xflag*v[i][0]*v[i][0];
      t[1] += massone * yflag*v[i][1]*v[i][1];
      t[2] += massone * zflag*v[i][2]*v[i][2];
      t[3] += massone * xflag*yflag*v[i][0]*v[i][1];
      t[4] += massone * xflag*zflag*v[i][0]*v[i][2];
      t[5] += massone * yflag*zflag*v[i][1]*v[i][2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempPartial::remove_bias(int /*i*/, double *v)
{
  if (!xflag) {
    vbias[0] = v[0];
    v[0] = 0.0;
  }
  if (!yflag) {
    vbias[1] = v[1];
    v[1] = 0.0;
  }
  if (!zflag) {
    vbias[2] = v[2];
    v[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempPartial::remove_bias_thr(int /*i*/, double *v, double *b)
{
  if (!xflag) {
    b[0] = v[0];
    v[0] = 0.0;
  }
  if (!yflag) {
    b[1] = v[1];
    v[1] = 0.0;
  }
  if (!zflag) {
    b[2] = v[2];
    v[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempPartial::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/partial:vbiasall");
  }

  if (!xflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        vbiasall[i][0] = v[i][0];
        v[i][0] = 0.0;
      }
  }
  if (!yflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        vbiasall[i][1] = v[i][1];
        v[i][1] = 0.0;
      }
  }
  if (!zflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        vbiasall[i][2] = v[i][2];
        v[i][2] = 0.0;
      }
  }
}

/* ----------------------------------------------------------------------
   reset thermal velocity of all atoms to be consistent with bias
   called from velocity command after it creates thermal velocities
   this re-zero components that should stay zero
------------------------------------------------------------------------- */

void ComputeTempPartial::reapply_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (!xflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) v[i][0] = 0.0;
  }
  if (!yflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) v[i][1] = 0.0;
  }
  if (!zflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) v[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempPartial::restore_bias(int /*i*/, double *v)
{
  if (!xflag) v[0] += vbias[0];
  if (!yflag) v[1] += vbias[1];
  if (!zflag) v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeTempPartial::restore_bias_thr(int /*i*/, double *v, double *b)
{
  if (!xflag) v[0] += b[0];
  if (!yflag) v[1] += b[1];
  if (!zflag) v[2] += b[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempPartial::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (!xflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        v[i][0] += vbiasall[i][0];
  }
  if (!yflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        v[i][1] += vbiasall[i][1];
  }
  if (!zflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        v[i][2] += vbiasall[i][2];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeTempPartial::memory_usage()
{
  double bytes = 3*maxbias * sizeof(double);
  return bytes;
}

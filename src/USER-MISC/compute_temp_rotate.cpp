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
   Contributing author: Laurent Joly (U Lyon, France), ljoly.ulyon@gmail.com
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "compute_temp_rotate.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "domain.h"
#include "lattice.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempRotate::ComputeTempRotate(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp/rotate command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  maxbias = 0;
  vbiasall = NULL;
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempRotate::~ComputeTempRotate()
{
  memory->destroy(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRotate::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeTempRotate::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempRotate::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempRotate::compute_scalar()
{
  double vthermal[3];
  double vcm[3],xcm[3],inertia[3][3],angmom[3],omega[3];
  double dx,dy,dz;
  double unwrap[3];

  invoked_scalar = update->ntimestep;

  if (dynamic) masstotal = group->mass(igroup);
  group->vcm(igroup,masstotal,vcm);
  group->xcm(igroup,masstotal,xcm);
  group->inertia(igroup,xcm,inertia);
  group->angmom(igroup,xcm,angmom);
  group->omega(angmom,inertia,omega);

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/rotate:vbiasall");
  }

  double t = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      vbiasall[i][0] = vcm[0] + dz*omega[1]-dy*omega[2];
      vbiasall[i][1] = vcm[1] + dx*omega[2]-dz*omega[0];
      vbiasall[i][2] = vcm[2] + dy*omega[0]-dx*omega[1];
      vthermal[0] = v[i][0] - vbiasall[i][0];
      vthermal[1] = v[i][1] - vbiasall[i][1];
      vthermal[2] = v[i][2] - vbiasall[i][2];
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

void ComputeTempRotate::compute_vector()
{
  double vthermal[3];
  double vcm[3],xcm[3],inertia[3][3],angmom[3],omega[3];
  double dx,dy,dz;
  double unwrap[3];

  invoked_vector = update->ntimestep;

  if (dynamic) masstotal = group->mass(igroup);
  group->vcm(igroup,masstotal,vcm);
  group->xcm(igroup,masstotal,xcm);
  group->inertia(igroup,xcm,inertia);
  group->angmom(igroup,xcm,angmom);
  group->omega(angmom,inertia,omega);

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/rotate:vbiasall");
  }

  double massone,t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      vbiasall[i][0] = vcm[0] + dz*omega[1]-dy*omega[2];
      vbiasall[i][1] = vcm[1] + dx*omega[2]-dz*omega[0];
      vbiasall[i][2] = vcm[2] + dy*omega[0]-dx*omega[1];
      vthermal[0] = v[i][0] - vbiasall[i][0];
      vthermal[1] = v[i][1] - vbiasall[i][1];
      vthermal[2] = v[i][2] - vbiasall[i][2];
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
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRotate::remove_bias(int i, double *v)
{
  v[0] -= vbiasall[i][0];
  v[1] -= vbiasall[i][1];
  v[2] -= vbiasall[i][2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRotate::remove_bias_thr(int i, double *v, double *)
{
  v[0] -= vbiasall[i][0];
  v[1] -= vbiasall[i][1];
  v[2] -= vbiasall[i][2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRotate::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vbiasall[i][0];
      v[i][1] -= vbiasall[i][1];
      v[i][2] -= vbiasall[i][2];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempRotate::restore_bias(int i, double *v)
{
  v[0] += vbiasall[i][0];
  v[1] += vbiasall[i][1];
  v[2] += vbiasall[i][2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called
------------------------------------------------------------------------- */

void ComputeTempRotate::restore_bias_thr(int i, double *v, double *)
{
  v[0] += vbiasall[i][0];
  v[1] += vbiasall[i][1];
  v[2] += vbiasall[i][2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempRotate::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] += vbiasall[i][0];
      v[i][1] += vbiasall[i][1];
      v[i][2] += vbiasall[i][2];
    }
}

/* ---------------------------------------------------------------------- */

double ComputeTempRotate::memory_usage()
{
  double bytes = maxbias * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "compute_temp_partial_cuda.h"
#include "compute_temp_partial_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "user_cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempPartialCuda::ComputeTempPartialCuda(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

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
  if (zflag && domain->dimension == 2)
    error->all(FLERR,"Compute temp/partial cannot use vz for 2d systemx");

  maxbias = 0;
  vbiasall = NULL;

  vector = new double[6];
  cu_t_vector = 0;
  cu_t_scalar = 0;
  cu_vbiasall=NULL;
  cudable=true;

}

/* ---------------------------------------------------------------------- */

ComputeTempPartialCuda::~ComputeTempPartialCuda()
{
  memory->destroy(vbiasall);
  delete [] vector;
  delete cu_t_vector;
  delete cu_t_scalar;
  delete cu_vbiasall;
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartialCuda::setup()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartialCuda::dof_compute()
{
  double natoms = group->count(igroup);
  int nper = xflag+yflag+zflag;
  dof = nper * natoms;
  dof -= (1.0*nper/domain->dimension)*fix_dof + extra_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputeTempPartialCuda::dof_remove(int i)
{
  int nper = xflag+yflag+zflag;
  return (domain->dimension - nper);
}

/* ---------------------------------------------------------------------- */

double ComputeTempPartialCuda::compute_scalar()
{
  if(cuda->begin_setup)
  {
          if(not cu_t_vector) cu_t_vector = new cCudaData<double, ENERGY_FLOAT, x> (t_vector,6);
          if(not cu_t_scalar) cu_t_scalar = new cCudaData<double, ENERGY_FLOAT, x> (&t_scalar,1);
    invoked_scalar = update->ntimestep;
    Cuda_ComputeTempPartialCuda_Scalar(&cuda->shared_data,groupbit,(ENERGY_FLOAT*) cu_t_scalar->dev_data(),xflag,yflag,zflag);
    cu_t_scalar->download();
  }
  else
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
        t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + zflag*v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + zflag*v[i][2]*v[i][2]) *
          mass[type[i]];
  }
  t_scalar=t;
  }

  MPI_Allreduce(&t_scalar,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  scalar *= tfactor;
  if(scalar>1e15)
  {
          cuda->cu_v->download();
          cuda->cu_x->download();
          cuda->cu_type->download();
    double **v = atom->v;
    double **x = atom->x;
    printf("Out of v-range atoms:  \n");
          for(int i=0;i<atom->nlocal;i++)
          if((v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])>1e5)
          printf("%i %i // %lf %lf %lf // %lf %lf %lf\n",atom->tag[i],atom->type[i],x[i][0], x[i][1], x[i][2],v[i][0], v[i][1], v[i][2]);
          error->all(FLERR,"Temperature out of range. Simulations will be abortet.\n");
  }
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartialCuda::compute_vector()
{
  int i;
  if(cuda->begin_setup)
  {
  if(not cu_t_vector) cu_t_vector = new cCudaData<double, ENERGY_FLOAT, x> (t_vector,6);
  if(not cu_t_scalar) cu_t_scalar = new cCudaData<double, ENERGY_FLOAT, x> (&t_scalar,1);

  invoked_vector = update->ntimestep;

  Cuda_ComputeTempPartialCuda_Vector(&cuda->shared_data,groupbit,(ENERGY_FLOAT*) cu_t_vector->dev_data(),xflag,yflag,zflag);
  cu_t_vector->download();
  }
  else
  {

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

  for (i = 0; i < 6; i++) t_vector[i]=t[i];
  }
  MPI_Allreduce(t_vector,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempPartialCuda::remove_bias(int i, double *v)
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
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempPartialCuda::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/partial:vbiasall");
        delete cu_vbiasall;
        cu_vbiasall = new cCudaData<double, V_FLOAT, yx> ((double*)vbiasall, atom->nmax, 3);
  }
  if(cuda->begin_setup)
  {
                  Cuda_ComputeTempPartialCuda_RemoveBiasAll(&cuda->shared_data,groupbit,xflag,yflag,zflag,cu_vbiasall->dev_data());
  }
  else
  {
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
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempPartialCuda::restore_bias(int i, double *v)
{
  if (!xflag) v[0] += vbias[0];
  if (!yflag) v[1] += vbias[1];
  if (!zflag) v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempPartialCuda::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if(cuda->begin_setup)
  {
                  Cuda_ComputeTempPartialCuda_RestoreBiasAll(&cuda->shared_data,groupbit,xflag,yflag,zflag,cu_vbiasall->dev_data());
  }
  else
  {

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
}

/* ---------------------------------------------------------------------- */

double ComputeTempPartialCuda::memory_usage()
{
  double bytes = maxbias * sizeof(double);
  return bytes;
}

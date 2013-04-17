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
#include "compute_temp_cuda.h"
#include "compute_temp_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "error.h"
#include "cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempCuda::ComputeTempCuda(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg != 3) error->all(FLERR,"Illegal compute temp/cuda command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[6];
  cu_t_vector = 0;
  cu_t_scalar = 0;
  cudable=true;

}

/* ---------------------------------------------------------------------- */

ComputeTempCuda::~ComputeTempCuda()
{
  delete [] vector;
  delete cu_t_vector;
  delete cu_t_scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCuda::setup()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempCuda::dof_compute()
{
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCuda::compute_scalar()
{
  if(cuda->begin_setup)
  {
          if(not cu_t_vector) cu_t_vector = new cCudaData<double, ENERGY_FLOAT, x> (t_vector,6);
          if(not cu_t_scalar) cu_t_scalar = new cCudaData<double, ENERGY_FLOAT, x> (&t_scalar,1);
    invoked_scalar = update->ntimestep;
    Cuda_ComputeTempCuda_Scalar(&cuda->shared_data,groupbit,(ENERGY_FLOAT*) cu_t_scalar->dev_data());
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
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
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

void ComputeTempCuda::compute_vector()
{
  int i;
  if(cuda->begin_setup)
  {
  if(not cu_t_vector) cu_t_vector = new cCudaData<double, ENERGY_FLOAT, x> (t_vector,6);
  if(not cu_t_scalar) cu_t_scalar = new cCudaData<double, ENERGY_FLOAT, x> (&t_scalar,1);

  invoked_vector = update->ntimestep;

  Cuda_ComputeTempCuda_Vector(&cuda->shared_data,groupbit,(ENERGY_FLOAT*) cu_t_vector->dev_data());
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
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  for (i = 0; i < 6; i++) t_vector[i]=t[i];
  }
  MPI_Allreduce(t_vector,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

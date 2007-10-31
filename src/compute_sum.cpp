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
#include "compute_sum.h"
#include "atom.h"
#include "modify.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSum::ComputeSum(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal compute sum command");

  // store pre-compute IDs

  npre = narg - 3;
  id_pre = new char*[npre];
  for (int i = 0; i < npre; i++) {
    int iarg = i + 3;
    int n = strlen(arg[iarg]) + 1;
    id_pre[i] = new char[n];
    strcpy(id_pre[i],arg[iarg]);
  }

  compute = new Compute*[npre];

  // all sub-computes must be peratom
  // check consistency of sub-computes for scalar & vector output

  int icompute;
  for (int i = 0; i < npre; i++) {
    icompute = modify->find_compute(id_pre[i]);
    if (icompute < 0)
      error->all("Could not find compute sum/atom pre-compute ID");
    if (modify->compute[icompute]->peratom_flag == 0)
      error->all("Compute sum compute is not a per-atom compute");
  }

  peratom_flag = 0;
  extensive = 0;

  icompute = modify->find_compute(id_pre[0]);
  int size = modify->compute[icompute]->size_peratom;
  if (size == 0) {
    scalar_flag = 1;
    vector = NULL;
  } else {
    vector_flag = 1;
    size_vector = size;
    vector = new double[size_vector];
  }

  for (int i = 1; i < npre; i++) {
    icompute = modify->find_compute(id_pre[i]);
    if (modify->compute[icompute]->size_peratom != size)
      error->all("Inconsistent sizes of compute sum compute quantities");
  }
}

/* ---------------------------------------------------------------------- */

ComputeSum::~ComputeSum()
{
  delete [] compute;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeSum::init()
{
  // set ptrs to Computes used as pre-computes by this compute

  for (int i = 0; i < npre; i++) {
    int icompute = modify->find_compute(id_pre[i]);
    if (icompute < 0)
      error->all("Could not find compute sum/atom pre-compute ID");
    compute[i] = modify->compute[icompute];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeSum::compute_scalar()
{
  int i;

  // invoke all the pre-computes
  // this is the only compute that does this
  // done b/c pre-computes are per-atom and this compute is not

  for (int icompute = 0; icompute < npre; icompute++)
    compute[icompute]->compute_peratom();

  // compute scalar quantity by summing over atom scalars

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  scalar = 0.0;

  for (int icompute = 0; icompute < npre; icompute++) {
    double *scalar_atom = compute[icompute]->scalar_atom;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) scalar += scalar_atom[i];
  }

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeSum::compute_vector()
{
  int i,j;

  // invoke all the pre-computes
  // this is the only compute that does this
  // done b/c pre-computes are per-atom and this compute is not

  for (int icompute = 0; icompute < npre; icompute++)
    compute[icompute]->compute_peratom();

  // compute vector quantity by summing over atom vectors

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (j = 0; j < size_vector; j++) vector[j] = 0.0;

  for (int icompute = 0; icompute < npre; icompute++) {
    double **vector_atom = compute[icompute]->vector_atom;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) 
	for (j = 0; j < size_vector; j++)
	  vector[j] += vector_atom[i][j];
  }
}

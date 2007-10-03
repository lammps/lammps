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
#include "compute_sum_atom.h"
#include "atom.h"
#include "modify.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSumAtom::ComputeSumAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal compute sum/atom command");

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

  // check consistency of set of pre-computes for scalar & vector output

  peratom_flag = 1;
  int icompute = modify->find_compute(id_pre[0]);
  size_peratom = modify->compute[icompute]->size_peratom;

  for (int i = 1; i < npre; i++) {
    icompute = modify->find_compute(id_pre[i]);
    if (icompute < 0)
      error->all("Could not find compute sum/atom pre-compute ID");
    if (modify->compute[icompute]->peratom_flag == 0)
      error->all("Compute sum/atom compute does not compute vector per atom");
    if (modify->compute[icompute]->size_peratom != size_peratom)
      error->all("Inconsistent sizes of compute sum/atom compute vectors");
  }

  nmax = 0;
  s_value = NULL;
  v_value = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSumAtom::~ComputeSumAtom()
{
  delete [] compute;
  memory->sfree(s_value);
  memory->destroy_2d_double_array(v_value);
}

/* ---------------------------------------------------------------------- */

void ComputeSumAtom::init()
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

void ComputeSumAtom::compute_peratom()
{
  int i,j,m;

  // grow sum array if necessary

  if (atom->nlocal > nmax) {
    nmax = atom->nmax;
    if (size_peratom == 0) {
      memory->sfree(s_value);
      s_value = (double *) 
	memory->smalloc(nmax*sizeof(double),"compute/sum/atom:s_value");
      scalar_atom = s_value;
    } else {
      memory->destroy_2d_double_array(v_value);
      v_value = memory->create_2d_double_array(nmax,size_peratom,
					       "compute/sum/atom:v_value");
      vector_atom = v_value;
    }
  }

  // sum over pre-computes
  // pre-computes of the pre-computes are not invoked

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (size_peratom == 0) {
    double *scalar = compute[0]->scalar_atom;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_value[i] = scalar[i];
      else s_value[i] = 0.0;
    
    for (m = 1; m < npre; m++) {
      scalar = compute[m]->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) s_value[i] += scalar[i];
    }

  } else {
    double **vector = compute[0]->vector_atom;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) 
	for (j = 0; j < size_peratom; j++)
	  v_value[i][j] = vector[i][j];
      else 
	for (j = 0; j < size_peratom; j++)
	  v_value[i][j] = 0.0;

    for (m = 1; m < npre; m++) {
      vector = compute[m]->vector_atom;
      for (j = 0; j < size_peratom; j++)
	if (mask[i] & groupbit) v_value[i][j] += vector[i][j];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

int ComputeSumAtom::memory_usage()
{
  int bytes = 0;
  if (size_peratom == 0) bytes = nmax * sizeof(double);
  else bytes = nmax*size_peratom * sizeof(double);
  return bytes;
}

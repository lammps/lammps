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
#include "compute_sum.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_PERATOM 4

enum{X,V,F,COMPUTE,FIX,VARIABLE};

/* ---------------------------------------------------------------------- */

ComputeSum::ComputeSum(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal compute sum command");

  // parse remaining values

  which = new int[narg-3];
  argindex = new int[narg-3];
  ids = new char*[narg-3];
  value2index = new int[narg-3];
  nvalues = 0;

  int iarg = 3;
  while (iarg < narg) {
    ids[nvalues] = NULL;

    if (strcmp(arg[iarg],"x") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"y") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"z") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 2;

    } else if ((strncmp(arg[iarg],"c_",2) == 0) || 
	       (strncmp(arg[iarg],"f_",2) == 0) || 
	       (strncmp(arg[iarg],"v_",2) == 0)) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Illegal compute sum command");
	argindex[nvalues] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all("Illegal compute sum command");

    iarg++;
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for compute sum does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
	error->all("Compute sum compute does not calculate per-atom values");
      if (argindex[i] == 0 && modify->compute[icompute]->size_peratom != 0)
	error->all("Compute sum compute does not calculate a per-atom scalar");
      if (argindex[i] && modify->compute[icompute]->size_peratom == 0)
	error->all("Compute sum compute does not calculate a per-atom vector");
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for compute sum does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
	error->all("Compute sum fix does not calculate per-atom values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_peratom != 0)
	error->all("Compute sum fix does not calculate a per-atom scalar");
      if (argindex[i] && modify->fix[ifix]->size_peratom == 0)
	error->all("Compute sum fix does not calculate a per-atom vector");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all("Variable name for compute sum does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
	error->all("Compute sum variable is not atom-style variable");
    }
  }

  // this compute produces either a scalar or vector

  if (nvalues == 1) {
    scalar_flag = 1;
    extscalar = 1;
    vector = NULL;
    onevec = NULL;
  } else {
    vector_flag = 1;
    size_vector = nvalues;
    extvector = 1;
    vector = new double[size_vector];
    onevec = new double[size_vector];
  }

  maxatom = 0;
  varatom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSum::~ComputeSum()
{
  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  delete [] vector;
  delete [] onevec;

  memory->sfree(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeSum::init()
{
  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
	error->all("Compute ID for compute sum does not exist");
      value2index[m] = icompute;
      
    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0) 
	error->all("Fix ID for compute sum does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0) 
	error->all("Variable name for compute sum does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ---------------------------------------------------------------------- */

double ComputeSum::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  double one = compute_one(0);
  MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeSum::compute_vector()
{
  invoked |= INVOKED_VECTOR;

  for (int m = 0; m < nvalues; m++) onevec[m] = compute_one(m);
  MPI_Allreduce(onevec,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

double ComputeSum::compute_one(int m)
{
  int i;

  // invoke the appropriate attribute,compute,fix,variable
  // compute scalar quantity by summing over atom scalars
  // only include atoms in group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = value2index[m];
  int j = argindex[m];
  double one = 0.0;

  if (which[m] == X) {
    double **x = atom->x;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += x[i][j];
  } else if (which[m] == V) {
    double **v = atom->v;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += v[i][j];
  } else if (which[m] == F) {
    double **f = atom->f;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += f[i][j];
    
  // invoke compute if not previously invoked

  } else if (which[m] == COMPUTE) {
    if (!(modify->compute[n]->invoked & INVOKED_PERATOM))
      modify->compute[n]->compute_peratom();

    if (j == 0) {
      double *compute_scalar = modify->compute[n]->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) one += compute_scalar[i];
    } else {
      double **compute_vector = modify->compute[n]->vector_atom;
      int jm1 = j - 1;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) one += compute_vector[i][jm1];
    }

  // access fix fields, check if frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[n]->peratom_freq)
      error->all("Fix used in compute sum not computed at compatible time");

    if (j == 0) {
      double *fix_scalar = modify->fix[n]->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) one += fix_scalar[i];
    } else {
      double **fix_vector = modify->fix[n]->vector_atom;
      int jm1 = j - 1;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) one += fix_vector[i][jm1];
    }
    
  // evaluate atom-style variable

  } else if (which[m] == VARIABLE) {
    if (nlocal > maxatom) {
      maxatom = atom->nmax;
      memory->sfree(varatom);
      varatom =	(double *) 
	memory->smalloc(maxatom*sizeof(double),"compute/sum:varatom");
    }

    input->variable->compute_atom(n,igroup,varatom,1,0);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += varatom[i];
  }

  return one;
}

/* ----------------------------------------------------------------------
   memory usage of varatom
------------------------------------------------------------------------- */

double ComputeSum::memory_usage()
{
  double bytes = maxatom * sizeof(double);
  return bytes;
}

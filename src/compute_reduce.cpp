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
#include "compute_reduce.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SUM,MINN,MAXX};
enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduce::ComputeReduce(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  int iarg;
  if (strcmp(style,"reduce") == 0) {
    if (narg < 5) error->all("Illegal compute reduce command");
    iarg = 3;
  } else if (strcmp(style,"reduce/region") == 0) {
    if (narg < 6) error->all("Illegal compute reduce/region command");
    iregion = domain->find_region(arg[3]);
    if (iregion == -1) error->all("Compute reduce region ID does not exist");
    iarg = 4;
  }

  if (strcmp(arg[iarg],"sum") == 0) mode = SUM;
  else if (strcmp(arg[iarg],"min") == 0) mode = MINN;
  else if (strcmp(arg[iarg],"max") == 0) mode = MAXX;
  else error->all("Illegal compute reduce command");
  iarg++;

  // parse remaining values

  which = new int[narg-4];
  argindex = new int[narg-4];
  ids = new char*[narg-4];
  value2index = new int[narg-4];
  nvalues = 0;

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
	  error->all("Illegal compute reduce command");
	argindex[nvalues] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all("Illegal compute reduce command");

    iarg++;
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for compute reduce does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
	error->all("Compute reduce compute does not calculate per-atom values");
      if (argindex[i] == 0 && modify->compute[icompute]->size_peratom != 0)
	error->all("Compute reduce compute does not calculate a per-atom scalar");
      if (argindex[i] && modify->compute[icompute]->size_peratom == 0)
	error->all("Compute reduce compute does not calculate a per-atom vector");
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for compute reduce does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
	error->all("Compute reduce fix does not calculate per-atom values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_peratom != 0)
	error->all("Compute reduce fix does not calculate a per-atom scalar");
      if (argindex[i] && modify->fix[ifix]->size_peratom == 0)
	error->all("Compute reduce fix does not calculate a per-atom vector");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all("Variable name for compute reduce does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
	error->all("Compute reduce variable is not atom-style variable");
    }
  }

  // this compute produces either a scalar or vector

  if (nvalues == 1) {
    scalar_flag = 1;
    if (mode == SUM) extscalar = 1;
    else extscalar = 0;
    vector = NULL;
    onevec = NULL;
  } else {
    vector_flag = 1;
    size_vector = nvalues;
    if (mode == SUM) extvector = 1;
    else extvector = 0;
    vector = new double[size_vector];
    onevec = new double[size_vector];
  }

  maxatom = 0;
  varatom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeReduce::~ComputeReduce()
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

void ComputeReduce::init()
{
  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
	error->all("Compute ID for compute reduce does not exist");
      value2index[m] = icompute;
      
    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0) 
	error->all("Fix ID for compute reduce does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0) 
	error->all("Variable name for compute reduce does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ---------------------------------------------------------------------- */

double ComputeReduce::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double one = compute_one(0);
  if (mode == SUM)
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  else if (mode == MINN)
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MIN,world);
  else if (mode == MAXX)
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MAX,world);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::compute_vector()
{
  invoked_vector = update->ntimestep;

  for (int m = 0; m < nvalues; m++) onevec[m] = compute_one(m);
  if (mode == SUM)
    MPI_Allreduce(onevec,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);
  else if (mode == MINN)
    MPI_Allreduce(onevec,vector,size_vector,MPI_DOUBLE,MPI_MIN,world);
  else if (mode == MAXX)
    MPI_Allreduce(onevec,vector,size_vector,MPI_DOUBLE,MPI_MAX,world);
}

/* ---------------------------------------------------------------------- */

double ComputeReduce::compute_one(int m)
{
  int i;

  // invoke the appropriate attribute,compute,fix,variable
  // compute scalar quantity by summing over atom scalars
  // only include atoms in group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = value2index[m];
  int j = argindex[m];

  double one;
  if (mode == SUM) one = 0.0;
  else if (mode == MINN) one = BIG;
  else if (mode == MAXX) one = -BIG;

  if (which[m] == X) {
    double **x = atom->x;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) combine(one,x[i][j]);
  } else if (which[m] == V) {
    double **v = atom->v;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) combine(one,v[i][j]);
  } else if (which[m] == F) {
    double **f = atom->f;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) combine(one,f[i][j]);
    
  // invoke compute if not previously invoked

  } else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[n];
    if (!(compute->invoked_flag & INVOKED_PERATOM)) {
      compute->compute_peratom();
      compute->invoked_flag |= INVOKED_PERATOM;
    }

    if (j == 0) {
      double *compute_scalar = compute->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,compute_scalar[i]);
    } else {
      double **compute_vector = compute->vector_atom;
      int jm1 = j - 1;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,compute_vector[i][jm1]);
    }

  // access fix fields, check if frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[n]->peratom_freq)
      error->all("Fix used in compute reduce not computed at compatible time");

    if (j == 0) {
      double *fix_scalar = modify->fix[n]->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,fix_scalar[i]);
    } else {
      double **fix_vector = modify->fix[n]->vector_atom;
      int jm1 = j - 1;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,fix_vector[i][jm1]);
    }
    
  // evaluate atom-style variable

  } else if (which[m] == VARIABLE) {
    if (nlocal > maxatom) {
      maxatom = atom->nmax;
      memory->sfree(varatom);
      varatom =	(double *) 
	memory->smalloc(maxatom*sizeof(double),"compute/reduce:varatom");
    }

    input->variable->compute_atom(n,igroup,varatom,1,0);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) combine(one,varatom[i]);
  }

  return one;
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
------------------------------------------------------------------------- */

void ComputeReduce::combine(double &one, double two)
{
  if (mode == SUM) one += two;
  else if (mode == MINN) one = MIN(one,two);
  else if (mode == MAXX) one = MAX(one,two);
}

/* ----------------------------------------------------------------------
   memory usage of varatom
------------------------------------------------------------------------- */

double ComputeReduce::memory_usage()
{
  double bytes = maxatom * sizeof(double);
  return bytes;
}

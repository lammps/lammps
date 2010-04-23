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
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SUM,MINN,MAXX,AVE};
enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{GLOBAL,PERATOM,LOCAL};

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

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
  else if (strcmp(arg[iarg],"ave") == 0) mode = AVE;
  else error->all("Illegal compute reduce command");
  iarg++;

  MPI_Comm_rank(world,&me);

  // parse remaining values until one isn't recognized

  which = new int[narg-4];
  argindex = new int[narg-4];
  flavor = new int[narg-4];
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

    } else if (strncmp(arg[iarg],"c_",2) == 0 || 
	       strncmp(arg[iarg],"f_",2) == 0 || 
	       strncmp(arg[iarg],"v_",2) == 0) {
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

    } else break;

    iarg++;
  }

  // optional args

  replace = new int[nvalues];
  for (int i = 0; i < nvalues; i++) replace[i] = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"replace") == 0) {
      if (iarg+3 > narg) error->all("Illegal compute reduce command");
      if (mode != MINN && mode != MAXX)
	error->all("Compute reduce replace requires min or max mode");
      int col1 = atoi(arg[iarg+1]) - 1;
      int col2 = atoi(arg[iarg+2]) - 1;
      if (col1 < 0 || col1 >= nvalues || col2 < 0 || col2 >= nvalues)
	error->all("Illegal compute reduce command");
      if (col1 == col2)	error->all("Illegal compute reduce command");
      if (replace[col1] >= 0 || replace[col2] >= 0)
	error->all("Invalid replace values in compute reduce"); 
      replace[col1] = col2;
      iarg += 3;
    } else error->all("Illegal compute reduce command");
  }

  // delete replace if not set

  int flag = 0;
  for (int i = 0; i < nvalues; i++)
    if (replace[i] >= 0) flag = 1;
  if (!flag) {
    delete [] replace;
    replace = NULL;
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for compute reduce does not exist");
      if (modify->compute[icompute]->vector_flag || 
	  modify->compute[icompute]->array_flag) {
	flavor[i] = GLOBAL;
	if (argindex[i] == 0 && 
	    modify->compute[icompute]->vector_flag == 0)
	  error->all("Compute reduce compute does not "
		     "calculate a global vector");
	if (argindex[i] && modify->compute[icompute]->array_flag == 0)
	  error->all("Compute reduce compute does not "
		     "calculate a global array");
	if (argindex[i] && 
	    argindex[i] > modify->compute[icompute]->size_array_cols)
	  error->all("Compute reduce compute array is accessed out-of-range");
      } else if (modify->compute[icompute]->peratom_flag) {
	flavor[i] = PERATOM;
	if (argindex[i] == 0 && 
	    modify->compute[icompute]->size_peratom_cols != 0)
	  error->all("Compute reduce compute does not "
		     "calculate a per-atom vector");
	if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
	  error->all("Compute reduce compute does not "
		     "calculate a per-atom array");
	if (argindex[i] && 
	    argindex[i] > modify->compute[icompute]->size_peratom_cols)
	  error->all("Compute reduce compute array is accessed out-of-range");
      } else if (modify->compute[icompute]->local_flag) {
	flavor[i] = LOCAL;
	if (argindex[i] == 0 && 
	    modify->compute[icompute]->size_local_cols != 0)
	  error->all("Compute reduce compute does not "
		     "calculate a local vector");
	if (argindex[i] && modify->compute[icompute]->size_local_cols == 0)
	  error->all("Compute reduce compute does not "
		     "calculate a local array");
	if (argindex[i] && 
	    argindex[i] > modify->compute[icompute]->size_local_cols)
	  error->all("Compute reduce compute array is accessed out-of-range");
      }

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for compute reduce does not exist");
      if (modify->fix[ifix]->vector_flag || 
	  modify->fix[ifix]->array_flag) {
	flavor[i] = GLOBAL;
	if (argindex[i] == 0 && 
	    modify->fix[ifix]->vector_flag == 0)
	  error->all("Compute reduce fix does not "
		     "calculate a global vector");
	if (argindex[i] && modify->fix[ifix]->array_flag == 0)
	  error->all("Compute reduce fix does not "
		     "calculate a global array");
	if (argindex[i] && 
	    argindex[i] > modify->fix[ifix]->size_array_cols)
	  error->all("Compute reduce fix array is accessed out-of-range");
      } else if (modify->fix[ifix]->peratom_flag) {
	flavor[i] = PERATOM;
	if (argindex[i] == 0 && 
	    modify->fix[ifix]->size_peratom_cols != 0)
	  error->all("Compute reduce fix does not "
		     "calculate a per-atom vector");
	if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
	  error->all("Compute reduce fix does not "
		     "calculate a per-atom array");
	if (argindex[i] && 
	    argindex[i] > modify->fix[ifix]->size_peratom_cols)
	  error->all("Compute reduce fix array is accessed out-of-range");
      } else if (modify->fix[ifix]->local_flag) {
	flavor[i] = LOCAL;
	if (argindex[i] == 0 && 
	    modify->fix[ifix]->size_local_cols != 0)
	  error->all("Compute reduce fix does not "
		     "calculate a local vector");
	if (argindex[i] && modify->fix[ifix]->size_local_cols == 0)
	  error->all("Compute reduce fix does not "
		     "calculate a local array");
	if (argindex[i] && 
	    argindex[i] > modify->fix[ifix]->size_local_cols)
	  error->all("Compute reduce fix array is accessed out-of-range");
      }

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all("Variable name for compute reduce does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
	error->all("Compute reduce variable is not atom-style variable");
      flavor[i] = PERATOM;
    }
  }

  // this compute produces either a scalar or vector

  if (nvalues == 1) {
    scalar_flag = 1;
    if (mode == SUM) extscalar = 1;
    else extscalar = 0;
    vector = onevec = NULL;
    indices = owner = NULL;
  } else {
    vector_flag = 1;
    size_vector = nvalues;
    if (mode == SUM) extvector = 1;
    else extvector = 0;
    vector = new double[size_vector];
    onevec = new double[size_vector];
    indices = new int[size_vector];
    owner = new int[size_vector];
  }

  maxatom = 0;
  varatom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeReduce::~ComputeReduce()
{
  delete [] which;
  delete [] argindex;
  delete [] flavor;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;
  delete [] replace;

  delete [] vector;
  delete [] onevec;
  delete [] indices;
  delete [] owner;

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

  double one = compute_one(0,-1);

  if (mode == SUM) {
    if (flavor[0] == GLOBAL)
      scalar = one;
    else
      MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  } else if (mode == MINN) {
    if (flavor[0] == GLOBAL)
      scalar = one;
    else
      MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MIN,world);
  } else if (mode == MAXX) {
    if (flavor[0] == GLOBAL)
      scalar = one;
    else
      MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MAX,world);
  } else if (mode == AVE) {
    if (flavor[0] == GLOBAL)
      scalar = one;
    else
      MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    scalar /= count(0);
  }

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::compute_vector()
{
  invoked_vector = update->ntimestep;

  for (int m = 0; m < nvalues; m++)
    if (!replace || replace[m] < 0) {
      onevec[m] = compute_one(m,-1);
      indices[m] = index;
    }

  if (mode == SUM) {
    for (int m = 0; m < nvalues; m++) {
      if (flavor[m] == GLOBAL)
	vector[m] = onevec[m];
      else
	MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
    }

  } else if (mode == MINN) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++) {
	if (flavor[m] == GLOBAL)
	  vector[m] = onevec[m];
	else
	  MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_MIN,world);
      }
    } else {
      for (int m = 0; m < nvalues; m++)
	if (replace[m] < 0) {
	  if (flavor[m] == GLOBAL) {
	    vector[m] = onevec[m];
	    owner[m] = me;
	  } else {
	    pairme.value = onevec[m];
	    pairme.proc = me;
	    MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MINLOC,world);
	    vector[m] = pairall.value;
	    owner[m] = pairall.proc;
	  }
	}
      for (int m = 0; m < nvalues; m++)
	if (replace[m] >= 0) {
	  if (me == owner[replace[m]])
	    vector[m] = compute_one(m,indices[replace[m]]);
	  MPI_Bcast(&vector[m],1,MPI_DOUBLE,owner[replace[m]],world);
	}
    }

  } else if (mode == MAXX) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++) {
	if (flavor[m] == GLOBAL)
	  vector[m] = onevec[m];
	else
	  MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_MAX,world);
      }
    } else {
      for (int m = 0; m < nvalues; m++)
	if (replace[m] < 0) {
	  if (flavor[m] == GLOBAL) {
	    vector[m] = onevec[m];
	    owner[m] = me;
	  } else {
	    pairme.value = onevec[m];
	    pairme.proc = me;
	    MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MAXLOC,world);
	    vector[m] = pairall.value;
	    owner[m] = pairall.proc;
	  }
	}
      for (int m = 0; m < nvalues; m++)
	if (replace[m] >= 0) {
	  if (me == owner[replace[m]])
	    vector[m] = compute_one(m,indices[replace[m]]);
	  MPI_Bcast(&vector[m],1,MPI_DOUBLE,owner[replace[m]],world);
	}
    }

  } else if (mode == AVE) {
    for (int m = 0; m < nvalues; m++) {
      if (flavor[m] == GLOBAL)
	vector[m] = onevec[m];
      else
	MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
      vector[m] /= count(m);
    }
  }
}

/* ----------------------------------------------------------------------
   calculate reduced value for one input M and return it
   if flag = -1:
     sum/min/max/ave all values in vector
     for per-atom quantities, limit to atoms in group
     if mode = MIN or MAX, also set index to which vector value wins
   if flag >= 0: simply return vector[flag]
------------------------------------------------------------------------- */

double ComputeReduce::compute_one(int m, int flag)
{
  int i;

  // invoke the appropriate attribute,compute,fix,variable
  // for flag = -1, compute scalar quantity by scanning over atom properties
  // only include atoms in group for atom properties and per-atom quantities

  index = -1;
  int n = value2index[m];
  int j = argindex[m];

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double one;
  if (mode == SUM) one = 0.0;
  else if (mode == MINN) one = BIG;
  else if (mode == MAXX) one = -BIG;
  else if (mode == AVE) one = 0.0;

  if (which[m] == X) {
    double **x = atom->x;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,x[i][j],i);
    } else one = x[flag][j];
  } else if (which[m] == V) {
    double **v = atom->v;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,v[i][j],i);
    } else one = v[flag][j];
  } else if (which[m] == F) {
    double **f = atom->f;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,f[i][j],i);
    } else one = f[flag][j];
    
  // invoke compute if not previously invoked

  } else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[n];
    
    if (flavor[m] == GLOBAL) {
      if (j == 0) {
	if (!(compute->invoked_flag & INVOKED_VECTOR)) {
	  compute->compute_vector();
	  compute->invoked_flag |= INVOKED_VECTOR;
	}
	double *compute_vector = compute->vector;
	int n = compute->size_vector;
	if (flag < 0) 
	  for (i = 0; i < n; i++)
	    combine(one,compute_vector[i],i);
	else one = compute_vector[flag];
      } else {
	if (!(compute->invoked_flag & INVOKED_ARRAY)) {
	  compute->compute_array();
	  compute->invoked_flag |= INVOKED_ARRAY;
	}
	double **compute_array = compute->array;
	int n = compute->size_array_rows;
	int jm1 = j - 1;
	if (flag < 0)
	  for (i = 0; i < n; i++)
	    combine(one,compute_array[i][jm1],i);
	else one = compute_array[flag][jm1];
      }

    } else if (flavor[m] == PERATOM) {
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
	compute->compute_peratom();
	compute->invoked_flag |= INVOKED_PERATOM;
      }

      if (j == 0) {
	double *compute_vector = compute->vector_atom;
	int n = nlocal;
	if (flag < 0) {
	  for (i = 0; i < n; i++)
	    if (mask[i] & groupbit) combine(one,compute_vector[i],i);
	} else one = compute_vector[flag];
      } else {
	double **compute_array = compute->array_atom;
	int n = nlocal;
	int jm1 = j - 1;
	if (flag < 0) {
	  for (i = 0; i < n; i++)
	    if (mask[i] & groupbit) combine(one,compute_array[i][jm1],i);
	} else one = compute_array[flag][jm1];
      }

    } else if (flavor[m] == LOCAL) {
      if (!(compute->invoked_flag & INVOKED_LOCAL)) {
	compute->compute_local();
	compute->invoked_flag |= INVOKED_LOCAL;
      }

      if (j == 0) {
	double *compute_vector = compute->vector_local;
	int n = compute->size_local_rows;
	if (flag < 0)
	  for (i = 0; i < n; i++)
	    combine(one,compute_vector[i],i);
	else one = compute_vector[flag];
      } else {
	double **compute_array = compute->array_local;
	int n = compute->size_local_rows;
	int jm1 = j - 1;
	if (flag < 0)
	  for (i = 0; i < n; i++)
	    combine(one,compute_array[i][jm1],i);
	else one = compute_array[flag][jm1];
      }
    }

  // access fix fields, check if fix frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[n]->peratom_freq)
      error->all("Fix used in compute reduce not computed at compatible time");
    Fix *fix = modify->fix[n];

    if (flavor[m] == GLOBAL) {
      if (j == 0) {
	int n = fix->size_vector;
	if (flag < 0)
	  for (i = 0; i < n; i++)
	    combine(one,fix->compute_vector(i),i);
	else one = fix->compute_vector(flag);
      } else {
	int n = fix->size_array_rows;
	int jm1 = j - 1;
	if (flag < 0)
	  for (i = 0; i < nlocal; i++)
	    combine(one,fix->compute_array(i,jm1),i);
	else one = fix->compute_array(flag,jm1);
      }

    } else if (flavor[m] == PERATOM) {
      if (j == 0) {
	double *fix_vector = fix->vector_atom;
	int n = nlocal;
	if (flag < 0) {
	  for (i = 0; i < n; i++)
	    if (mask[i] & groupbit) combine(one,fix_vector[i],i);
	} else one = fix_vector[flag];
      } else {
	double **fix_array = fix->array_atom;
	int n = nlocal;
	int jm1 = j - 1;
	if (flag < 0) {
	  for (i = 0; i < nlocal; i++)
	    if (mask[i] & groupbit) combine(one,fix_array[i][jm1],i);
	} else one = fix_array[flag][jm1];
      }

    } else if (flavor[m] == LOCAL) {
      if (j == 0) {
	double *fix_vector = fix->vector_local;
	int n = fix->size_local_rows;
	if (flag < 0)
	  for (i = 0; i < n; i++)
	    combine(one,fix_vector[i],i);
	else one = fix_vector[flag];
      } else {
	double **fix_array = fix->array_local;
	int n = fix->size_local_rows;
	int jm1 = j - 1;
	if (flag < 0)
	  for (i = 0; i < n; i++)
	    combine(one,fix_array[i][jm1],i);
	else one = fix_array[flag][jm1];
      }
    }
    
  // evaluate atom-style variable

  } else if (which[m] == VARIABLE) {
    if (nlocal > maxatom) {
      maxatom = atom->nmax;
      memory->sfree(varatom);
      varatom =	(double *) 
	memory->smalloc(maxatom*sizeof(double),"reduce:varatom");
    }

    input->variable->compute_atom(n,igroup,varatom,1,0);
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) combine(one,varatom[i],i);
    } else one = varatom[flag];
  }

  return one;
}

/* ---------------------------------------------------------------------- */

double ComputeReduce::count(int m)
{
  int n = value2index[m];
  int j = argindex[m];

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (which[m] == X || which[m] == V || which[m] == F)
    return group->count(igroup);
  else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[n];
    if (flavor[m] == GLOBAL) {
      if (j == 0) return(1.0*compute->size_vector);
      else return(1.0*compute->size_array_rows);
    } else if (flavor[m] == PERATOM) {
      return group->count(igroup);
    } else if (flavor[m] == LOCAL) {
      double ncount = compute->size_local_rows;
      double ncountall;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_DOUBLE,MPI_SUM,world);
      return ncountall;
    }
  } else if (which[m] == FIX) {
    Fix *fix = modify->fix[n];
    if (flavor[m] == GLOBAL) {
      if (j == 0) return(1.0*fix->size_vector);
      else return(1.0*fix->size_array_rows);
    } else if (flavor[m] == PERATOM) {
      return group->count(igroup);
    } else if (flavor[m] == LOCAL) {
      double ncount = fix->size_local_rows;
      double ncountall;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_DOUBLE,MPI_SUM,world);
      return ncountall;
    }
  } else if (which[m] == VARIABLE)
    return group->count(igroup);

  return 0.0;
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
   for MIN/MAX, also update index with winner
------------------------------------------------------------------------- */

void ComputeReduce::combine(double &one, double two, int i)
{
  if (mode == SUM || mode == AVE) one += two;
  else if (mode == MINN) {
    if (two < one) {
      one = two;
      index = i;
    }
  } else if (mode == MAXX) {
    if (two > one) {
      one = two;
      index = i;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of varatom
------------------------------------------------------------------------- */

double ComputeReduce::memory_usage()
{
  double bytes = maxatom * sizeof(double);
  return bytes;
}

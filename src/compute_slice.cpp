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

#include "stdlib.h"
#include "string.h"
#include "compute_slice.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

enum{COMPUTE,FIX};

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

/* ---------------------------------------------------------------------- */

ComputeSlice::ComputeSlice(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal compute slice command");

  MPI_Comm_rank(world,&me);

  nstart = force->inumeric(FLERR,arg[3]);
  nstop = force->inumeric(FLERR,arg[4]);
  nskip = force->inumeric(FLERR,arg[5]);

  if (nstart < 1 || nstop < nstart || nskip < 1)
    error->all(FLERR,"Illegal compute slice command");

  // parse remaining values until one isn't recognized

  which = new int[narg-6];
  argindex = new int[narg-6];
  ids = new char*[narg-6];
  value2index = new int[narg-6];
  nvalues = 0;

  for (int iarg = 6; iarg < narg; iarg++) {
    if (strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal compute slice command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all(FLERR,"Illegal compute slice command");
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute slice does not exist");
      if (modify->compute[icompute]->vector_flag) {
        if (argindex[i])
          error->all(FLERR,"Compute slice compute does not calculate a global array");
        if (nstop > modify->compute[icompute]->size_vector)
          error->all(FLERR,"Compute slice compute vector is accessed out-of-range");
      } else if (modify->compute[icompute]->array_flag) {
        if (argindex[i] == 0)
          error->all(FLERR,"Compute slice compute does not calculate a global vector");
        if (argindex[i] > modify->compute[icompute]->size_array_cols)
          error->all(FLERR,"Compute slice compute array is accessed out-of-range");
        if (nstop > modify->compute[icompute]->size_array_rows)
          error->all(FLERR,"Compute slice compute array is accessed out-of-range");
      } else error->all(FLERR,"Compute slice compute does not calculate "
                        "global vector or array");
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute slice does not exist");
      if (modify->fix[ifix]->vector_flag) {
        if (argindex[i])
          error->all(FLERR,"Compute slice fix does not calculate a global array");
        if (nstop > modify->fix[ifix]->size_vector)
          error->all(FLERR,"Compute slice fix vector is accessed out-of-range");
      } else if (modify->fix[ifix]->array_flag) {
        if (argindex[i] == 0)
          error->all(FLERR,"Compute slice fix does not calculate a global vector");
        if (argindex[i] > modify->fix[ifix]->size_array_cols)
          error->all(FLERR,"Compute slice fix array is accessed out-of-range");
        if (nstop > modify->fix[ifix]->size_array_rows)
          error->all(FLERR,"Compute slice fix array is accessed out-of-range");
      } else error->all(FLERR,"Compute slice fix does not calculate "
                        "global vector or array");
    }
  }

  // this compute produces either a vector or array
  // for vector, set intensive/extensive to mirror input values
  // for array, set intensive if all input values are intensive, else extensive

  vector = NULL;
  array = NULL;
  extlist = NULL;

  if (nvalues == 1) {
    vector_flag = 1;
    size_vector = (nstop-nstart) / nskip;
    memory->create(vector,size_vector,"slice:vector");

    if (which[0] == COMPUTE) {
      int icompute = modify->find_compute(ids[0]);
      if (argindex[0] == 0) {
        extvector = modify->compute[icompute]->extvector;
        if (modify->compute[icompute]->extvector == -1) {
          extlist = new int[size_vector];
          int j = 0;
          for (int i = nstart; i < nstop; i += nskip)
            extlist[j++] = modify->compute[icompute]->extlist[i-1];
        }
      } else extvector = modify->compute[icompute]->extarray;
    } else if (which[0] == FIX) {
      int ifix = modify->find_fix(ids[0]);
      if (argindex[0] == 0) {
        extvector = modify->fix[ifix]->extvector;
        if (modify->fix[ifix]->extvector == -1) {
          extlist = new int[size_vector];
          int j = 0;
          for (int i = nstart; i < nstop; i += nskip)
            extlist[j++] = modify->fix[ifix]->extlist[i-1];
        }
      } else extvector = modify->fix[ifix]->extarray;
    }

  } else {
    array_flag = 1;
    size_array_rows = (nstop-nstart) / nskip;
    size_array_cols = nvalues;
    memory->create(array,size_array_rows,size_array_cols,"slice:array");

    extarray = 0;
    for (int i = 0; i < nvalues; i++) {
      if (which[i] == COMPUTE) {
        int icompute = modify->find_compute(ids[i]);
        if (argindex[i] == 0) {
          if (modify->compute[icompute]->extvector == 1) extarray = 1;
          if (modify->compute[icompute]->extvector == -1) {
            for (int j = 0; j < modify->compute[icompute]->size_vector; j++)
              if (modify->compute[icompute]->extlist[j]) extarray = 1;
          }
        } else {
          if (modify->compute[icompute]->extarray) extarray = 1;
        }
      } else if (which[i] == FIX) {
        int ifix = modify->find_fix(ids[i]);
        if (argindex[i] == 0) {
          if (modify->fix[ifix]->extvector == 1) extarray = 1;
          if (modify->fix[ifix]->extvector == -1) {
            for (int j = 0; j < modify->fix[ifix]->size_vector; j++)
              if (modify->fix[ifix]->extlist[j]) extarray = 1;
          }
        } else {
          if (modify->fix[ifix]->extarray) extarray = 1;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

ComputeSlice::~ComputeSlice()
{
  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::init()
{
  // set indices and check validity of all computes,fixes

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute slice does not exist");
      value2index[m] = icompute;
    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute slice does not exist");
      value2index[m] = ifix;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::compute_vector()
{
  invoked_vector = update->ntimestep;

  extract_one(0,vector,1);
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::compute_array()
{
  invoked_array = update->ntimestep;

  for (int m = 0; m < nvalues; m++)
    extract_one(0,&array[m][0],nvalues);
}

/* ----------------------------------------------------------------------
   calculate sliced value for one input M and return it in vec
   vec may be array so that returned values are with stride
------------------------------------------------------------------------- */

void ComputeSlice::extract_one(int m, double *vec, int stride)
{
  int i,j;

  // invoke the appropriate compute if needed

  if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[value2index[m]];

    if (argindex[m] == 0) {
      if (!(compute->invoked_flag & INVOKED_VECTOR)) {
        compute->compute_vector();
        compute->invoked_flag |= INVOKED_VECTOR;
      }
      double *cvector = compute->vector;
      j = 0;
      for (i = nstart; i < nstop; i += nskip) {
        vec[j] = cvector[i-1];
        j += stride;
      }

    } else {
      if (!(compute->invoked_flag & INVOKED_ARRAY)) {
        compute->compute_array();
        compute->invoked_flag |= INVOKED_ARRAY;
      }
      double **carray = compute->array;
      int icol = argindex[m]-1;
      j = 0;
      for (i = nstart; i < nstop; i += nskip) {
        vec[j] = carray[i-1][icol];
        j += stride;
      }
    }

  // access fix fields, check if fix frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[value2index[m]]->global_freq)
      error->all(FLERR,"Fix used in compute slice not computed at compatible time");
    Fix *fix = modify->fix[value2index[m]];

    if (argindex[m] == 0) {
      j = 0;
      for (i = nstart; i < nstop; i += nskip) {
        vec[j] = fix->compute_vector(i-1);
        j += stride;
      }
    } else {
      int icol = argindex[m]-1;
      j = 0;
      for (i = nstart; i < nstop; i += nskip) {
        vec[j] = fix->compute_array(i-1,icol);
        j += stride;
      }
    }
  }
}

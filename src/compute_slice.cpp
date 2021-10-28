// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_slice.h"

#include "arg_info.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSlice::ComputeSlice(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  nvalues(0), which(nullptr), argindex(nullptr), value2index(nullptr), ids(nullptr)
{
  if (narg < 7) error->all(FLERR,"Illegal compute slice command");

  MPI_Comm_rank(world,&me);

  nstart = utils::inumeric(FLERR,arg[3],false,lmp);
  nstop = utils::inumeric(FLERR,arg[4],false,lmp);
  nskip = utils::inumeric(FLERR,arg[5],false,lmp);

  if (nstart < 1 || nstop < nstart || nskip < 1)
    error->all(FLERR,"Illegal compute slice command");

  // parse remaining values until one isn't recognized

  which = new int[narg-6];
  argindex = new int[narg-6];
  ids = new char*[narg-6];
  value2index = new int[narg-6];
  nvalues = 0;

  for (int iarg = 6; iarg < narg; iarg++) {
    ArgInfo argi(arg[iarg]);

    which[nvalues] = argi.get_type();
    argindex[nvalues] = argi.get_index1();
    ids[nvalues] = argi.copy_name();

    if ((which[nvalues] == ArgInfo::UNKNOWN) || (which[nvalues] == ArgInfo::NONE)
        || (argi.get_dim() > 1))
      error->all(FLERR,"Illegal compute slice command");

    nvalues++;
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute slice does not exist");
      if (modify->compute[icompute]->vector_flag) {
        if (argindex[i])
          error->all(FLERR,"Compute slice compute does not "
                     "calculate a global array");
        if (nstop > modify->compute[icompute]->size_vector)
          error->all(FLERR,"Compute slice compute vector is "
                     "accessed out-of-range");
      } else if (modify->compute[icompute]->array_flag) {
        if (argindex[i] == 0)
          error->all(FLERR,"Compute slice compute does not "
                     "calculate a global vector");
        if (argindex[i] > modify->compute[icompute]->size_array_cols)
          error->all(FLERR,"Compute slice compute array is "
                     "accessed out-of-range");
        if (nstop > modify->compute[icompute]->size_array_rows)
          error->all(FLERR,"Compute slice compute array is "
                     "accessed out-of-range");
      } else error->all(FLERR,"Compute slice compute does not calculate "
                        "global vector or array");

    } else if (which[i] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute slice does not exist");
      if (modify->fix[ifix]->vector_flag) {
        if (argindex[i])
          error->all(FLERR,"Compute slice fix does not "
                     "calculate a global array");
        if (nstop > modify->fix[ifix]->size_vector)
          error->all(FLERR,"Compute slice fix vector is accessed out-of-range");
      } else if (modify->fix[ifix]->array_flag) {
        if (argindex[i] == 0)
          error->all(FLERR,"Compute slice fix does not "
                     "calculate a global vector");
        if (argindex[i] > modify->fix[ifix]->size_array_cols)
          error->all(FLERR,"Compute slice fix array is accessed out-of-range");
        if (nstop > modify->fix[ifix]->size_array_rows)
          error->all(FLERR,"Compute slice fix array is accessed out-of-range");
      } else error->all(FLERR,"Compute slice fix does not calculate "
                        "global vector or array");

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute slice does not exist");
      if (argindex[i] == 0 && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Compute slice variable is not vector-style variable");
      if (argindex[i])
        error->all(FLERR,"Compute slice vector variable cannot be indexed");
    }
  }

  // this compute produces either a vector or array
  // for vector, set intensive/extensive to mirror input values
  // for array, set intensive if all input values are intensive, else extensive

  if (nvalues == 1) {
    vector_flag = 1;
    size_vector = (nstop-nstart) / nskip;
    memory->create(vector,size_vector,"slice:vector");

    if (which[0] == ArgInfo::COMPUTE) {
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
    } else if (which[0] == ArgInfo::FIX) {
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
    } else if (which[0] == ArgInfo::VARIABLE) {
      extvector = 0;
    }

  } else {
    array_flag = 1;
    size_array_rows = (nstop-nstart) / nskip;
    size_array_cols = nvalues;
    memory->create(array,size_array_rows,size_array_cols,"slice:array");

    extarray = 0;
    for (int i = 0; i < nvalues; i++) {
      if (which[i] == ArgInfo::COMPUTE) {
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
      } else if (which[i] == ArgInfo::FIX) {
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
      } else if (which[i] == ArgInfo::VARIABLE) {
        // variable is always intensive, does not change extarray
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
  delete [] extlist;

  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::init()
{
  // set indices and check validity of all computes,fixes

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute slice does not exist");
      value2index[m] = icompute;
    } else if (which[m] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute slice does not exist");
      value2index[m] = ifix;
    } else if (which[m] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute slice does not exist");
      value2index[m] = ivariable;
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

  if (which[m] == ArgInfo::COMPUTE) {
    Compute *compute = modify->compute[value2index[m]];

    if (argindex[m] == 0) {
      if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
        compute->compute_vector();
        compute->invoked_flag |= Compute::INVOKED_VECTOR;
      }
      double *cvector = compute->vector;
      j = 0;
      for (i = nstart; i < nstop; i += nskip) {
        vec[j] = cvector[i-1];
        j += stride;
      }

    } else {
      if (!(compute->invoked_flag & Compute::INVOKED_ARRAY)) {
        compute->compute_array();
        compute->invoked_flag |= Compute::INVOKED_ARRAY;
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

  } else if (which[m] == ArgInfo::FIX) {
    if (update->ntimestep % modify->fix[value2index[m]]->global_freq)
      error->all(FLERR,"Fix used in compute slice not "
                 "computed at compatible time");
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

    // invoke vector-style variable

  } else if (which[m] == ArgInfo::VARIABLE) {
    double *varvec;
    int nvec = input->variable->compute_vector(value2index[m],&varvec);
    if (nvec < nstop)
      error->all(FLERR,"Compute slice variable is not long enough");
    j = 0;
    for (i = nstart; i < nstop; i += nskip) {
      vec[j] = varvec[i-1];
      j += stride;
    }
  }
}

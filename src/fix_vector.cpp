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

#include "fix_vector.h"

#include "arg_info.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ONE,RUNNING,WINDOW};
enum{SCALAR,VECTOR};


/* ---------------------------------------------------------------------- */

FixVector::FixVector(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0), which(nullptr), argindex(nullptr), value2index(nullptr), ids(nullptr), vector(nullptr), array(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix vector command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix vector command");

  nvalues = narg-4;
  which = new int[nvalues];
  argindex = new int[nvalues];
  value2index = new int[nvalues];
  ids = new char*[nvalues];

  nvalues = 0;
  for (int iarg = 4; iarg < narg; iarg++) {
    ArgInfo argi(arg[iarg]);

    which[nvalues] = argi.get_type();
    argindex[nvalues] = argi.get_index1();
    ids[nvalues] = argi.copy_name();

    if ((argi.get_type() == ArgInfo::UNKNOWN)
        || (argi.get_type() == ArgInfo::NONE)
        || (argi.get_dim() > 1))
      error->all(FLERR,"Illegal fix vector command");

    nvalues++;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix vector does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
        error->all(FLERR,"Fix vector compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,"Fix vector compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
        error->all(FLERR,
                   "Fix vector compute vector is accessed out-of-range");

    } else if (which[i] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix vector does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
        error->all(FLERR,"Fix vector fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,"Fix vector fix does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
        error->all(FLERR,"Fix vector fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix vector not computed at compatible time");

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix vector does not exist");
      if (argindex[i] == 0 && input->variable->equalstyle(ivariable) == 0)
        error->all(FLERR,"Fix vector variable is not equal-style variable");
      if (argindex[i] && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Fix vector variable is not vector-style variable");
    }
  }

  // this fix produces either a global vector or array
  // intensive/extensive flags set by compute,fix,variable that produces value

  int value,finalvalue;
  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      Compute *compute = modify->compute[modify->find_compute(ids[i])];
      if (argindex[0] == 0) value = compute->extscalar;
      else if (compute->extvector >= 0) value = compute->extvector;
      else value = compute->extlist[argindex[0]-1];
    } else if (which[i] == ArgInfo::FIX) {
      Fix *fix = modify->fix[modify->find_fix(ids[i])];
      if (argindex[i] == 0) value = fix->extvector;
      else value = fix->extarray;
    } else if (which[i] == ArgInfo::VARIABLE) value = 0;
    if (i == 0) finalvalue = value;
    else if (value != finalvalue)
      error->all(FLERR,"Fix vector cannot set output array "
                 "intensive/extensive from these inputs");
  }

  if (nvalues == 1) {
    vector_flag = 1;
    extvector = finalvalue;
  } else {
    array_flag = 1;
    size_array_cols = nvalues;
    extarray = finalvalue;
  }

  global_freq = nevery;
  time_depend = 1;

  // ncount = current size of vector or array

  vector = nullptr;
  array = nullptr;
  ncount = ncountmax = 0;
  if (nvalues == 1) size_vector = 0;
  else size_array_rows = 0;

  // nextstep = next step on which end_of_step does something
  // add nextstep to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nextstep = (update->ntimestep/nevery)*nevery;
  if (nextstep < update->ntimestep) nextstep += nevery;
  modify->addstep_compute_all(nextstep);

  // initialstep = first step the vector/array will store values for

  initialstep = nextstep;
}

/* ---------------------------------------------------------------------- */

FixVector::~FixVector()
{
  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixVector::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixVector::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix vector does not exist");
      value2index[i] = icompute;

    } else if (which[i] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix vector does not exist");
      value2index[i] = ifix;

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix vector does not exist");
      value2index[i] = ivariable;
    }
  }

  // reallocate vector or array for accumulated size at end of run
  // use endstep to allow for subsequent runs with "pre no"
  // nsize = # of entries from initialstep to finalstep

  bigint finalstep = update->endstep/nevery * nevery;
  if (finalstep > update->endstep) finalstep -= nevery;
  ncountmax = (finalstep-initialstep)/nevery + 1;
  if (nvalues == 1) memory->grow(vector,ncountmax,"vector:vector");
  else memory->grow(array,ncountmax,nvalues,"vector:array");
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixVector::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixVector::end_of_step()
{
  // skip if not step which requires doing something

  if (update->ntimestep != nextstep) return;
  if (ncount == ncountmax)
    error->all(FLERR,"Overflow of allocated fix vector storage");

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  double *result;
  if (nvalues == 1) result = &vector[ncount];
  else result = array[ncount];

  modify->clearstep_compute();

  for (int i = 0; i < nvalues; i++) {
    int m = value2index[i];

    // invoke compute if not previously invoked

    if (which[i] == ArgInfo::COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[i] == 0) {
        if (!(compute->invoked_flag & Compute::INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        result[i] = compute->scalar;
      } else {
        if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        result[i] = compute->vector[argindex[i]-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == ArgInfo::FIX) {
      if (argindex[i] == 0)
        result[i] = modify->fix[m]->compute_scalar();
      else
        result[i] = modify->fix[m]->compute_vector(argindex[i]-1);

    // evaluate equal-style or vector-style variable

    } else if (which[i] == ArgInfo::VARIABLE) {
      if (argindex[i] == 0)
        result[i] = input->variable->compute_equal(m);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(m,&varvec);
        int index = argindex[i];
        if (nvec < index) result[i] = 0.0;
        else result[i] = varvec[index-1];
      }
    }
  }

  // trigger computes on next needed step

  nextstep += nevery;
  modify->addstep_compute(nextstep);

  // update size of vector or array

  ncount++;
  if (nvalues == 1) size_vector++;
  else size_array_rows++;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixVector::compute_vector(int i)
{
  return vector[i];
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixVector::compute_array(int i, int j)
{
  return array[i][j];
}

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

/* ---------------------------------------------------------------------- */

FixVector::FixVector(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), vector(nullptr), array(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix vector", error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Invalid fix vector every argument: {}", nevery);

  nmaxval = MAXSMALLINT;
  nindex = 0;

  // parse values

  int iarg = 4;
  values.clear();

  while (iarg < narg) {
    ArgInfo argi(arg[iarg]);

    value_t val;
    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();
    val.val.c = nullptr;

    if ((val.which == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
      error->all(FLERR, "Invalid fix vector argument: {}", arg[iarg]);

    if (val.which == ArgInfo::NONE) break;

    values.push_back(val);
    ++iarg;
  }

  while (iarg < narg) {

    if (strcmp(arg[iarg], "nmax") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix vector nmax", error);
      nmaxval = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      if (nmaxval < 1) error->all(FLERR, "Invalid nmax value");
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix vector keyword: {}", arg[iarg]);
    }
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable
  // this fix produces either a global vector or array
  // intensive/extensive flags set by compute,fix,variable that produces value

  int value, finalvalue;
  bool first = true;
  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      auto icompute = modify->get_compute_by_id(val.id);
      if (!icompute) error->all(FLERR, "Compute ID {} for fix vector does not exist", val.id);
      if (val.argindex == 0 && icompute->scalar_flag == 0)
        error->all(FLERR, "Fix vector compute {} does not calculate a scalar", val.id);
      if (val.argindex && icompute->vector_flag == 0)
        error->all(FLERR, "Fix vector compute {} does not calculate a vector", val.id);
      if (val.argindex && (val.argindex > icompute->size_vector))
        error->all(FLERR, "Fix vector compute {} vector is accessed out-of-range", val.id);

      if (val.argindex == 0)
        value = icompute->extscalar;
      else if (icompute->extvector >= 0)
        value = icompute->extvector;
      else
        value = icompute->extlist[val.argindex - 1];
      val.val.c = icompute;

    } else if (val.which == ArgInfo::FIX) {
      auto ifix = modify->get_fix_by_id(val.id);
      if (!ifix) error->all(FLERR, "Fix ID {} for fix vector does not exist", val.id);
      if (val.argindex == 0 && ifix->scalar_flag == 0)
        error->all(FLERR, "Fix vector fix {} does not calculate a scalar", val.id);
      if (val.argindex && ifix->vector_flag == 0)
        error->all(FLERR, "Fix vector fix {} does not calculate a vector", val.id);
      if (val.argindex && val.argindex > ifix->size_vector)
        error->all(FLERR, "Fix vector fix {} vector is accessed out-of-range", val.id);
      if (nevery % ifix->global_freq)
        error->all(FLERR, "Fix for fix {} vector not computed at compatible time", val.id);

      if (val.argindex == 0)
        value = ifix->extvector;
      else
        value = ifix->extarray;
      val.val.f = ifix;

    } else if (val.which == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(val.id.c_str());
      if (ivariable < 0)
        error->all(FLERR, "Variable name {} for fix vector does not exist", val.id);
      if (val.argindex == 0 && input->variable->equalstyle(ivariable) == 0)
        error->all(FLERR, "Fix vector variable {} is not equal-style variable", val.id);
      if (val.argindex && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR, "Fix vector variable {} is not vector-style variable", val.id);
      value = 0;
      val.val.v = ivariable;
    }

    if (first) {
      finalvalue = value;
      first = false;
    } else if (value != finalvalue)
      error->all(FLERR, "Fix vector cannot set output array intensive/extensive from these inputs");
  }

  if (values.size() == 1) {
    vector_flag = 1;
    extvector = finalvalue;
  } else {
    array_flag = 1;
    size_array_cols = values.size();
    extarray = finalvalue;
  }

  global_freq = nevery;
  time_depend = 1;

  // ncount = current size of vector or array

  vector = nullptr;
  array = nullptr;
  ncount = ncountmax = nindex = 0;
  if (values.size() == 1)
    size_vector = 0;
  else
    size_array_rows = 0;

  // nextstep = next step on which end_of_step does something
  // add nextstep to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nextstep = (update->ntimestep / nevery) * nevery;
  if (nextstep < update->ntimestep) nextstep += nevery;
  modify->addstep_compute_all(nextstep);

  // initialstep = first step the vector/array will store values for

  initialstep = nextstep;
}

/* ---------------------------------------------------------------------- */

FixVector::~FixVector()
{
  values.clear();
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

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR, "Compute ID {} for fix vector does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for fix vector does not exist", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(val.id.c_str());
      if (ivariable < 0) error->all(FLERR, "Variable name for fix vector does not exist");
      val.val.v = ivariable;
    }
  }

  // reallocate vector or array for accumulated size at end of run
  // use endstep to allow for subsequent runs with "pre no"
  // nsize = # of entries from initialstep to finalstep

  bigint finalstep = update->endstep / nevery * nevery;
  if (finalstep > update->endstep) finalstep -= nevery;
  ncountmax = (finalstep - initialstep) / nevery + 1;
  if (ncountmax > nmaxval) ncountmax = nmaxval;
  if (values.size() == 1)
    memory->grow(vector, ncountmax, "vector:vector");
  else
    memory->grow(array, ncountmax, values.size(), "vector:array");
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

  // wrap around when vector/array is full
  nindex = ncount % ncountmax;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  double *result;
  if (values.size() == 1)
    result = &vector[nindex];
  else
    result = array[nindex];

  modify->clearstep_compute();

  int i = 0;
  for (auto &val : values) {

    // invoke compute if not previously invoked

    if (val.which == ArgInfo::COMPUTE) {

      if (val.argindex == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_SCALAR)) {
          val.val.c->compute_scalar();
          val.val.c->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        result[i] = val.val.c->scalar;
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        result[i] = val.val.c->vector[val.argindex - 1];
      }

      // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      if (val.argindex == 0)
        result[i] = val.val.f->compute_scalar();
      else
        result[i] = val.val.f->compute_vector(val.argindex - 1);

      // evaluate equal-style or vector-style variable

    } else if (val.which == ArgInfo::VARIABLE) {
      if (val.argindex == 0)
        result[i] = input->variable->compute_equal(val.val.v);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(val.val.v, &varvec);
        int index = val.argindex;
        if (nvec < index)
          result[i] = 0.0;
        else
          result[i] = varvec[index - 1];
      }
    }
    ++i;
  }

  // trigger computes on next needed step

  nextstep += nevery;
  modify->addstep_compute(nextstep);

  // update size of vector or array

  ncount++;
  if (values.size() == 1)
    size_vector = MIN(size_vector + 1, ncountmax);
  else
    size_array_rows = MIN(size_array_rows + 1, ncountmax);
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixVector::compute_vector(int i)
{
  int idx = i;
  if (ncount >= ncountmax) idx = (i + ncount) % ncountmax;
  return vector[idx];
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixVector::compute_array(int i, int j)
{
  int idx = i;
  if (ncount >= ncountmax) idx = (i + ncount) % ncountmax;
  return array[idx][j];
}

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

ComputeSlice::ComputeSlice(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "compute slice", error);

  nstart = utils::inumeric(FLERR, arg[3], false, lmp);
  nstop = utils::inumeric(FLERR, arg[4], false, lmp);
  nskip = utils::inumeric(FLERR, arg[5], false, lmp);

  if (nstart < 1) error->all(FLERR, "Invalid compute slice nstart value {} < 1", nstart);
  if (nstop < nstart) error->all(FLERR, "Invalid compute slice nstop value {} < {}", nstop, nstart);
  if (nskip < 1) error->all(FLERR, "Invalid compute slice nskip value < 1: {}", nskip);

  // parse values

  values.clear();
  for (int iarg = 6; iarg < narg; iarg++) {
    ArgInfo argi(arg[iarg]);

    value_t val;
    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();
    val.val.c = nullptr;

    if ((val.which == ArgInfo::UNKNOWN) || (val.which == ArgInfo::NONE) || (argi.get_dim() > 1))
      error->all(FLERR, "Illegal compute slice argument: {}", arg[iarg]);

    values.push_back(val);
  }

  // setup and error check

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR, "Compute ID {} for compute slice does not exist", val.id);
      if (val.val.c->vector_flag) {
        if (val.argindex)
          error->all(FLERR, "Compute slice compute {} does not calculate a global array", val.id);
        if (nstop > val.val.c->size_vector)
          error->all(FLERR, "Compute slice compute {} vector is accessed out-of-range", val.id);
      } else if (val.val.c->array_flag) {
        if (val.argindex == 0)
          error->all(FLERR, "Compute slice compute {} does not calculate a global vector", val.id);
        if (val.argindex > val.val.c->size_array_cols)
          error->all(FLERR, "Compute slice compute {} array is accessed out-of-range", val.id);
        if (nstop > val.val.c->size_array_rows)
          error->all(FLERR, "Compute slice compute {} array is accessed out-of-range", val.id);
      } else {
        error->all(FLERR, "Compute slice compute {} does not calculate global vector or array",
                   val.id);
      }
    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for compute slice does not exist", val.id);
      if (val.val.f->vector_flag) {
        if (val.argindex)
          error->all(FLERR, "Compute slice fix {} does not calculate a global array", val.id);
        if (nstop > val.val.f->size_vector)
          error->all(FLERR, "Compute slice fix {} vector is accessed out-of-range", val.id);
      } else if (val.val.f->array_flag) {
        if (val.argindex == 0)
          error->all(FLERR, "Compute slice fix {} does not calculate a global vector", val.id);
        if (val.argindex > val.val.f->size_array_cols)
          error->all(FLERR, "Compute slice fix {} array is accessed out-of-range", val.id);
        if (nstop > val.val.f->size_array_rows)
          error->all(FLERR, "Compute slice fix {} array is accessed out-of-range", val.id);
      } else {
        error->all(FLERR, "Compute slice fix {} does not calculate global vector or array", val.id);
      }
    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, "Variable name {} for compute slice does not exist", val.id);
      if (val.argindex == 0 && input->variable->vectorstyle(val.val.v) == 0)
        error->all(FLERR, "Compute slice variable {} is not vector-style variable", val.id);
      if (val.argindex)
        error->all(FLERR, "Compute slice vector variable {} cannot be indexed", val.id);
    }
  }

  // this compute produces either a vector or array
  // for vector, set intensive/extensive to mirror input values
  // for array, set intensive if all input values are intensive, else extensive

  if (values.size() == 1) {
    auto &val = values[0];
    vector_flag = 1;
    size_vector = (nstop - nstart) / nskip;
    memory->create(vector, size_vector, "slice:vector");

    if (val.which == ArgInfo::COMPUTE) {
      if (val.argindex == 0) {
        extvector = val.val.c->extvector;
        if (val.val.c->extvector == -1) {
          extlist = new int[size_vector];
          int j = 0;
          for (int i = nstart; i < nstop; i += nskip) extlist[j++] = val.val.c->extlist[i - 1];
        }
      } else
        extvector = val.val.c->extarray;
    } else if (val.which == ArgInfo::FIX) {
      if (val.argindex == 0) {
        extvector = val.val.f->extvector;
        if (val.val.f->extvector == -1) {
          extlist = new int[size_vector];
          int j = 0;
          for (int i = nstart; i < nstop; i += nskip) extlist[j++] = val.val.f->extlist[i - 1];
        }
      } else
        extvector = val.val.f->extarray;
    } else if (val.which == ArgInfo::VARIABLE) {
      extvector = 0;
    }

  } else {
    array_flag = 1;
    size_array_rows = (nstop - nstart) / nskip;
    size_array_cols = values.size();
    memory->create(array, size_array_rows, size_array_cols, "slice:array");

    extarray = 0;
    for (auto &val : values) {
      if (val.which == ArgInfo::COMPUTE) {
        if (val.argindex == 0) {
          if (val.val.c->extvector == 1) extarray = 1;
          if (val.val.c->extvector == -1) {
            for (int j = 0; j < val.val.c->size_vector; j++)
              if (val.val.c->extlist[j]) extarray = 1;
          }
        } else {
          if (val.val.c->extarray) extarray = 1;
        }
      } else if (val.which == ArgInfo::FIX) {
        if (val.argindex == 0) {
          if (val.val.f->extvector == 1) extarray = 1;
          if (val.val.f->extvector == -1) {
            for (int j = 0; j < val.val.f->size_vector; j++)
              if (val.val.f->extlist[j]) extarray = 1;
          }
        } else {
          if (val.val.f->extarray) extarray = 1;
        }
      } else if (val.which == ArgInfo::VARIABLE) {
        // variable is always intensive, does not change extarray
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

ComputeSlice::~ComputeSlice()
{
  delete[] extlist;

  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::init()
{
  // set indices and check validity of all computes,fixes

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR, "Compute ID {} for compute slice does not exist", val.id);
    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for compute slice does not exist", val.id);
    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, "Variable name {} for compute slice does not exist", val.id);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::compute_vector()
{
  invoked_vector = update->ntimestep;

  extract_one(0, vector, 1);
}

/* ---------------------------------------------------------------------- */

void ComputeSlice::compute_array()
{
  invoked_array = update->ntimestep;

  for (std::size_t m = 0; m < values.size(); m++) extract_one(0, &array[m][0], values.size());
}

/* ----------------------------------------------------------------------
   calculate sliced value for one input M and return it in vec
   vec may be array so that returned values are with stride
------------------------------------------------------------------------- */

void ComputeSlice::extract_one(int m, double *vec, int stride)
{
  auto &val = values[m];

  // invoke the appropriate compute if needed

  if (val.which == ArgInfo::COMPUTE) {
    if (val.argindex == 0) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
        val.val.c->compute_vector();
        val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
      }
      double *cvector = val.val.c->vector;
      int j = 0;
      for (int i = nstart; i < nstop; i += nskip) {
        vec[j] = cvector[i - 1];
        j += stride;
      }

    } else {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_ARRAY)) {
        val.val.c->compute_array();
        val.val.c->invoked_flag |= Compute::INVOKED_ARRAY;
      }
      double **carray = val.val.c->array;
      int icol = val.argindex - 1;
      int j = 0;
      for (int i = nstart; i < nstop; i += nskip) {
        vec[j] = carray[i - 1][icol];
        j += stride;
      }
    }

    // access fix fields, check if fix frequency is a match

  } else if (val.which == ArgInfo::FIX) {
    if (update->ntimestep % val.val.f->global_freq)
      error->all(FLERR, "Fix {} used in compute slice not computed at compatible time", val.id);

    if (val.argindex == 0) {
      int j = 0;
      for (int i = nstart; i < nstop; i += nskip) {
        vec[j] = val.val.f->compute_vector(i - 1);
        j += stride;
      }
    } else {
      int icol = val.argindex - 1;
      int j = 0;
      for (int i = nstart; i < nstop; i += nskip) {
        vec[j] = val.val.f->compute_array(i - 1, icol);
        j += stride;
      }
    }

    // invoke vector-style variable

  } else if (val.which == ArgInfo::VARIABLE) {
    double *varvec;
    int nvec = input->variable->compute_vector(val.val.v, &varvec);
    if (nvec < nstop) error->all(FLERR, "Compute slice variable {} is not long enough", val.id);
    int j = 0;
    for (int i = nstart; i < nstop; i += nskip) {
      vec[j] = varvec[i - 1];
      j += stride;
    }
  }
}

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

#include "compute_reduce_region.h"

#include "arg_info.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

ComputeReduceRegion::ComputeReduceRegion(LAMMPS *lmp, int narg, char **arg) :
    ComputeReduce(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   calculate reduced value for one input M and return it
   if flag = -1:
     sum/min/max/ave all values in vector
     for per-atom quantities, limit to atoms in group and region
     if mode = MIN or MAX, also set index to which vector value wins
   if flag >= 0: simply return vector[flag]
------------------------------------------------------------------------- */

double ComputeReduceRegion::compute_one(int m, int flag)
{
  region->prematch();

  // invoke the appropriate attribute,compute,fix,variable
  // compute scalar quantity by summing over atom scalars
  // only include atoms in group

  index = -1;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = value2index[m];

  // initialization in case it has not yet been run,
  // e.g. when invoked
  if (n == ArgInfo::UNKNOWN) {
    init();
    n = value2index[m];
  }

  int j = argindex[m];

  double one = 0.0;
  if (mode == MINN) one = BIG;
  if (mode == MAXX) one = -BIG;

  if (which[m] == ArgInfo::X) {
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, x[i][j], i);
    } else
      one = x[flag][j];
  } else if (which[m] == ArgInfo::V) {
    double **v = atom->v;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, v[i][j], i);
    } else
      one = v[flag][j];
  } else if (which[m] == ArgInfo::F) {
    double **f = atom->f;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, f[i][j], i);
    } else
      one = f[flag][j];

    // invoke compute if not previously invoked

  } else if (which[m] == ArgInfo::COMPUTE) {
    Compute *compute = modify->compute[n];

    if (flavor[m] == PERATOM) {
      if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= Compute::INVOKED_PERATOM;
      }

      if (j == 0) {
        double *compute_vector = compute->vector_atom;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, compute_vector[i], i);
        } else
          one = compute_vector[flag];
      } else {
        double **compute_array = compute->array_atom;
        int jm1 = j - 1;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, compute_array[i][jm1], i);
        } else
          one = compute_array[flag][jm1];
      }

    } else if (flavor[m] == LOCAL) {
      if (!(compute->invoked_flag & Compute::INVOKED_LOCAL)) {
        compute->compute_local();
        compute->invoked_flag |= Compute::INVOKED_LOCAL;
      }

      if (j == 0) {
        double *compute_vector = compute->vector_local;
        if (flag < 0)
          for (int i = 0; i < compute->size_local_rows; i++) combine(one, compute_vector[i], i);
        else
          one = compute_vector[flag];
      } else {
        double **compute_array = compute->array_local;
        int jm1 = j - 1;
        if (flag < 0)
          for (int i = 0; i < compute->size_local_rows; i++) combine(one, compute_array[i][jm1], i);
        else
          one = compute_array[flag][jm1];
      }
    }

    // check if fix frequency is a match

  } else if (which[m] == ArgInfo::FIX) {
    if (update->ntimestep % modify->fix[n]->peratom_freq)
      error->all(FLERR, "Fix used in compute reduce not computed at compatible time");
    Fix *fix = modify->fix[n];

    if (flavor[m] == PERATOM) {
      if (j == 0) {
        double *fix_vector = fix->vector_atom;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, fix_vector[i], i);
        } else
          one = fix_vector[flag];
      } else {
        double **fix_array = fix->array_atom;
        int jm1 = j - 1;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, fix_array[i][jm1], i);
        } else
          one = fix_array[flag][jm1];
      }

    } else if (flavor[m] == LOCAL) {
      if (j == 0) {
        double *fix_vector = fix->vector_local;
        if (flag < 0)
          for (int i = 0; i < fix->size_local_rows; i++) combine(one, fix_vector[i], i);
        else
          one = fix_vector[flag];
      } else {
        double **fix_array = fix->array_local;
        int jm1 = j - 1;
        if (flag < 0)
          for (int i = 0; i < fix->size_local_rows; i++) combine(one, fix_array[i][jm1], i);
        else
          one = fix_array[flag][jm1];
      }
    }

    // evaluate atom-style variable

  } else if (which[m] == ArgInfo::VARIABLE) {
    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom, maxatom, "reduce/region:varatom");
    }

    input->variable->compute_atom(n, igroup, varatom, 1, 0);
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, varatom[i], i);
    } else
      one = varatom[flag];
  }

  return one;
}

/* ---------------------------------------------------------------------- */

bigint ComputeReduceRegion::count(int m)
{
  int n = value2index[m];

  if (which[m] == ArgInfo::X || which[m] == ArgInfo::V || which[m] == ArgInfo::F)
    return group->count(igroup, region);
  else if (which[m] == ArgInfo::COMPUTE) {
    Compute *compute = modify->compute[n];
    if (flavor[m] == PERATOM) {
      return group->count(igroup, region);
    } else if (flavor[m] == LOCAL) {
      bigint ncount = compute->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_DOUBLE, MPI_SUM, world);
      return ncountall;
    }
  } else if (which[m] == ArgInfo::FIX) {
    Fix *fix = modify->fix[n];
    if (flavor[m] == PERATOM) {
      return group->count(igroup, region);
    } else if (flavor[m] == LOCAL) {
      bigint ncount = fix->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_DOUBLE, MPI_SUM, world);
      return ncountall;
    }
  } else if (which[m] == ArgInfo::VARIABLE)
    return group->count(igroup, region);

  bigint dummy = 0;
  return dummy;
}

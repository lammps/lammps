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

#include "compute_reduce_region.h"

#include "arg_info.h"
#include "atom.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "memory.h"
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
  auto &val = values[m];

  // initialization in case it has not yet been run, e.g. when
  // the compute was invoked right after it has been created
  if ((val.which == ArgInfo::COMPUTE) || (val.which == ArgInfo::FIX)) {
    if (val.val.c == nullptr) init();
  }

  int aidx = val.argindex;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double one = 0.0;
  if (mode == MINN) one = BIG;
  if (mode == MAXX) one = -BIG;

  if (val.which == ArgInfo::X) {
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, x[i][aidx], i);
    } else
      one = x[flag][aidx];
  } else if (val.which == ArgInfo::V) {
    double **v = atom->v;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, v[i][aidx], i);
    } else
      one = v[flag][aidx];
  } else if (val.which == ArgInfo::F) {
    double **f = atom->f;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
          combine(one, f[i][aidx], i);
    } else
      one = f[flag][aidx];

    // invoke compute if not previously invoked

  } else if (val.which == ArgInfo::COMPUTE) {
    if (val.flavor == PERATOM) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
        val.val.c->compute_peratom();
        val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
      }

      if (aidx == 0) {
        double *compute_vector = val.val.c->vector_atom;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, compute_vector[i], i);
        } else
          one = compute_vector[flag];
      } else {
        double **compute_array = val.val.c->array_atom;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, compute_array[i][aidxm1], i);
        } else
          one = compute_array[flag][aidxm1];
      }

    } else if (val.flavor == LOCAL) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_LOCAL)) {
        val.val.c->compute_local();
        val.val.c->invoked_flag |= Compute::INVOKED_LOCAL;
      }

      if (aidx == 0) {
        double *compute_vector = val.val.c->vector_local;
        if (flag < 0)
          for (int i = 0; i < val.val.c->size_local_rows; i++) combine(one, compute_vector[i], i);
        else
          one = compute_vector[flag];
      } else {
        double **compute_array = val.val.c->array_local;
        int aidxm1 = aidx - 1;
        if (flag < 0)
          for (int i = 0; i < val.val.c->size_local_rows; i++)
            combine(one, compute_array[i][aidxm1], i);
        else
          one = compute_array[flag][aidxm1];
      }
    }

    // check if fix frequency is a match

  } else if (val.which == ArgInfo::FIX) {
    if (update->ntimestep % val.val.f->peratom_freq)
      error->all(FLERR, "Fix {} used in compute {} not computed at compatible time", val.id, style);

    if (val.flavor == PERATOM) {
      if (aidx == 0) {
        double *fix_vector = val.val.f->vector_atom;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, fix_vector[i], i);
        } else
          one = fix_vector[flag];
      } else {
        double **fix_array = val.val.f->array_atom;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
              combine(one, fix_array[i][aidxm1], i);
        } else
          one = fix_array[flag][aidxm1];
      }

    } else if (val.flavor == LOCAL) {
      if (aidx == 0) {
        double *fix_vector = val.val.f->vector_local;
        if (flag < 0)
          for (int i = 0; i < val.val.f->size_local_rows; i++) combine(one, fix_vector[i], i);
        else
          one = fix_vector[flag];
      } else {
        double **fix_array = val.val.f->array_local;
        int aidxm1 = aidx - 1;
        if (flag < 0)
          for (int i = 0; i < val.val.f->size_local_rows; i++)
            combine(one, fix_array[i][aidxm1], i);
        else
          one = fix_array[flag][aidxm1];
      }
    }

    // evaluate atom-style variable

  } else if (val.which == ArgInfo::VARIABLE) {
    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom, maxatom, "reduce/region:varatom");
    }

    input->variable->compute_atom(val.val.v, igroup, varatom, 1, 0);
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
  auto &val = values[m];

  if (val.which == ArgInfo::X || val.which == ArgInfo::V || val.which == ArgInfo::F)
    return group->count(igroup, region);
  else if (val.which == ArgInfo::COMPUTE) {
    if (val.flavor == PERATOM) {
      return group->count(igroup, region);
    } else if (val.flavor == LOCAL) {
      bigint ncount = val.val.c->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_DOUBLE, MPI_SUM, world);
      return ncountall;
    }
  } else if (val.which == ArgInfo::FIX) {
    if (val.flavor == PERATOM) {
      return group->count(igroup, region);
    } else if (val.flavor == LOCAL) {
      bigint ncount = val.val.f->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_DOUBLE, MPI_SUM, world);
      return ncountall;
    }
  } else if (val.which == ArgInfo::VARIABLE)
    return group->count(igroup, region);

  bigint dummy = 0;
  return dummy;
}

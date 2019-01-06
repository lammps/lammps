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

#include <cstring>
#include <cstdlib>
#include "compute_reduce_region.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "group.h"
#include "region.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SUM,SUMSQ,MINN,MAXX,AVE,AVESQ};             // also in ComputeReduce
enum{UNKNOWN=-1,X,V,F,COMPUTE,FIX,VARIABLE};
enum{PERATOM,LOCAL};

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduceRegion::ComputeReduceRegion(LAMMPS *lmp, int narg, char **arg) :
  ComputeReduce(lmp, narg, arg) {}

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
  int i;

  Region *region = domain->regions[iregion];
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
  if (n == UNKNOWN) {
    init();
    n = value2index[m];
  }

  int aidx = argindex[m];
  int j = argindex[m];

  double one = 0.0;
  if (mode == MINN) one = BIG;
  if (mode == MAXX) one = -BIG;

  if (which[m] == X) {
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
          combine(one,x[i][j],i);
    } else one = x[flag][j];
  } else if (which[m] == V) {
    double **v = atom->v;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
          combine(one,v[i][j],i);
    } else one = v[flag][j];
  } else if (which[m] == F) {
    double **f = atom->f;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
          combine(one,f[i][j],i);
    } else one = f[flag][j];

  // invoke compute if not previously invoked

  } else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[n];

    if (flavor[m] == PERATOM) {
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
      }

      if (j == 0) {
        double *compute_vector = compute->vector_atom;
        int n = nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
              combine(one,compute_vector[i],i);
        } else one = compute_vector[flag];
      } else {
        double **compute_array = compute->array_atom;
        int n = nlocal;
        int jm1 = j - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
              combine(one,compute_array[i][jm1],i);
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

  // check if fix frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[n]->peratom_freq)
      error->all(FLERR,"Fix used in compute reduce not computed at "
                 "compatible time");
    Fix *fix = modify->fix[n];

    if (flavor[m] == PERATOM) {
      if (j == 0) {
        double *fix_vector = fix->vector_atom;
        int n = nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
              combine(one,fix_vector[i],i);
        } else one = fix_vector[flag];
      } else {
        double **fix_array = fix->array_atom;
        int jm1 = j - 1;
        if (flag < 0) {
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
              combine(one,fix_array[i][jm1],i);
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
    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom,maxatom,"reduce/region:varatom");
    }

    input->variable->compute_atom(n,igroup,varatom,1,0);
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
          combine(one,varatom[i],i);
    } else one = varatom[flag];
  }

  return one;
}

/* ---------------------------------------------------------------------- */

bigint ComputeReduceRegion::count(int m)
{
  int n = value2index[m];

  if (which[m] == X || which[m] == V || which[m] == F)
    return group->count(igroup,iregion);
  else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[n];
    if (flavor[m] == PERATOM) {
      return group->count(igroup,iregion);
    } else if (flavor[m] == LOCAL) {
      bigint ncount = compute->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_DOUBLE,MPI_SUM,world);
      return ncountall;
    }
  } else if (which[m] == FIX) {
    Fix *fix = modify->fix[n];
    if (flavor[m] == PERATOM) {
      return group->count(igroup,iregion);
    } else if (flavor[m] == LOCAL) {
      bigint ncount = fix->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_DOUBLE,MPI_SUM,world);
      return ncountall;
    }
  } else if (which[m] == VARIABLE)
    return group->count(igroup,iregion);

  bigint dummy = 0;
  return dummy;
}

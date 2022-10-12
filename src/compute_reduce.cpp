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

#include "compute_reduce.h"

#include "arg_info.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduce::ComputeReduce(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), nvalues(0), which(nullptr), argindex(nullptr), flavor(nullptr),
    value2index(nullptr), ids(nullptr), onevec(nullptr), replace(nullptr), indices(nullptr),
    owner(nullptr), idregion(nullptr), region(nullptr), varatom(nullptr)
{
  int iarg = 0;
  if (strcmp(style, "reduce") == 0) {
    if (narg < 5) error->all(FLERR, "Illegal compute reduce command");
    iarg = 3;
  } else if (strcmp(style, "reduce/region") == 0) {
    if (narg < 6) error->all(FLERR, "Illegal compute reduce/region command");
    if (!domain->get_region_by_id(arg[3]))
      error->all(FLERR, "Region {} for compute reduce/region does not exist", arg[3]);
    idregion = utils::strdup(arg[3]);
    iarg = 4;
  }

  if (strcmp(arg[iarg], "sum") == 0)
    mode = SUM;
  else if (strcmp(arg[iarg], "sumsq") == 0)
    mode = SUMSQ;
  else if (strcmp(arg[iarg], "sumabs") == 0)
    mode = SUMABS;
  else if (strcmp(arg[iarg], "min") == 0)
    mode = MINN;
  else if (strcmp(arg[iarg], "max") == 0)
    mode = MAXX;
  else if (strcmp(arg[iarg], "ave") == 0)
    mode = AVE;
  else if (strcmp(arg[iarg], "avesq") == 0)
    mode = AVESQ;
  else if (strcmp(arg[iarg], "aveabs") == 0)
    mode = AVEABS;
  else
    error->all(FLERR, "Illegal compute {} operation {}", style, arg[iarg]);
  iarg++;

  MPI_Comm_rank(world, &me);

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR, narg - iarg, &arg[iarg], 1, earg, lmp);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  which = new int[nargnew];
  argindex = new int[nargnew];
  flavor = new int[nargnew];
  ids = new char *[nargnew];
  value2index = new int[nargnew];
  for (int i = 0; i < nargnew; ++i) {
    which[i] = argindex[i] = flavor[i] = value2index[i] = ArgInfo::UNKNOWN;
    ids[i] = nullptr;
  }
  nvalues = 0;

  iarg = 0;
  while (iarg < nargnew) {
    ids[nvalues] = nullptr;

    if (strcmp(arg[iarg], "x") == 0) {
      which[nvalues] = ArgInfo::X;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg], "y") == 0) {
      which[nvalues] = ArgInfo::X;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg], "z") == 0) {
      which[nvalues] = ArgInfo::X;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg], "vx") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg], "vy") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg], "vz") == 0) {
      which[nvalues] = ArgInfo::V;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg], "fx") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg], "fy") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg], "fz") == 0) {
      which[nvalues] = ArgInfo::F;
      argindex[nvalues++] = 2;

    } else {

      ArgInfo argi(arg[iarg]);

      which[nvalues] = argi.get_type();
      argindex[nvalues] = argi.get_index1();
      ids[nvalues] = argi.copy_name();

      if ((which[nvalues] == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR, "Illegal compute reduce command");

      if (which[nvalues] == ArgInfo::NONE) break;
      nvalues++;
    }

    iarg++;
  }

  // optional args

  replace = new int[nvalues];
  for (int i = 0; i < nvalues; i++) replace[i] = -1;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg], "replace") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal compute reduce command");
      if (mode != MINN && mode != MAXX)
        error->all(FLERR, "Compute reduce replace requires min or max mode");
      int col1 = utils::inumeric(FLERR, arg[iarg + 1], false, lmp) - 1;
      int col2 = utils::inumeric(FLERR, arg[iarg + 2], false, lmp) - 1;
      if (col1 < 0 || col1 >= nvalues || col2 < 0 || col2 >= nvalues)
        error->all(FLERR, "Illegal compute reduce command");
      if (col1 == col2) error->all(FLERR, "Illegal compute reduce command");
      if (replace[col1] >= 0 || replace[col2] >= 0)
        error->all(FLERR, "Invalid replace values in compute reduce");
      replace[col1] = col2;
      iarg += 3;
    } else
      error->all(FLERR, "Illegal compute reduce command");
  }

  // delete replace if not set

  int flag = 0;
  for (int i = 0; i < nvalues; i++)
    if (replace[i] >= 0) flag = 1;
  if (!flag) {
    delete[] replace;
    replace = nullptr;
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::X || which[i] == ArgInfo::V || which[i] == ArgInfo::F)
      flavor[i] = PERATOM;

    else if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0) error->all(FLERR, "Compute ID for compute reduce does not exist");
      if (modify->compute[icompute]->peratom_flag) {
        flavor[i] = PERATOM;
        if (argindex[i] == 0 && modify->compute[icompute]->size_peratom_cols != 0)
          error->all(FLERR,
                     "Compute reduce compute does not "
                     "calculate a per-atom vector");
        if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
          error->all(FLERR,
                     "Compute reduce compute does not "
                     "calculate a per-atom array");
        if (argindex[i] && argindex[i] > modify->compute[icompute]->size_peratom_cols)
          error->all(FLERR, "Compute reduce compute array is accessed out-of-range");
      } else if (modify->compute[icompute]->local_flag) {
        flavor[i] = LOCAL;
        if (argindex[i] == 0 && modify->compute[icompute]->size_local_cols != 0)
          error->all(FLERR,
                     "Compute reduce compute does not "
                     "calculate a local vector");
        if (argindex[i] && modify->compute[icompute]->size_local_cols == 0)
          error->all(FLERR,
                     "Compute reduce compute does not "
                     "calculate a local array");
        if (argindex[i] && argindex[i] > modify->compute[icompute]->size_local_cols)
          error->all(FLERR, "Compute reduce compute array is accessed out-of-range");
      } else
        error->all(FLERR, "Compute reduce compute calculates global values");

    } else if (which[i] == ArgInfo::FIX) {
      auto ifix = modify->get_fix_by_id(ids[i]);
      if (!ifix) error->all(FLERR, "Fix ID {} for compute reduce does not exist", ids[i]);
      if (ifix->peratom_flag) {
        flavor[i] = PERATOM;
        if (argindex[i] == 0 && (ifix->size_peratom_cols != 0))
          error->all(FLERR, "Compute reduce fix {} does not calculate a per-atom vector", ids[i]);
        if (argindex[i] && (ifix->size_peratom_cols == 0))
          error->all(FLERR, "Compute reduce fix {} does not calculate a per-atom array", ids[i]);
        if (argindex[i] && (argindex[i] > ifix->size_peratom_cols))
          error->all(FLERR, "Compute reduce fix {} array is accessed out-of-range", ids[i]);
      } else if (ifix->local_flag) {
        flavor[i] = LOCAL;
        if (argindex[i] == 0 && (ifix->size_local_cols != 0))
          error->all(FLERR, "Compute reduce fix {} does not calculate a local vector", ids[i]);
        if (argindex[i] && (ifix->size_local_cols == 0))
          error->all(FLERR, "Compute reduce fix {} does not calculate a local array", ids[i]);
        if (argindex[i] && (argindex[i] > ifix->size_local_cols))
          error->all(FLERR, "Compute reduce fix {} array is accessed out-of-range", ids[i]);
      } else
        error->all(FLERR, "Compute reduce fix {} calculates global values", ids[i]);

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0) error->all(FLERR, "Variable name for compute reduce does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR, "Compute reduce variable is not atom-style variable");
      flavor[i] = PERATOM;
    }
  }

  // this compute produces either a scalar or vector

  if (nvalues == 1) {
    scalar_flag = 1;
    if (mode == SUM || mode == SUMSQ || mode == SUMABS)
      extscalar = 1;
    else
      extscalar = 0;
    vector = onevec = nullptr;
    indices = owner = nullptr;
  } else {
    vector_flag = 1;
    size_vector = nvalues;
    if (mode == SUM || mode == SUMSQ || mode == SUMABS)
      extvector = 1;
    else
      extvector = 0;
    vector = new double[size_vector];
    onevec = new double[size_vector];
    indices = new int[size_vector];
    owner = new int[size_vector];
  }

  maxatom = 0;
  varatom = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeReduce::~ComputeReduce()
{
  delete[] which;
  delete[] argindex;
  delete[] flavor;
  for (int m = 0; m < nvalues; m++) delete[] ids[m];
  delete[] ids;
  delete[] value2index;
  delete[] replace;
  delete[] idregion;

  delete[] vector;
  delete[] onevec;
  delete[] indices;
  delete[] owner;

  memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::init()
{
  // set indices of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0) error->all(FLERR, "Compute ID for compute reduce does not exist");
      value2index[m] = icompute;

    } else if (which[m] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0) error->all(FLERR, "Fix ID for compute reduce does not exist");
      value2index[m] = ifix;

    } else if (which[m] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0) error->all(FLERR, "Variable name for compute reduce does not exist");
      value2index[m] = ivariable;

    } else
      value2index[m] = ArgInfo::UNKNOWN;
  }

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for compute reduce/region does not exist", idregion);
  }
}

/* ---------------------------------------------------------------------- */

double ComputeReduce::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double one = compute_one(0, -1);

  if (mode == SUM || mode == SUMSQ || mode == SUMABS) {
    MPI_Allreduce(&one, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  } else if (mode == MINN) {
    MPI_Allreduce(&one, &scalar, 1, MPI_DOUBLE, MPI_MIN, world);
  } else if (mode == MAXX) {
    MPI_Allreduce(&one, &scalar, 1, MPI_DOUBLE, MPI_MAX, world);
  } else if (mode == AVE || mode == AVESQ || mode == AVEABS) {
    MPI_Allreduce(&one, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
    bigint n = count(0);
    if (n) scalar /= n;
  }

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::compute_vector()
{
  invoked_vector = update->ntimestep;

  for (int m = 0; m < nvalues; m++)
    if (!replace || replace[m] < 0) {
      onevec[m] = compute_one(m, -1);
      indices[m] = index;
    }

  if (mode == SUM || mode == SUMSQ || mode == AVEABS) {
    for (int m = 0; m < nvalues; m++)
      MPI_Allreduce(&onevec[m], &vector[m], 1, MPI_DOUBLE, MPI_SUM, world);

  } else if (mode == MINN) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++)
        MPI_Allreduce(&onevec[m], &vector[m], 1, MPI_DOUBLE, MPI_MIN, world);

    } else {
      for (int m = 0; m < nvalues; m++)
        if (replace[m] < 0) {
          pairme.value = onevec[m];
          pairme.proc = me;
          MPI_Allreduce(&pairme, &pairall, 1, MPI_DOUBLE_INT, MPI_MINLOC, world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (me == owner[replace[m]]) vector[m] = compute_one(m, indices[replace[m]]);
          MPI_Bcast(&vector[m], 1, MPI_DOUBLE, owner[replace[m]], world);
        }
    }

  } else if (mode == MAXX) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++)
        MPI_Allreduce(&onevec[m], &vector[m], 1, MPI_DOUBLE, MPI_MAX, world);

    } else {
      for (int m = 0; m < nvalues; m++)
        if (replace[m] < 0) {
          pairme.value = onevec[m];
          pairme.proc = me;
          MPI_Allreduce(&pairme, &pairall, 1, MPI_DOUBLE_INT, MPI_MAXLOC, world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (me == owner[replace[m]]) vector[m] = compute_one(m, indices[replace[m]]);
          MPI_Bcast(&vector[m], 1, MPI_DOUBLE, owner[replace[m]], world);
        }
    }

  } else if (mode == AVE || mode == AVESQ || mode == AVEABS) {
    for (int m = 0; m < nvalues; m++) {
      MPI_Allreduce(&onevec[m], &vector[m], 1, MPI_DOUBLE, MPI_SUM, world);
      bigint n = count(m);
      if (n) vector[m] /= n;
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
  int vidx = value2index[m];

  // initialization in case it has not yet been run, e.g. when
  // the compute was invoked right after it has been created

  if (vidx == ArgInfo::UNKNOWN) {
    init();
    vidx = value2index[m];
  }

  int aidx = argindex[m];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double one = 0.0;
  if (mode == MINN) one = BIG;
  if (mode == MAXX) one = -BIG;

  if (which[m] == ArgInfo::X) {
    double **x = atom->x;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, x[i][aidx], i);
    } else
      one = x[flag][aidx];
  } else if (which[m] == ArgInfo::V) {
    double **v = atom->v;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, v[i][aidx], i);
    } else
      one = v[flag][aidx];
  } else if (which[m] == ArgInfo::F) {
    double **f = atom->f;
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, f[i][aidx], i);
    } else
      one = f[flag][aidx];

    // invoke compute if not previously invoked

  } else if (which[m] == ArgInfo::COMPUTE) {
    Compute *compute = modify->compute[vidx];

    if (flavor[m] == PERATOM) {
      if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= Compute::INVOKED_PERATOM;
      }

      if (aidx == 0) {
        double *comp_vec = compute->vector_atom;
        int n = nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            if (mask[i] & groupbit) combine(one, comp_vec[i], i);
        } else
          one = comp_vec[flag];
      } else {
        double **carray_atom = compute->array_atom;
        int n = nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            if (mask[i] & groupbit) combine(one, carray_atom[i][aidxm1], i);
        } else
          one = carray_atom[flag][aidxm1];
      }

    } else if (flavor[m] == LOCAL) {
      if (!(compute->invoked_flag & Compute::INVOKED_LOCAL)) {
        compute->compute_local();
        compute->invoked_flag |= Compute::INVOKED_LOCAL;
      }

      if (aidx == 0) {
        double *comp_vec = compute->vector_local;
        int n = compute->size_local_rows;
        if (flag < 0)
          for (i = 0; i < n; i++) combine(one, comp_vec[i], i);
        else
          one = comp_vec[flag];
      } else {
        double **carray_local = compute->array_local;
        int n = compute->size_local_rows;
        int aidxm1 = aidx - 1;
        if (flag < 0)
          for (i = 0; i < n; i++) combine(one, carray_local[i][aidxm1], i);
        else
          one = carray_local[flag][aidxm1];
      }
    }

    // access fix fields, check if fix frequency is a match

  } else if (which[m] == ArgInfo::FIX) {
    if (update->ntimestep % modify->fix[vidx]->peratom_freq)
      error->all(FLERR,
                 "Fix used in compute reduce not "
                 "computed at compatible time");
    Fix *fix = modify->fix[vidx];

    if (flavor[m] == PERATOM) {
      if (aidx == 0) {
        double *fix_vector = fix->vector_atom;
        int n = nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            if (mask[i] & groupbit) combine(one, fix_vector[i], i);
        } else
          one = fix_vector[flag];
      } else {
        double **fix_array = fix->array_atom;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) combine(one, fix_array[i][aidxm1], i);
        } else
          one = fix_array[flag][aidxm1];
      }

    } else if (flavor[m] == LOCAL) {
      if (aidx == 0) {
        double *fix_vector = fix->vector_local;
        int n = fix->size_local_rows;
        if (flag < 0)
          for (i = 0; i < n; i++) combine(one, fix_vector[i], i);
        else
          one = fix_vector[flag];
      } else {
        double **fix_array = fix->array_local;
        int n = fix->size_local_rows;
        int aidxm1 = aidx - 1;
        if (flag < 0)
          for (i = 0; i < n; i++) combine(one, fix_array[i][aidxm1], i);
        else
          one = fix_array[flag][aidxm1];
      }
    }

    // evaluate atom-style variable

  } else if (which[m] == ArgInfo::VARIABLE) {
    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom, maxatom, "reduce:varatom");
    }

    input->variable->compute_atom(vidx, igroup, varatom, 1, 0);
    if (flag < 0) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, varatom[i], i);
    } else
      one = varatom[flag];
  }

  return one;
}

/* ---------------------------------------------------------------------- */

bigint ComputeReduce::count(int m)
{
  int vidx = value2index[m];

  if (which[m] == ArgInfo::X || which[m] == ArgInfo::V || which[m] == ArgInfo::F)
    return group->count(igroup);
  else if (which[m] == ArgInfo::COMPUTE) {
    Compute *compute = modify->compute[vidx];
    if (flavor[m] == PERATOM) {
      return group->count(igroup);
    } else if (flavor[m] == LOCAL) {
      bigint ncount = compute->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
      return ncountall;
    }
  } else if (which[m] == ArgInfo::FIX) {
    Fix *fix = modify->fix[vidx];
    if (flavor[m] == PERATOM) {
      return group->count(igroup);
    } else if (flavor[m] == LOCAL) {
      bigint ncount = fix->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
      return ncountall;
    }
  } else if (which[m] == ArgInfo::VARIABLE)
    return group->count(igroup);

  bigint dummy = 0;
  return dummy;
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
   for MIN/MAX, also update index with winner
------------------------------------------------------------------------- */

void ComputeReduce::combine(double &one, double two, int i)
{
  if (mode == SUM || mode == AVE)
    one += two;
  else if (mode == SUMSQ || mode == AVESQ)
    one += two * two;
  else if (mode == SUMABS || mode == AVEABS)
    one += std::fabs(two);
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
  double bytes = (double) maxatom * sizeof(double);
  return bytes;
}

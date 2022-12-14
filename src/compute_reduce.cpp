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

#include "compute_reduce.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduce::ComputeReduce(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), nvalues(0), onevec(nullptr), replace(nullptr), indices(nullptr),
    owner(nullptr), idregion(nullptr), region(nullptr), varatom(nullptr)
{
  int iarg = 0;
  if (strcmp(style, "reduce") == 0) {
    if (narg < 5) utils::missing_cmd_args(FLERR, "compute reduce", error);
    iarg = 3;
  } else if (strcmp(style, "reduce/region") == 0) {
    if (narg < 6) utils::missing_cmd_args(FLERR, "compute reduce/region", error);
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
    error->all(FLERR, "Unknown compute {} mode: {}", style, arg[iarg]);
  iarg++;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR, narg - iarg, &arg[iarg], 1, earg, lmp);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values

  values.clear();
  nvalues = 0;
  for (int iarg = 0; iarg < nargnew; ++iarg) {
    value_t val;

    val.id = "";
    val.flavor = 0;
    val.val.c = nullptr;

    if (strcmp(arg[iarg], "x") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 0;
    } else if (strcmp(arg[iarg], "y") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 1;
    } else if (strcmp(arg[iarg], "z") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 2;

    } else if (strcmp(arg[iarg], "vx") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 0;
    } else if (strcmp(arg[iarg], "vy") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 1;
    } else if (strcmp(arg[iarg], "vz") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 2;

    } else if (strcmp(arg[iarg], "fx") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 0;
    } else if (strcmp(arg[iarg], "fy") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 1;
    } else if (strcmp(arg[iarg], "fz") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 2;

    } else {

      ArgInfo argi(arg[iarg]);

      val.which = argi.get_type();
      val.argindex = argi.get_index1();
      val.id = argi.get_name();

      if ((val.which == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR, "Illegal compute {} argument: {}", style, arg[iarg]);

      if (val.which == ArgInfo::NONE) break;
    }
    values.push_back(val);
  }

  // optional args

  nvalues = values.size();
  replace = new int[nvalues];
  for (int i = 0; i < nvalues; ++i) replace[i] = -1;
  std::string mycmd = "compute ";
  mycmd += style;

  for (int iarg = nvalues; iarg < nargnew; iarg++) {
    if (strcmp(arg[iarg], "replace") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, mycmd + " replace", error);
      if (mode != MINN && mode != MAXX)
        error->all(FLERR, "Compute {} replace requires min or max mode", style);
      int col1 = utils::inumeric(FLERR, arg[iarg + 1], false, lmp) - 1;
      int col2 = utils::inumeric(FLERR, arg[iarg + 2], false, lmp) - 1;
      if ((col1 < 0) || (col1 >= nvalues))
        error->all(FLERR, "Invalid compute {} replace first column index {}", style, col1);
      if ((col2 < 0) || (col2 >= nvalues))
        error->all(FLERR, "Invalid compute {} replace second column index {}", style, col2);
      if (col1 == col2) error->all(FLERR, "Compute {} replace columns must be different");
      if ((replace[col1] >= 0) || (replace[col2] >= 0))
        error->all(FLERR, "Compute {} replace column already used for another replacement");
      replace[col1] = col2;
      iarg += 2;
    } else
      error->all(FLERR, "Unknown compute {} keyword: {}", style, arg[iarg]);
  }

  // delete replace list if not set

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

  for (auto &val : values) {
    if (val.which == ArgInfo::X || val.which == ArgInfo::V || val.which == ArgInfo::F)
      val.flavor = PERATOM;

    else if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR, "Compute ID {} for compute {} does not exist", val.id, style);
      if (val.val.c->peratom_flag) {
        val.flavor = PERATOM;
        if (val.argindex == 0 && val.val.c->size_peratom_cols != 0)
          error->all(FLERR, "Compute {} compute {} does not calculate a per-atom vector", style,
                     val.id);
        if (val.argindex && val.val.c->size_peratom_cols == 0)
          error->all(FLERR, "Compute {} compute {} does not calculate a per-atom array", style,
                     val.id);
        if (val.argindex && val.argindex > val.val.c->size_peratom_cols)
          error->all(FLERR, "Compute {} compute {} array is accessed out-of-range", style, val.id);
      } else if (val.val.c->local_flag) {
        val.flavor = LOCAL;
        if (val.argindex == 0 && val.val.c->size_local_cols != 0)
          error->all(FLERR, "Compute {} compute {} does not calculate a local vector", style,
                     val.id);
        if (val.argindex && val.val.c->size_local_cols == 0)
          error->all(FLERR, "Compute {} compute {} does not calculate a local array", style,
                     val.id);
        if (val.argindex && val.argindex > val.val.c->size_local_cols)
          error->all(FLERR, "Compute {} compute {} array is accessed out-of-range", style, val.id);
      } else
        error->all(FLERR, "Compute {} compute {} calculates global values", style, val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for compute {} does not exist", val.id, style);
      if (val.val.f->peratom_flag) {
        val.flavor = PERATOM;
        if (val.argindex == 0 && (val.val.f->size_peratom_cols != 0))
          error->all(FLERR, "Compute {} fix {} does not calculate a per-atom vector", style,
                     val.id);
        if (val.argindex && (val.val.f->size_peratom_cols == 0))
          error->all(FLERR, "Compute {} fix {} does not calculate a per-atom array", style, val.id);
        if (val.argindex && (val.argindex > val.val.f->size_peratom_cols))
          error->all(FLERR, "Compute {} fix {} array is accessed out-of-range", style, val.id);
      } else if (val.val.f->local_flag) {
        val.flavor = LOCAL;
        if (val.argindex == 0 && (val.val.f->size_local_cols != 0))
          error->all(FLERR, "Compute {} fix {} does not calculate a local vector", style, val.id);
        if (val.argindex && (val.val.f->size_local_cols == 0))
          error->all(FLERR, "Compute {} fix {} does not calculate a local array", style, val.id);
        if (val.argindex && (val.argindex > val.val.f->size_local_cols))
          error->all(FLERR, "Compute {} fix {} array is accessed out-of-range", style, val.id);
      } else
        error->all(FLERR, "Compute {} fix {} calculates global values", style, val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, "Variable name {} for compute {} does not exist", val.id, style);
      if (input->variable->atomstyle(val.val.v) == 0)
        error->all(FLERR, "Compute {} variable {} is not atom-style variable", style, val.id);
      val.flavor = PERATOM;
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

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR, "Compute ID {} for compute {} does not exist", val.id, style);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR, "Fix ID {} for compute {} does not exist", val.id, style);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, "Variable name {} for compute {} does not exist", val.id, style);
    }
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
          pairme.proc = comm->me;
          MPI_Allreduce(&pairme, &pairall, 1, MPI_DOUBLE_INT, MPI_MINLOC, world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (comm->me == owner[replace[m]]) vector[m] = compute_one(m, indices[replace[m]]);
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
          pairme.proc = comm->me;
          MPI_Allreduce(&pairme, &pairall, 1, MPI_DOUBLE_INT, MPI_MAXLOC, world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (comm->me == owner[replace[m]]) vector[m] = compute_one(m, indices[replace[m]]);
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
  // invoke the appropriate attribute,compute,fix,variable
  // for flag = -1, compute scalar quantity by scanning over atom properties
  // only include atoms in group for atom properties and per-atom quantities

  index = -1;
  auto &val = values[m];

  // initialization in case it has not yet been run, e.g. when
  // the compute was invoked right after it has been created

  if ((val.which == ArgInfo::COMPUTE) || (val.which == ArgInfo::FIX)) {
    if (val.val.c == nullptr) init();
  }

  int aidx = val.argindex;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double one = 0.0;
  if (mode == MINN) one = BIG;
  if (mode == MAXX) one = -BIG;

  if (val.which == ArgInfo::X) {
    double **x = atom->x;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, x[i][aidx], i);
    } else
      one = x[flag][aidx];
  } else if (val.which == ArgInfo::V) {
    double **v = atom->v;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, v[i][aidx], i);
    } else
      one = v[flag][aidx];
  } else if (val.which == ArgInfo::F) {
    double **f = atom->f;
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, f[i][aidx], i);
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
        double *comp_vec = val.val.c->vector_atom;
        int n = nlocal;
        if (flag < 0) {
          for (int i = 0; i < n; i++)
            if (mask[i] & groupbit) combine(one, comp_vec[i], i);
        } else
          one = comp_vec[flag];
      } else {
        double **carray_atom = val.val.c->array_atom;
        int n = nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (int i = 0; i < n; i++)
            if (mask[i] & groupbit) combine(one, carray_atom[i][aidxm1], i);
        } else
          one = carray_atom[flag][aidxm1];
      }

    } else if (val.flavor == LOCAL) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_LOCAL)) {
        val.val.c->compute_local();
        val.val.c->invoked_flag |= Compute::INVOKED_LOCAL;
      }

      if (aidx == 0) {
        double *comp_vec = val.val.c->vector_local;
        int n = val.val.c->size_local_rows;
        if (flag < 0)
          for (int i = 0; i < n; i++) combine(one, comp_vec[i], i);
        else
          one = comp_vec[flag];
      } else {
        double **carray_local = val.val.c->array_local;
        int n = val.val.c->size_local_rows;
        int aidxm1 = aidx - 1;
        if (flag < 0)
          for (int i = 0; i < n; i++) combine(one, carray_local[i][aidxm1], i);
        else
          one = carray_local[flag][aidxm1];
      }
    }

    // access fix fields, check if fix frequency is a match

  } else if (val.which == ArgInfo::FIX) {
    if (update->ntimestep % val.val.f->peratom_freq)
      error->all(FLERR, "Fix {} used in compute {} not computed at compatible time", val.id, style);

    if (val.flavor == PERATOM) {
      if (aidx == 0) {
        double *fix_vector = val.val.f->vector_atom;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) combine(one, fix_vector[i], i);
        } else
          one = fix_vector[flag];
      } else {
        double **fix_array = val.val.f->array_atom;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) combine(one, fix_array[i][aidxm1], i);
        } else
          one = fix_array[flag][aidxm1];
      }

    } else if (val.flavor == LOCAL) {
      if (aidx == 0) {
        double *fix_vector = val.val.f->vector_local;
        int n = val.val.f->size_local_rows;
        if (flag < 0)
          for (int i = 0; i < n; i++) combine(one, fix_vector[i], i);
        else
          one = fix_vector[flag];
      } else {
        double **fix_array = val.val.f->array_local;
        int n = val.val.f->size_local_rows;
        int aidxm1 = aidx - 1;
        if (flag < 0)
          for (int i = 0; i < n; i++) combine(one, fix_array[i][aidxm1], i);
        else
          one = fix_array[flag][aidxm1];
      }
    }

    // evaluate atom-style variable

  } else if (val.which == ArgInfo::VARIABLE) {
    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom, maxatom, "reduce:varatom");
    }

    input->variable->compute_atom(val.val.v, igroup, varatom, 1, 0);
    if (flag < 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) combine(one, varatom[i], i);
    } else
      one = varatom[flag];
  }

  return one;
}

/* ---------------------------------------------------------------------- */

bigint ComputeReduce::count(int m)
{
  auto &val = values[m];
  if ((val.which == ArgInfo::X) || (val.which == ArgInfo::V) || (val.which == ArgInfo::F))
    return group->count(igroup);
  else if (val.which == ArgInfo::COMPUTE) {
    if (val.flavor == PERATOM) {
      return group->count(igroup);
    } else if (val.flavor == LOCAL) {
      bigint ncount = val.val.c->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
      return ncountall;
    }
  } else if (val.which == ArgInfo::FIX) {
    if (val.flavor == PERATOM) {
      return group->count(igroup);
    } else if (val.flavor == LOCAL) {
      bigint ncount = val.val.f->size_local_rows;
      bigint ncountall;
      MPI_Allreduce(&ncount, &ncountall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
      return ncountall;
    }
  } else if (val.which == ArgInfo::VARIABLE)
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

// clang-format off
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

#include "compute_global_atom.h"

#include "arg_info.h"
#include "atom.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGlobalAtom::ComputeGlobalAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), indices(nullptr), varatom(nullptr), vecglobal(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR,"compute global/atom", error);

  peratom_flag = 1;

  // process index arg

  int iarg = 3;
  ArgInfo argi(arg[iarg]);

  reference.which = argi.get_type();
  reference.argindex = argi.get_index1();
  reference.id = argi.get_name();

  if ((reference.which == ArgInfo::UNKNOWN) || (reference.which == ArgInfo::NONE)
      || (argi.get_dim() > 1))
    error->all(FLERR,"Illegal compute global/atom index property: {}", arg[iarg]);

  iarg++;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR,narg-iarg,&arg[iarg],1,earg,lmp);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  values.clear();

  for (iarg = 0; iarg < nargnew; iarg++) {
    ArgInfo argi2(arg[iarg]);

    value_t val;
    val.which = argi2.get_type();
    val.argindex = argi2.get_index1();
    val.id = argi2.get_name();
    val.val.c = nullptr;

    if ((val.which == ArgInfo::UNKNOWN) || (val.which == ArgInfo::NONE)
        || (argi2.get_dim() > 1))
      error->all(FLERR,"Illegal compute global/atom global property: {}", arg[iarg]);

    values.push_back(val);
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check for both, index arg and values

  if (reference.which == ArgInfo::COMPUTE) {
    reference.val.c = modify->get_compute_by_id(reference.id);
    if (!reference.val.c)
      error->all(FLERR,"Compute ID {} for compute global/atom index", reference.id);

    if (!reference.val.c->peratom_flag)
      error->all(FLERR,"Compute global/atom compute {} does not calculate a per-atom "
                 "vector or array", reference.id);
    if ((reference.argindex == 0) && (reference.val.c->size_peratom_cols != 0))
      error->all(FLERR,"Compute global/atom compute {} does not calculate a per-atom "
                 "vector", reference.id);
    if (reference.argindex && (reference.val.c->size_peratom_cols == 0))
      error->all(FLERR,"Compute global/atom compute does not calculate a per-atom "
                 "array", reference.id);
    if (reference.argindex && (reference.argindex > reference.val.c->size_peratom_cols))
      error->all(FLERR, "Compute global/atom compute array {} is accessed out-of-range",
        reference.id);

  } else if (reference.which == ArgInfo::FIX) {
    reference.val.f =modify->get_fix_by_id(reference.id);
    if (!reference.val.f)
      error->all(FLERR,"Fix ID {} for compute global/atom does not exist", reference.id);
    if (!reference.val.f->peratom_flag)
      error->all(FLERR,"Compute global/atom fix {} does not calculate a per-atom vector "
                 "or array", reference.id);
    if (reference.argindex == 0 && (reference.val.f->size_peratom_cols != 0))
      error->all(FLERR,"Compute global/atom fix {} does not calculate a per-atom vector",
                 reference.id);
    if (reference.argindex && (reference.val.f->size_peratom_cols == 0))
      error->all(FLERR,"Compute global/atom fix {} does not calculate a per-atom array",
                 reference.id);
    if (reference.argindex && (reference.argindex > reference.val.f->size_peratom_cols))
      error->all(FLERR, "Compute global/atom fix {} array is accessed out-of-range", reference.id);

  } else if (reference.which == ArgInfo::VARIABLE) {
    reference.val.v = input->variable->find(reference.id.c_str());
    if (reference.val.v < 0)
      error->all(FLERR,"Variable name {} for compute global/atom index does not exist",
                 reference.id);
    if (input->variable->atomstyle(reference.val.v) == 0)
      error->all(FLERR,"Compute global/atom index variable {} is not atom-style variable",
        reference.id);
  }

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for compute global/atom does not exist", val.id);
      if (val.argindex == 0) {
        if (!val.val.c->vector_flag)
          error->all(FLERR,"Compute ID {} for global/atom compute does not calculate "
                     "a global vector", val.id);
      } else {
        if (!val.val.c->array_flag)
          error->all(FLERR,"Compute ID {} for global/atom compute does not calculate "
                     "a global array", val.id);
        if (val.argindex > val.val.c->size_array_cols)
          error->all(FLERR,"Compute global/atom compute {} array is accessed out-of-range", val.id);
      }

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for compute global/atom does not exist", val.id);
      if (val.argindex == 0) {
        if (!val.val.f->vector_flag)
          error->all(FLERR,"Fix ID {} for compute global/atom compute does not calculate "
                     "a global vector", val.id);
      } else {
        if (!val.val.f->array_flag)
          error->all(FLERR,"Fix ID {} for compute global/atom compute does not calculate "
                     "a global array", val.id);
        if (val.argindex > val.val.f->size_array_cols)
          error->all(FLERR,"Compute global/atom fix {} array is accessed out-of-range", val.id);
      }

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for compute global/atom does not exist", val.id);
      if (input->variable->vectorstyle(val.val.v) == 0)
        error->all(FLERR,"Compute global/atom variable {} is not vector-style variable", val.id);
    }
  }

  // this compute produces either a peratom vector or array

  if (values.size() == 1) size_peratom_cols = 0;
  else size_peratom_cols = values.size();

  nmax = maxvector = 0;
  vector_atom = nullptr;
  array_atom = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeGlobalAtom::~ComputeGlobalAtom()
{
  memory->destroy(indices);
  memory->destroy(varatom);
  memory->destroy(vecglobal);
  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeGlobalAtom::init()
{
  // store references of all computes, fixes, or variables

  if (reference.which == ArgInfo::COMPUTE) {
    reference.val.c = modify->get_compute_by_id(reference.id);
    if (!reference.val.c)
      error->all(FLERR,"Compute ID {} for compute global/atom index does not exist", reference.id);
  } else if (reference.which == ArgInfo::FIX) {
    reference.val.f = modify->get_fix_by_id(reference.id);
    if (reference.val.f)
      error->all(FLERR,"Fix ID {} for compute global/atom index does not exist", reference.id);
  } else if (reference.which == ArgInfo::VARIABLE) {
    reference.val.v = input->variable->find(reference.id.c_str());
    if (reference.val.v < 0)
      error->all(FLERR,"Variable name {} for compute global/atom index does not exist",
                 reference.id);
  }

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for compute global/atom does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR,"Fix ID {} for compute global/atom does not exist", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for compute global/atom does not exist", val.id);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeGlobalAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow indices and output vector or array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(indices);
    memory->create(indices,nmax,"global/atom:indices");
    if (reference.which == ArgInfo::VARIABLE) {
      memory->destroy(varatom);
      memory->create(varatom,nmax,"global/atom:varatom");
    }
    if (values.size() == 1) {
      memory->destroy(vector_atom);
      memory->create(vector_atom,nmax,"global/atom:vector_atom");
    } else {
      memory->destroy(array_atom);
      memory->create(array_atom,nmax,values.size(),"global/atom:array_atom");
    }
  }

  // setup current peratom indices
  // integer indices are rounded down from double values
  // indices are decremented from 1 to N -> 0 to N-1

  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  if (reference.which == ArgInfo::COMPUTE) {

    if (!(reference.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
      reference.val.c->compute_peratom();
      reference.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
    }

    if (reference.argindex == 0) {
      double *compute_vector = reference.val.c->vector_atom;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (compute_vector[i]) - 1;
    } else {
      double **compute_array = reference.val.c->array_atom;
      int im1 = reference.argindex - 1;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (compute_array[i][im1]) - 1;
    }

  } else if (reference.which == ArgInfo::FIX) {
    if (update->ntimestep % reference.val.f->peratom_freq)
      error->all(FLERR,"Fix {} used in compute global/atom not computed at compatible time",
                 reference.id);

    if (reference.argindex == 0) {
      double *fix_vector = reference.val.f->vector_atom;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (fix_vector[i]) - 1;
    } else {
      double **fix_array = reference.val.f->array_atom;
      int im1 = reference.argindex - 1;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (fix_array[i][im1]) - 1;
    }

  } else if (reference.which == ArgInfo::VARIABLE) {
    input->variable->compute_atom(reference.val.v, igroup, varatom, 1, 0);
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        indices[i] = static_cast<int> (varatom[i]) - 1;
  }

  // loop over values to fill output vector or array

  int m = 0;
  for (auto &val : values) {

    // output = vector

    if (val.argindex == 0) {
      int vmax;
      double *source;

      if (val.which == ArgInfo::COMPUTE) {

        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }

        source = val.val.c->vector;
        vmax = val.val.c->size_vector;

      } else if (val.which == ArgInfo::FIX) {
        if (update->ntimestep % val.val.f->peratom_freq)
          error->all(FLERR,"Fix {} used in compute global/atom not computed at compatible time",
            val.id);
        vmax = reference.val.f->size_vector;

        if (vmax > maxvector) {
          maxvector = vmax;
          memory->destroy(vecglobal);
          memory->create(vecglobal,maxvector,"global/atom:vecglobal");
        }

        for (int i = 0; i < vmax; i++)
          vecglobal[i] = val.val.f->compute_vector(i);

        source = vecglobal;

      } else if (val.which == ArgInfo::VARIABLE) {
        vmax = input->variable->compute_vector(val.val.v, &source);
      }

      if (values.size() == 1) {
        for (int i = 0; i < nlocal; i++) {
          vector_atom[i] = 0.0;
          if (mask[i] & groupbit) {
            int j = indices[i];
            if (j >= 0 && j < vmax) vector_atom[i] = source[j];
          }
        }
      } else {
        for (int i = 0; i < nlocal; i++) {
          array_atom[i][m] = 0.0;
          if (mask[i] & groupbit) {
            int j = indices[i];
            if (j >= 0 && j < vmax) array_atom[i][m] = source[j];
          }
        }
      }

    // output = array

    } else {
      int vmax;
      double *source;
      int col = val.argindex - 1;

      if (val.which == ArgInfo::COMPUTE) {

        if (!(val.val.c->invoked_flag & Compute::INVOKED_ARRAY)) {
          val.val.c->compute_array();
          val.val.c->invoked_flag |= Compute::INVOKED_ARRAY;
        }

        double **compute_array = val.val.c->array;
        vmax = val.val.c->size_array_rows;

        if (vmax > maxvector) {
          maxvector = vmax;
          memory->destroy(vecglobal);
          memory->create(vecglobal,maxvector,"global/atom:vecglobal");
        }

        for (int i = 0; i < vmax; i++)
          vecglobal[i] = compute_array[i][col];

        source = vecglobal;

      } else if (val.which == ArgInfo::FIX) {
        if (update->ntimestep % val.val.f->peratom_freq)
          error->all(FLERR,"Fix {} used in compute global/atom not computed at compatible time",
            val.id);
        vmax = val.val.f->size_array_rows;

        if (vmax > maxvector) {
          maxvector = vmax;
          memory->destroy(vecglobal);
          memory->create(vecglobal,maxvector,"global/atom:vecglobal");
        }

        for (int i = 0; i < vmax; i++)
          vecglobal[i] = val.val.f->compute_array(i,col);

        source = vecglobal;

      } else if (val.which == ArgInfo::VARIABLE) {
        vmax = input->variable->compute_vector(val.val.v, &source);
      }

      if (values.size() == 1) {
        for (int i = 0; i < nlocal; i++) {
          vector_atom[i] = 0.0;
          if (mask[i] & groupbit) {
            int j = indices[i];
            if (j >= 0 && j < vmax) vector_atom[i] = source[j];
          }
        }
      } else {
        for (int i = 0; i < nlocal; i++) {
          array_atom[i][m] = 0.0;
          if (mask[i] & groupbit) {
            int j = indices[i];
            if (j >= 0 && j < vmax) array_atom[i][m] = source[j];
          }
        }
      }
    }
  }
  ++m;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeGlobalAtom::memory_usage()
{
  double bytes = (double)nmax*values.size() * sizeof(double);
  bytes += (double)nmax * sizeof(int);                    // indices
  if (varatom) bytes += (double)nmax * sizeof(double);    // varatom
  bytes += (double)maxvector * sizeof(double);            // vecglobal
  return bytes;
}

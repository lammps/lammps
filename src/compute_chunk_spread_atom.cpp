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

#include "compute_chunk_spread_atom.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeChunkSpreadAtom::
ComputeChunkSpreadAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), idchunk(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal compute chunk/spread/atom command");

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);
  init_chunk();

  // expand args if any have wildcard character "*"

  int iarg = 4;
  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR,narg-iarg,&arg[iarg],1,earg,lmp);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values

  values.clear();
  for (iarg = 0; iarg < nargnew; ++iarg) {
    ArgInfo argi(arg[iarg], ArgInfo::COMPUTE|ArgInfo::FIX);

    value_t val;
    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();
    val.val.c = nullptr;

    if ((val.which == ArgInfo::UNKNOWN) || (val.which == ArgInfo::NONE) || (argi.get_dim() > 1))
      error->all(FLERR,"Illegal compute chunk/spread/atom argument: {}", arg[iarg]);

    values.push_back(val);
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // for compute, must calculate per-chunk values, i.e. style ends in "/chunk"
  // for fix, assume a global vector or array is per-chunk data

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      auto icompute = modify->get_compute_by_id(val.id);
      if (!icompute)
        error->all(FLERR,"Compute ID {} for compute chunk/spread/atom does not exist", val.id);

      if (!utils::strmatch(icompute->style,"/chunk$"))
        error->all(FLERR,"Compute chunk/spread/atom compute {} does not calculate per-chunk values",
                   val.id);

      if (val.argindex == 0) {
        if (!icompute->vector_flag)
          error->all(FLERR,"Compute chunk/spread/atom compute {} does not calculate global vector",
                     val.id);
      } else {
        if (!icompute->array_flag)
          error->all(FLERR,"Compute chunk/spread/atom compute {} does not calculate global array",
                     val.id);
        if (val.argindex > icompute->size_array_cols)
          error->all(FLERR,"Compute chunk/spread/atom compute {} array is accessed out-of-range",
                     val.id);
      }
      val.val.c = icompute;

    } else if (val.which == ArgInfo::FIX) {
      auto ifix = modify->get_fix_by_id(val.id);
      if (!ifix)
        error->all(FLERR,"Fix ID {} for compute chunk/spread/atom does not exist", val.id);
      if (val.argindex == 0) {
        if (!ifix->vector_flag)
          error->all(FLERR,"Compute chunk/spread/atom {} fix does not calculate global vector",
                     val.id);
      } else {
        if (!ifix->array_flag)
          error->all(FLERR,"Compute chunk/spread/atom {} fix does not calculate global array",
                     val.id);
        if (val.argindex > ifix->size_array_cols)
          error->all(FLERR,"Compute chunk/spread/atom fix {} array is accessed out-of-range",
                     val.id);
      }
    }
  }

  // this compute produces a peratom vector or array

  peratom_flag = 1;
  if (values.size() == 1) size_peratom_cols = 0;
  else size_peratom_cols = values.size();

  // per-atom vector or array

  nmax = 0;
  vector_atom = nullptr;
  array_atom = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeChunkSpreadAtom::~ComputeChunkSpreadAtom()
{
  delete[] idchunk;

  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeChunkSpreadAtom::init()
{
  init_chunk();

  // store references of all computes and fixes

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for compute chunk/spread/atom does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR,"Fix ID {} for compute chunk/spread/atom does not exist", val.id);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeChunkSpreadAtom::init_chunk()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Chunk/atom compute {} does not exist for compute chunk/spread/atom "
               "or is of invalid style", idchunk);
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute chunk/spread/atom {} does not use chunk/atom compute", idchunk);
}

/* ---------------------------------------------------------------------- */

void ComputeChunkSpreadAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local vector_atom or array_atom if necessary

  if (atom->nmax > nmax) {
    if (values.size() == 1) {
      memory->destroy(vector_atom);
      nmax = atom->nmax;
      memory->create(vector_atom,nmax,"chunk/spread/atom:vector_atom");
    } else {
      memory->destroy(array_atom);
      nmax = atom->nmax;
      memory->create(array_atom,nmax,values.size(),"chunk/spread/atom:array_atom");
    }
  }

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  int nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  // loop over values, access compute or fix
  // loop over atoms, use chunk ID of each atom to store value from compute/fix

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int index,nstride;
  double *ptr;

  int m = 0;
  for (auto &val : values) {

    // copy compute/fix values into vector_atom or array_atom
    // nstride between values for each atom

    if (values.size() == 1) {
      ptr = vector_atom;
      nstride = 1;
    } else {
      ptr = &array_atom[0][m];
      nstride = values.size();
    }

    // invoke compute if not previously invoked

    if (val.which == ArgInfo::COMPUTE) {
      Compute *compute = val.val.c;

      if (val.argindex == 0) {
        if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        double *cvector = compute->vector;
        for (int i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk) continue;
          *ptr = cvector[index];
        }

      } else {
        if (!(compute->invoked_flag & Compute::INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= Compute::INVOKED_ARRAY;
        }
        int icol = val.argindex-1;
        double **carray = compute->array;
        for (int i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk) continue;
          *ptr = carray[index][icol];
        }
      }

    // access fix data, check if fix frequency is a match
    // are assuming the fix global vector/array is per-chunk data
    // check if index exceeds fix output length/rows

    } else if (val.which == ArgInfo::FIX) {
      Fix *fix = val.val.f;
      if (update->ntimestep % fix->global_freq)
        error->all(FLERR,"Fix {} used in compute chunk/spread/atom not computed at compatible time",
                   val.id);

      if (val.argindex == 0) {
        int nfix = fix->size_vector;
        for (int i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk || index >= nfix) continue;
          *ptr = fix->compute_vector(index);
        }

      } else {
        int icol = val.argindex-1;
        int nfix = fix->size_array_rows;
        for (int i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk || index >= nfix) continue;
          *ptr = fix->compute_array(index,icol);
        }
      }
    }
    ++m;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeChunkSpreadAtom::memory_usage()
{
  double bytes = (double)nmax*values.size() * sizeof(double);
  bytes += values.size() * sizeof(value_t);
  return bytes;
}

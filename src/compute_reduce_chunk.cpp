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

#include "compute_reduce_chunk.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{ SUM, MINN, MAXX };

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduceChunk::ComputeReduceChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), idchunk(nullptr), vlocal(nullptr), vglobal(nullptr),
  alocal(nullptr), aglobal(nullptr), varatom(nullptr), cchunk(nullptr), ichunk(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR,"compute reduce/chunk", error);

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);
  init_chunk();

  // mode

  if (strcmp(arg[4],"sum") == 0) mode = SUM;
  else if (strcmp(arg[4],"min") == 0) mode = MINN;
  else if (strcmp(arg[4],"max") == 0) mode = MAXX;
  else error->all(FLERR,"Unknown compute reduce/chunk mode: {}", arg[4]);

  int iarg = 5;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR,narg-iarg,&arg[iarg],1,earg,lmp);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values

  values.clear();
  for (iarg = 0; iarg < nargnew; iarg++) {
    ArgInfo argi(arg[iarg]);

    value_t val;
    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();
    val.val.c = nullptr;

    if ((val.which == ArgInfo::UNKNOWN) || (val.which == ArgInfo::NONE) || (argi.get_dim() > 1))
      error->all(FLERR,"Illegal compute reduce/chunk argument: {}", arg[iarg]);

    values.push_back(val);
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // error check

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for compute reduce/chunk does not exist", val.id);
      if (!val.val.c->peratom_flag)
        error->all(FLERR,"Compute reduce/chunk compute {} does not calculate per-atom values",
                   val.id);
      if ((val.argindex == 0) && (val.val.c->size_peratom_cols != 0))
        error->all(FLERR,"Compute reduce/chunk compute {} does not calculate a per-atom vector",
                   val.id);
      if (val.argindex && (val.val.c->size_peratom_cols == 0))
        error->all(FLERR,"Compute reduce/chunk compute {} does not calculate a per-atom array",
                   val.id);
      if (val.argindex && (val.argindex > val.val.c->size_peratom_cols))
        error->all(FLERR, "Compute reduce/chunk compute array {} is accessed out-of-range", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR,"Fix ID {} for compute reduce/chunk does not exist", val.id);
      if (!val.val.f->peratom_flag)
        error->all(FLERR,"Compute reduce/chunk fix {} does not calculate per-atom values", val.id);
      if ((val.argindex == 0) && (val.val.f->size_peratom_cols != 0))
        error->all(FLERR,"Compute reduce/chunk fix {} does not calculate a per-atom vector", val.id);
      if (val.argindex && (val.val.f->size_peratom_cols == 0))
        error->all(FLERR,"Compute reduce/chunk fix {} does not calculate a per-atom array", val.id);
      if (val.argindex && (val.argindex > val.val.f->size_peratom_cols))
        error->all(FLERR,"Compute reduce/chunk fix {} array is accessed out-of-range", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for compute reduce/chunk does not exist", val.id);
      if (input->variable->atomstyle(val.val.v) == 0)
        error->all(FLERR,"Compute reduce/chunk variable is not atom-style variable");
    }
  }

  // this compute produces either a vector or array

  if (values.size() == 1) {
    vector_flag = 1;
    size_vector_variable = 1;
    extvector = 0;
  } else {
    array_flag = 1;
    size_array_rows_variable = 1;
    size_array_cols = values.size();
    extarray = 0;
  }

  // setup

  if (mode == SUM) initvalue = 0.0;
  else if (mode == MINN) initvalue = BIG;
  else if (mode == MAXX) initvalue = -BIG;

  maxchunk = 0;
  vlocal = vglobal = nullptr;
  alocal = aglobal = nullptr;

  maxatom = 0;
  varatom = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeReduceChunk::~ComputeReduceChunk()
{
  delete[] idchunk;

  memory->destroy(vlocal);
  memory->destroy(vglobal);
  memory->destroy(alocal);
  memory->destroy(aglobal);

  memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::init()
{
  init_chunk();

  // set indices of all computes,fixes,variables

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for compute reduce/chunk does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR,"Fix ID {} for compute reduce/chunk does not exist", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for compute reduce/chunk does not exist", val.id);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::init_chunk()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Compute chunk/atom {} does not exist or is incorrect style for "
               "compute reduce/chunk", idchunk);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::compute_vector()
{
  invoked_vector = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  ichunk = cchunk->ichunk;
  if (!nchunk) return;

  size_vector = nchunk;

  if (nchunk > maxchunk) {
    memory->destroy(vlocal);
    memory->destroy(vglobal);
    maxchunk = nchunk;
    memory->create(vlocal,maxchunk,"reduce/chunk:vlocal");
    memory->create(vglobal,maxchunk,"reduce/chunk:vglobal");
    vector = vglobal;
  }

  // perform local reduction of single peratom value

  compute_one(0,vlocal,1);

  // reduce the per-chunk values across all procs

  if (mode == SUM)
    MPI_Allreduce(vlocal,vglobal,nchunk,MPI_DOUBLE,MPI_SUM,world);
  else if (mode == MINN)
    MPI_Allreduce(vlocal,vglobal,nchunk,MPI_DOUBLE,MPI_MIN,world);
  else if (mode == MAXX)
    MPI_Allreduce(vlocal,vglobal,nchunk,MPI_DOUBLE,MPI_MAX,world);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::compute_array()
{
  invoked_array = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  ichunk = cchunk->ichunk;
  if (!nchunk) return;

  size_array_rows = nchunk;

  if (nchunk > maxchunk) {
    memory->destroy(alocal);
    memory->destroy(aglobal);
    maxchunk = nchunk;
    memory->create(alocal,maxchunk,values.size(),"reduce/chunk:alocal");
    memory->create(aglobal,maxchunk,values.size(),"reduce/chunk:aglobal");
    array = aglobal;
  }

  // perform local reduction of all peratom values

  for (std::size_t m = 0; m < values.size(); m++) compute_one(m,&alocal[0][m],values.size());

  // reduce the per-chunk values across all procs

  if (mode == SUM)
    MPI_Allreduce(&alocal[0][0],&aglobal[0][0],nchunk*values.size(),MPI_DOUBLE,MPI_SUM,world);
  else if (mode == MINN)
    MPI_Allreduce(&alocal[0][0],&aglobal[0][0],nchunk*values.size(),MPI_DOUBLE,MPI_MIN,world);
  else if (mode == MAXX)
    MPI_Allreduce(&alocal[0][0],&aglobal[0][0],nchunk*values.size(),MPI_DOUBLE,MPI_MAX,world);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::compute_one(int m, double *vchunk, int nstride)
{
  // initialize per-chunk values in accumulation vector

  for (std::size_t i = 0; i < values.size()*nchunk; i += nstride) vchunk[i] = initvalue;

  // loop over my atoms
  // use peratom input and chunk ID of each atom to update vector

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  auto &val = values[m];
  int index = -1;

  // initialization in case it has not yet been run, e.g. when
  // the compute was invoked right after it has been created

  if (val.val.c == nullptr) init();

  if (val.which == ArgInfo::COMPUTE) {
    if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
      val.val.c->compute_peratom();
      val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
    }

    if (val.argindex == 0) {
      double *vcompute = val.val.c->vector_atom;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],vcompute[i]);
      }
    } else {
      double **acompute = val.val.c->array_atom;
      int argindexm1 = val.argindex - 1;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],acompute[i][argindexm1]);
      }
    }

  // access fix fields, check if fix frequency is a match

  } else if (val.which == ArgInfo::FIX) {
    if (update->ntimestep % val.val.f->peratom_freq)
      error->all(FLERR,"Fix used in compute reduce/chunk not computed at compatible time");

    if (val.argindex == 0) {
      double *vfix = val.val.f->vector_atom;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],vfix[i]);
      }
    } else {
      double **afix = val.val.f->array_atom;
      int argindexm1 = val.argindex - 1;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],afix[i][argindexm1]);
      }
    }

  // evaluate atom-style variable

  } else if (val.which == ArgInfo::VARIABLE) {
    if (atom->nmax > maxatom) {
      memory->destroy(varatom);
      maxatom = atom->nmax;
      memory->create(varatom,maxatom,"reduce/chunk:varatom");
    }

    input->variable->compute_atom(val.val.v,igroup,varatom,1,0);
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      index = ichunk[i]-1;
      if (index < 0) continue;
      combine(vchunk[index*nstride],varatom[i]);
    }
  }
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
------------------------------------------------------------------------- */

void ComputeReduceChunk::combine(double &one, double two)
{
  if (mode == SUM) one += two;
  else if (mode == MINN) {
    if (two < one) one = two;
  } else if (mode == MAXX) {
    if (two > one) one = two;
  }
}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods ensure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeReduceChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, if it still exists
------------------------------------------------------------------------- */

void ComputeReduceChunk::lock_disable()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (cchunk) cchunk->lockcount--;
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

int ComputeReduceChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeReduceChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeReduceChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeReduceChunk::memory_usage()
{
  double bytes = (bigint) maxatom * sizeof(double);
  if (values.size() == 1) bytes += (double) maxchunk * 2 * sizeof(double);
  else bytes += (double) maxchunk * values.size() * 2 * sizeof(double);
  return bytes;
}

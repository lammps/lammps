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

#include "fix_ave_chunk.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { SCALAR, VECTOR };
enum { SAMPLE, ALL };
enum { NOSCALE, ATOM };
enum { ONE, RUNNING, WINDOW };

/* ---------------------------------------------------------------------- */

FixAveChunk::FixAveChunk(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), nrepeat(0), fp(nullptr), idchunk(nullptr), varatom(nullptr),
    count_one(nullptr), count_many(nullptr), count_sum(nullptr), values_one(nullptr),
    values_many(nullptr), values_sum(nullptr), count_total(nullptr), count_list(nullptr),
    values_total(nullptr), values_list(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix ave/chunk", error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[5], false, lmp);

  idchunk = utils::strdup(arg[6]);

  global_freq = nfreq;
  no_change_box = 1;
  time_depend = 1;

  char *group = arg[1];

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR, narg - 7, &arg[7], 1, earg, lmp);

  if (earg != &arg[7]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  densityflag = 0;

  int iarg = 0;
  values.clear();
  while (iarg < nargnew) {

    value_t val;
    val.id = "";
    val.val.c = nullptr;

    if (strcmp(arg[iarg],"vx") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 0;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 2;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 0;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 2;

    } else if (strcmp(arg[iarg],"mass") == 0) {
      val.which = ArgInfo::MASS;
      val.argindex = 0;
    } else if (strcmp(arg[iarg],"density/number") == 0) {
      densityflag = 1;
      val.which = ArgInfo::DENSITY_NUMBER;
      val.argindex = 0;
    } else if (strcmp(arg[iarg],"density/mass") == 0) {
      densityflag = 1;
      val.which = ArgInfo::DENSITY_MASS;
      val.argindex = 0;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      val.which = ArgInfo::TEMPERATURE;
      val.argindex = 0;

    } else {
      ArgInfo argi(arg[iarg]);

      if (argi.get_type() == ArgInfo::NONE) break;
      if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR,"Unknown fix ave/chunk data value: {}", arg[iarg]);

      val.which = argi.get_type();
      val.argindex = argi.get_index1();
      val.id = argi.get_name();
    }
    values.push_back(val);
    iarg++;
  }

  nvalues = values.size();
  if (nvalues == 0) error->all(FLERR, "No values in fix ave/chunk command");

  // optional args

  normflag = ALL;
  scaleflag = ATOM;
  ave = ONE;
  nwindow = 0;
  biasflag = 0;
  id_bias = nullptr;
  adof = domain->dimension;
  cdof = 0.0;
  overwrite = 0;
  format_user = nullptr;
  format = (char *) " %g";
  char *title1 = nullptr;
  char *title2 = nullptr;
  char *title3 = nullptr;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/chunk norm", error);
      if (strcmp(arg[iarg+1],"all") == 0) {
        normflag = ALL;
        scaleflag = ATOM;
      } else if (strcmp(arg[iarg+1],"sample") == 0) {
        normflag = SAMPLE;
        scaleflag = ATOM;
      } else if (strcmp(arg[iarg+1],"none") == 0) {
        normflag = SAMPLE;
        scaleflag = NOSCALE;
      } else error->all(FLERR,"Unknown fix ave/chunk norm mode: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/chunk ave", error);
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Unknown fix ave/chunk ave mode: {}", arg[iarg+1]);
      if (ave == WINDOW) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix ave/chunk ave window", error);
        nwindow = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/chunk number of windows: {}", nwindow);
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;

    } else if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/chunk bias", error);
      biasflag = 1;
      id_bias = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"adof") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/chunk adof", error);
      adof = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cdof") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix ave/chunk cdof", error);
      cdof = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if ((strcmp(arg[iarg],"file") == 0) || (strcmp(arg[iarg],"append") == 0)) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ave/chunk ")+arg[iarg], error);
      if (comm->me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == nullptr)
          error->one(FLERR, "Cannot open fix ave/chunk file {}: {}",
                     arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg)  utils::missing_cmd_args(FLERR, "fix ave/chunk format", error);
      delete[] format_user;
      format_user = utils::strdup(arg[iarg+1]);
      format = format_user;
      iarg += 2;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg)  utils::missing_cmd_args(FLERR, "fix ave/chunk title1", error);
      delete[] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg)  utils::missing_cmd_args(FLERR, "fix ave/chunk title2", error);
      delete[] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg)  utils::missing_cmd_args(FLERR, "fix ave/chunk title3", error);
      delete[] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Unknown fix ave/chunk keyword: {}", arg[iarg]);
  }

  // setup and error check

  if (nevery <= 0) error->all(FLERR,"Illegal fix ave/chunk nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR,"Illegal fix ave/chunk nrepeat value: {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR,"Illegal fix ave/chunk nfreq value: {}", nfreq);
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Inconsistent fix ave/chunk nevery/nrepeat/nfreq values");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Fix ave/chunk overwrite keyword requires ave running setting");

  if (biasflag) {
    tbias = modify->get_compute_by_id(id_bias);
    if (!tbias) error->all(FLERR,"Could not find compute ID {} for temperature bias", id_bias);
    if (tbias->tempflag == 0)
      error->all(FLERR,"Bias compute {} does not calculate temperature", id_bias);
    if (tbias->tempbias == 0)
      error->all(FLERR,"Bias compute {} does not calculate a velocity bias", id_bias);
  }

  for (auto &val : values) {

    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for fix ave/chunk does not exist",val.id);
      if (val.val.c->peratom_flag == 0)
        error->all(FLERR,"Fix ave/chunk compute {} does not calculate per-atom values",val.id);
      if (val.argindex == 0 && (val.val.c->size_peratom_cols != 0))
        error->all(FLERR,"Fix ave/chunk compute {} does not calculate a per-atom vector",val.id);
      if (val.argindex && (val.val.c->size_peratom_cols == 0))
        error->all(FLERR,"Fix ave/chunk compute {} does not calculate a per-atom array",val.id);
      if (val.argindex && (val.argindex > val.val.c->size_peratom_cols))
        error->all(FLERR,"Fix ave/chunk compute {} vector is accessed out-of-range",val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR, "Fix ID {} for fix ave/chunk does not exist",val.id);
      if (val.val.f->peratom_flag == 0)
        error->all(FLERR, "Fix ave/chunk fix {} does not calculate per-atom values",val.id);
      if (val.argindex == 0 && (val.val.f->size_peratom_cols != 0))
        error->all(FLERR, "Fix ave/chunk fix {} does not calculate a per-atom vector",val.id);
      if (val.argindex && (val.val.f->size_peratom_cols == 0))
        error->all(FLERR, "Fix ave/chunk fix {} does not calculate a per-atom array",val.id);
      if (val.argindex && val.argindex > val.val.f->size_peratom_cols)
        error->all(FLERR,"Fix ave/chunk fix {} vector is accessed out-of-range",val.id);
    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix ave/chunk does not exist",val.id);
      if (input->variable->atomstyle(val.val.v) == 0)
        error->all(FLERR,"Fix ave/chunk variable {} is not atom-style variable",val.id);
    }
  }

  // increment lock counter in compute chunk/atom
  // only if nrepeat > 1 or ave = RUNNING/WINDOW,
  //   so that locking spans multiple timesteps

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Chunk/atom compute {} does not exist or is "
               "incorrect style for fix ave/chunk",idchunk);

  if ((nrepeat > 1) || (ave == RUNNING) || (ave == WINDOW)) cchunk->lockcount++;
  lockforever = 0;

  // print file comment lines

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Chunk-averaged data for fix %s and group %s\n", id, group);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# Timestep Number-of-chunks Total-count\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else {
      int compress = cchunk->compress;
      int ncoord = cchunk->ncoord;
      if (!compress) {
        if (ncoord == 0) fprintf(fp,"# Chunk Ncount");
        else if (ncoord == 1) fprintf(fp,"# Chunk Coord1 Ncount");
        else if (ncoord == 2) fprintf(fp,"# Chunk Coord1 Coord2 Ncount");
        else if (ncoord == 3)
          fprintf(fp,"# Chunk Coord1 Coord2 Coord3 Ncount");
      } else {
        if (ncoord == 0) fprintf(fp,"# Chunk OrigID Ncount");
        else if (ncoord == 1) fprintf(fp,"# Chunk OrigID Coord1 Ncount");
        else if (ncoord == 2) fprintf(fp,"# Chunk OrigID Coord1 Coord2 Ncount");
        else if (ncoord == 3)
          fprintf(fp,"# Chunk OrigID Coord1 Coord2 Coord3 Ncount");
      }
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    }
    if (ferror(fp))
      error->one(FLERR,"Error writing file header");

    filepos = platform::ftell(fp);
  }

  delete[] title1;
  delete[] title2;
  delete[] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // this fix produces a global array
  // size_array_rows is variable and set by allocate()

  int compress = cchunk->compress;
  int ncoord = cchunk->ncoord;
  colextra = compress + ncoord;

  array_flag = 1;
  size_array_cols = colextra + 1 + nvalues;
  size_array_rows_variable = 1;
  extarray = 0;

  // initializations

  irepeat = 0;
  iwindow = window_limit = 0;
  normcount = 0;

  maxvar = 0;
  varatom = nullptr;

  count_one = count_many = count_sum = count_total = nullptr;
  count_list = nullptr;
  values_one = values_many = values_sum = values_total = nullptr;
  values_list = nullptr;

  maxchunk = 0;
  nchunk = 1;
  allocate();

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveChunk::~FixAveChunk()
{
  if (fp && comm->me == 0) fclose(fp);

  memory->destroy(varatom);
  memory->destroy(count_one);
  memory->destroy(count_many);
  memory->destroy(count_sum);
  memory->destroy(count_total);
  memory->destroy(count_list);
  memory->destroy(values_one);
  memory->destroy(values_many);
  memory->destroy(values_sum);
  memory->destroy(values_total);
  memory->destroy(values_list);

  // decrement lock counter in compute chunk/atom, it if still exists

  if (nrepeat > 1 || ave == RUNNING || ave == WINDOW) {
    cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
    if (cchunk) {
      if (ave == RUNNING || ave == WINDOW) cchunk->unlock(this);
      cchunk->lockcount--;
    }
  }

  delete[] idchunk;
  fp = nullptr;
  varatom = nullptr;
  count_one = nullptr;
  count_many = nullptr;
  count_sum = nullptr;
  count_total = nullptr;
  count_list = nullptr;
  values_one = nullptr;
  values_many = nullptr;
  values_sum = nullptr;
  values_total = nullptr;
  values_list = nullptr;
  idchunk = nullptr;
  cchunk = nullptr;
}

/* ---------------------------------------------------------------------- */

int FixAveChunk::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveChunk::init()
{
  // set indices and check validity of all computes,fixes,variables
  // check that fix frequency is acceptable

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Chunk/atom compute {} does not exist or is "
               "incorrect style for fix ave/chunk",idchunk);

  if (biasflag) {
    tbias = modify->get_compute_by_id(id_bias);
    if (!tbias)
      error->all(FLERR,"Could not find compute ID {} for temperature bias", id_bias);
  }

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR,"Compute ID {} for fix ave/chunk does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for fix ave/chunk does not exist", val.id);

      if (nevery % val.val.f->peratom_freq)
        error->all(FLERR, "Fix {} for fix ave/chunk not computed at compatible time", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix ave/chunk does not exist", val.id);
    }
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does averaging if nvalid = current timestep
   do not call setup_chunks(), even though fix ave/chunk called setup_bins()
   b/c could cause nchunk to change if Nfreq epoch crosses 2 runs
   does mean that if change_box is used between runs to change box size,
     that nchunk may not track it
------------------------------------------------------------------------- */

void FixAveChunk::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveChunk::end_of_step()
{
  int i,j,m,index;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // first sample within single Nfreq epoch
  // zero out arrays that accumulate over many samples, but not across epochs
  // invoke setup_chunks() to determine current nchunk
  //   re-allocate per-chunk arrays if needed
  // invoke lock() in two cases:
  //   if nrepeat > 1: so nchunk cannot change until Nfreq epoch is over,
  //     will be unlocked on last repeat of this Nfreq
  //   if ave = RUNNING/WINDOW and not yet locked:
  //     set forever, will be unlocked in fix destructor
  // wrap setup_chunks in clearstep/addstep b/c it may invoke computes
  //   both nevery and nfreq are future steps,
  //   since call below to cchunk->ichunk()
  //     does not re-invoke internal cchunk compute on this same step

  if (irepeat == 0) {
    if (cchunk->computeflag) modify->clearstep_compute();
    nchunk = cchunk->setup_chunks();
    if (cchunk->computeflag) {
      modify->addstep_compute(ntimestep+nevery);
      modify->addstep_compute(ntimestep+nfreq);
    }
    allocate();
    if (nrepeat > 1 && ave == ONE)
      cchunk->lock(this,ntimestep,ntimestep+((bigint)nrepeat-1)*nevery);
    else if ((ave == RUNNING || ave == WINDOW) && !lockforever) {
      cchunk->lock(this,update->ntimestep,-1);
      lockforever = 1;
    }
    for (m = 0; m < nchunk; m++) {
      count_many[m] = count_sum[m] = 0.0;
      for (i = 0; i < nvalues; i++) values_many[m][i] = 0.0;
    }

  // if any DENSITY requested, invoke setup_chunks() on each sampling step
  // nchunk will not change but bin volumes might, e.g. for NPT simulation

  } else if (densityflag) {
    cchunk->setup_chunks();
  }

  // zero out arrays for one sample

  for (m = 0; m < nchunk; m++) {
    count_one[m] = 0.0;
    for (i = 0; i < nvalues; i++) values_one[m][i] = 0.0;
  }

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms
  // wrap compute_ichunk in clearstep/addstep b/c it may invoke computes

  if (cchunk->computeflag) modify->clearstep_compute();

  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (cchunk->computeflag) modify->addstep_compute(ntimestep+nevery);

  // perform the computation for one sample
  // count # of atoms in each bin
  // accumulate results of attributes,computes,fixes,variables to local copy
  // sum within each chunk, only include atoms in fix group
  // compute/fix/variable may invoke computes so wrap with clear/add

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && ichunk[i] > 0)
      count_one[ichunk[i]-1]++;

  modify->clearstep_compute();

  m = 0;
  for (auto &val : values) {
    j = val.argindex;

    // V,F adds velocities,forces to values

    if (val.which == ArgInfo::V || val.which == ArgInfo::F) {
      double **attribute;
      if (val.which == ArgInfo::V) attribute = atom->v;
      else attribute = atom->f;

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && ichunk[i] > 0) {
          index = ichunk[i]-1;
          values_one[index][m] += attribute[i][j];
        }

    // DENSITY_NUMBER adds 1 to values

    } else if (val.which == ArgInfo::DENSITY_NUMBER) {

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && ichunk[i] > 0) {
          index = ichunk[i]-1;
          values_one[index][m] += 1.0;
        }

    // DENSITY_MASS or MASS adds mass to values

    } else if ((val.which == ArgInfo::DENSITY_MASS) || (val.which == ArgInfo::MASS)) {
      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;

      if (rmass) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && ichunk[i] > 0) {
            index = ichunk[i]-1;
            values_one[index][m] += rmass[i];
          }
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && ichunk[i] > 0) {
            index = ichunk[i]-1;
            values_one[index][m] += mass[type[i]];
          }
      }

    // TEMPERATURE adds KE to values
    // subtract and restore velocity bias if requested

    } else if (val.which == ArgInfo::TEMPERATURE) {

      if (biasflag) {
        if (tbias->invoked_scalar != ntimestep) tbias->compute_scalar();
        tbias->remove_bias_all();
      }

      double **v = atom->v;
      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;

      if (rmass) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && ichunk[i] > 0) {
            index = ichunk[i]-1;
            values_one[index][m] +=
              (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
          }
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && ichunk[i] > 0) {
            index = ichunk[i]-1;
            values_one[index][m] +=
              (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
              mass[type[i]];
          }
      }

      if (biasflag) tbias->restore_bias_all();

    // COMPUTE adds its scalar or vector component to values
    // invoke compute if not previously invoked

    } else if (val.which == ArgInfo::COMPUTE) {
      if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
        val.val.c->compute_peratom();
        val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
      }
      double *vector = val.val.c->vector_atom;
      double **array = val.val.c->array_atom;
      int jm1 = j - 1;

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && ichunk[i] > 0) {
          index = ichunk[i]-1;
          if (j == 0) values_one[index][m] += vector[i];
          else values_one[index][m] += array[i][jm1];
        }

    // FIX adds its scalar or vector component to values
    // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      double *vector = val.val.f->vector_atom;
      double **array = val.val.f->array_atom;
      int jm1 = j - 1;

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && ichunk[i] > 0) {
          index = ichunk[i]-1;
          if (j == 0) values_one[index][m] += vector[i];
          else values_one[index][m] += array[i][jm1];
        }

    // VARIABLE adds its per-atom quantities to values
    // evaluate atom-style variable

    } else if (val.which == ArgInfo::VARIABLE) {
      if (atom->nmax > maxvar) {
        maxvar = atom->nmax;
        memory->destroy(varatom);
        memory->create(varatom,maxvar,"ave/chunk:varatom");
      }

      input->variable->compute_atom(val.val.v,igroup,varatom,1,0);

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && ichunk[i] > 0) {
          index = ichunk[i]-1;
          values_one[index][m] += varatom[i];
        }
    }
    ++m;
  }

  // process the current sample
  // if normflag = ALL, accumulate values,count separately to many
  // if normflag = SAMPLE, one = value/count, accumulate one to many
  //   count is MPI summed here, value is MPI summed below across samples
  //   exception is TEMPERATURE: normalize by DOF
  //   exception is DENSITY_NUMBER:
  //     normalize by bin volume, not by atom count
  //   exception is DENSITY_MASS:
  //     scale by mv2d, normalize by bin volume, not by atom count
  //   exception is scaleflag = NOSCALE (norm = NONE):
  //     no normalize by atom count
  //     check last so other options can take precedence

  double mvv2e = force->mvv2e;
  double mv2d = force->mv2d;
  double boltz = force->boltz;

  if (normflag == ALL) {
    for (m = 0; m < nchunk; m++) {
      count_many[m] += count_one[m];
      for (j = 0; j < nvalues; j++)
        values_many[m][j] += values_one[m][j];
    }
  } else if (normflag == SAMPLE) {
    MPI_Allreduce(count_one,count_many,nchunk,MPI_DOUBLE,MPI_SUM,world);

    if (cchunk->chunk_volume_vec) {
      volflag = VECTOR;
      chunk_volume_vec = cchunk->chunk_volume_vec;
    } else {
      volflag = SCALAR;
      chunk_volume_scalar = cchunk->chunk_volume_scalar;
    }

    for (m = 0; m < nchunk; m++) {
      if (count_many[m] > 0.0)
        for (j = 0; j < nvalues; j++) {
          if (values[j].which == ArgInfo::TEMPERATURE) {
            values_many[m][j] += mvv2e*values_one[m][j] /
              ((cdof + adof*count_many[m]) * boltz);
          } else if (values[j].which == ArgInfo::DENSITY_NUMBER) {
            if (volflag == SCALAR) values_one[m][j] /= chunk_volume_scalar;
            else values_one[m][j] /= chunk_volume_vec[m];
            values_many[m][j] += values_one[m][j];
          } else if (values[j].which == ArgInfo::DENSITY_MASS) {
            if (volflag == SCALAR) values_one[m][j] /= chunk_volume_scalar;
            else values_one[m][j] /= chunk_volume_vec[m];
            values_many[m][j] += mv2d*values_one[m][j];
          } else if (scaleflag == NOSCALE) {
            values_many[m][j] += values_one[m][j];
          } else {
            values_many[m][j] += values_one[m][j]/count_many[m];
          }
        }
      count_sum[m] += count_many[m];
    }
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - ((bigint)nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // unlock compute chunk/atom at end of Nfreq epoch
  // do not unlock if ave = RUNNING or WINDOW

  if (nrepeat > 1 && ave == ONE) cchunk->unlock(this);

  // time average across samples
  // if normflag = ALL, final is total value / total count
  //   exception is TEMPERATURE: normalize by DOF for total count
  //   exception is DENSITY_NUMBER:
  //     normalize by final bin_volume and repeat, not by total count
  //   exception is DENSITY_MASS:
  //     scale by mv2d, normalize by bin volume and repeat, not by total count
  //   exception is scaleflag == NOSCALE:
  //     normalize by repeat, not by total count
  //     check last so other options can take precedence
  // if normflag = SAMPLE, final is sum of ave / repeat

  double repeat = nrepeat;

  if (normflag == ALL) {
    MPI_Allreduce(count_many,count_sum,nchunk,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nchunk*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);

    if (cchunk->chunk_volume_vec) {
      volflag = VECTOR;
      chunk_volume_vec = cchunk->chunk_volume_vec;
    } else {
      volflag = SCALAR;
      chunk_volume_scalar = cchunk->chunk_volume_scalar;
    }

    for (m = 0; m < nchunk; m++) {
      if (count_sum[m] > 0.0)
        for (j = 0; j < nvalues; j++) {
          if (values[j].which == ArgInfo::TEMPERATURE) {
            values_sum[m][j] *= mvv2e/((repeat*cdof + adof*count_sum[m])*boltz);
          } else if (values[j].which == ArgInfo::DENSITY_NUMBER) {
            if (volflag == SCALAR) values_sum[m][j] /= chunk_volume_scalar;
            else values_sum[m][j] /= chunk_volume_vec[m];
            values_sum[m][j] /= repeat;
          } else if (values[j].which == ArgInfo::DENSITY_MASS) {
            if (volflag == SCALAR) values_sum[m][j] /= chunk_volume_scalar;
            else values_sum[m][j] /= chunk_volume_vec[m];
            values_sum[m][j] *= mv2d/repeat;
          } else if (scaleflag == NOSCALE) {
            values_sum[m][j] /= repeat;
          } else {
            values_sum[m][j] /= count_sum[m];
          }
        }
      count_sum[m] /= repeat;
    }
  } else if (normflag == SAMPLE) {
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nchunk*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nchunk; m++) {
      for (j = 0; j < nvalues; j++) values_sum[m][j] /= repeat;
      count_sum[m] /= repeat;
    }
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, comine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (m = 0; m < nchunk; m++) {
      for (i = 0; i < nvalues; i++)
        values_total[m][i] = values_sum[m][i];
      count_total[m] = count_sum[m];
    }
    normcount = 1;

  } else if (ave == RUNNING) {
    for (m = 0; m < nchunk; m++) {
      for (i = 0; i < nvalues; i++)
        values_total[m][i] += values_sum[m][i];
      count_total[m] += count_sum[m];
    }
    normcount++;

  } else if (ave == WINDOW) {
    for (m = 0; m < nchunk; m++) {
      for (i = 0; i < nvalues; i++) {
        values_total[m][i] += values_sum[m][i];
        if (window_limit) values_total[m][i] -= values_list[iwindow][m][i];
        values_list[iwindow][m][i] = values_sum[m][i];
      }
      count_total[m] += count_sum[m];
      if (window_limit) count_total[m] -= count_list[iwindow][m];
      count_list[iwindow][m] = count_sum[m];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) normcount = nwindow;
    else normcount = iwindow;
  }

  // output result to file

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (overwrite) platform::fseek(fp,filepos);
    double count = 0.0;
    for (m = 0; m < nchunk; m++) count += count_total[m];
    fmt::print(fp,"{} {} {}\n",ntimestep,nchunk,count);

    int compress = cchunk->compress;
    int *chunkID = cchunk->chunkID;
    int ncoord = cchunk->ncoord;
    double **coord = cchunk->coord;

    if (!compress) {
      if (ncoord == 0) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp,"  %d %g",m+1,count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      } else if (ncoord == 1) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp,"  %d %g %g",m+1,coord[m][0],
                  count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      } else if (ncoord == 2) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp,"  %d %g %g %g",m+1,coord[m][0],coord[m][1],
                  count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      } else if (ncoord == 3) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp,"  %d %g %g %g %g",m+1,
                  coord[m][0],coord[m][1],coord[m][2],count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      }
    } else {
      if (ncoord == 0) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp,"  %d %d %g",m+1,chunkID[m],count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      } else if (ncoord == 1) {
        for (m = 0; m < nchunk; m++) {
          j = chunkID[m];
          fprintf(fp,"  %d %d %g %g",m+1,j,coord[j-1][0],
                  count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      } else if (ncoord == 2) {
        for (m = 0; m < nchunk; m++) {
          j = chunkID[m];
          fprintf(fp,"  %d %d %g %g %g",m+1,j,coord[j-1][0],coord[j-1][1],
                  count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      } else if (ncoord == 3) {
        for (m = 0; m < nchunk; m++) {
          j = chunkID[m];
          fprintf(fp,"  %d %d %g %g %g %g",m+1,j,coord[j-1][0],
                  coord[j-1][1],coord[j-1][2],count_total[m]/normcount);
          for (i = 0; i < nvalues; i++)
            fprintf(fp,format,values_total[m][i]/normcount);
          fprintf(fp,"\n");
        }
      }
    }
    if (ferror(fp))
      error->one(FLERR,"Error writing averaged chunk data");

    fflush(fp);

    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR,"Error while tuncating output: {}", utils::getsyserror());
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all per-chunk vectors
------------------------------------------------------------------------- */

void FixAveChunk::allocate()
{
  size_array_rows = nchunk;

  // reallocate chunk arrays if needed

  if (nchunk > maxchunk) {
    maxchunk = nchunk;
    memory->grow(count_one,nchunk,"ave/chunk:count_one");
    memory->grow(count_many,nchunk,"ave/chunk:count_many");
    memory->grow(count_sum,nchunk,"ave/chunk:count_sum");
    memory->grow(count_total,nchunk,"ave/chunk:count_total");

    memory->grow(values_one,nchunk,nvalues,"ave/chunk:values_one");
    memory->grow(values_many,nchunk,nvalues,"ave/chunk:values_many");
    memory->grow(values_sum,nchunk,nvalues,"ave/chunk:values_sum");
    memory->grow(values_total,nchunk,nvalues,"ave/chunk:values_total");

    // only allocate count and values list for ave = WINDOW

    if (ave == WINDOW) {
      memory->create(count_list,nwindow,nchunk,"ave/chunk:count_list");
      memory->create(values_list,nwindow,nchunk,nvalues,"ave/chunk:values_list");
    }

    // reinitialize regrown count/values total since they accumulate

    int i,m;
    for (m = 0; m < nchunk; m++) {
      for (i = 0; i < nvalues; i++) values_total[m][i] = 0.0;
      count_total[m] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   return I,J array value
   if I exceeds current nchunks, return 0.0 instead of generating an error
   columns 1 to colextra = chunkID + ncoord
   next column = count, remaining columns = Nvalues
------------------------------------------------------------------------- */

double FixAveChunk::compute_array(int i, int j)
{
  if (values_total == nullptr) return 0.0;
  if (i >= nchunk) return 0.0;
  if (j < colextra) {
    if (cchunk->compress) {
      if (j == 0) return (double) cchunk->chunkID[i];
      return cchunk->coord[i][j-1];
    } else return cchunk->coord[i][j];
  }
  j -= colextra + 1;
  if (!normcount) return 0.0;
  if (j < 0) return count_total[i]/normcount;
  return values_total[i][j]/normcount;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveChunk::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= ((bigint)nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ----------------------------------------------------------------------
   memory usage of varatom and bins
------------------------------------------------------------------------- */

double FixAveChunk::memory_usage()
{
  double bytes = (double)maxvar * sizeof(double);         // varatom
  bytes += (double)4*maxchunk * sizeof(double);           // count one,many,sum,total
  bytes += (double)nvalues*maxchunk * sizeof(double);     // values one,many,sum,total
  bytes += (double)nwindow*maxchunk * sizeof(double);          // count_list
  bytes += (double)nwindow*maxchunk*nvalues * sizeof(double);  // values_list
  return bytes;
}

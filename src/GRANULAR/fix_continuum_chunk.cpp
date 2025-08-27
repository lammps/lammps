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

/* ----------------------------------------------------------------------
   Contributing author: Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "fix_continuum_chunk.h"

#include "arg_info.h"
#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace NeighConst;
using namespace MathConst;

enum { OTHER, GRANULAR };
enum { DENSITY, VOLFRAC, MOMENTUM, VELOCITY, MGRAD, VGRAD, STRAINRATE, STRESS, STRESSKE, STRESSCON, IFD, FABRIC };
enum { BOUNDARY_NONE, BOUNDARY_FIX, BOUNDARY_ATOM, BOUNDARY_BOTH };
enum { SCALAR, VECTOR };
enum { SAMPLE, ALL };
enum { NOSCALE, ATOM };
enum { ONE, RUNNING, WINDOW };

static const char cite_continuum[] =
    "Coarse-graining procedure: doi:10.1007/s10035-010-0181-z\n\n"
    "@Article{Goldhirsch2010,\n"
    " author = {Goldhirsch, Isaac},\n"
    " title = {{Stress, stress asymmetry and couple stress: From discrete particles to continuous fields}},\n"
    " journal = {Granular Matter},\n"
    " year =    2010,\n"
    " volume =  12,\n"
    " number =  3,\n"
    " pages =   {239--252}\n"
    "}\n\n";

static const char cite_boundary[] =
    "Boundary corrections: doi:10.1007/s10035-012-0317-4\n\n"
    "@Article{Weinhart2012,\n"
    " author = {Weinhart, Thomas and Thornton, Anthony R. and Luding, Stefan and Bokhove, Onno},\n"
    " title = {{From discrete particles to continuum fields near a boundary}},\n"
    " journal = {Granular Matter},\n"
    " year =    2012,\n"
    " volume =  14,\n"
    " number =  2,\n"
    " pages =   {289--294}\n"
    "}\n\n";

inline double FixContinuumChunk::calc_w(double r) const
{
  return w_scale * exp(-(r * r) / (2 * w_sd_sq)) - w_offset;
}

inline double FixContinuumChunk::calc_dw(double r) const
{
  return -w_scale * exp(-(r * r) / (2 * w_sd_sq)) / w_sd_sq;
  // missing factor of x, added after called
}

inline double FixContinuumChunk::calc_w_int(double dr_dot_dr, double dr_dot_rij, double rij_dot_rij) const
{
  double sqrt_rij_dot_rij = sqrt(rij_dot_rij);
  double tmp = MY_SQRT2 * sqrt_rij_dot_rij * w_sd;
  double w_int = erf(dr_dot_rij / tmp) + erf((rij_dot_rij - dr_dot_rij) / tmp);
  w_int *= exp((dr_dot_rij * dr_dot_rij - dr_dot_dr * rij_dot_rij) / (tmp * tmp));
  w_int *= sqrt(0.5 * MY_PI) * w_sd / sqrt_rij_dot_rij;
  return w_scale * w_int - w_offset;
}

/* ---------------------------------------------------------------------- */

FixContinuumChunk::FixContinuumChunk(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), nrepeat(0), fp(nullptr), idchunk(nullptr), varatom(nullptr),
    count_one(nullptr), count_many(nullptr), count_sum(nullptr), values_one(nullptr),
    values_many(nullptr), values_sum(nullptr), count_total(nullptr), count_list(nullptr),
    values_total(nullptr), values_list(nullptr), density_one(nullptr), density_sum_now(nullptr),
    momentum_one(nullptr), momentum_sum_now(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix continuum/chunk", error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[5], false, lmp);

  idchunk = utils::strdup(arg[6]);
  w_cut = utils::numeric(FLERR, arg[7], false, lmp);
  w_sd = utils::numeric(FLERR, arg[8], false, lmp);

  global_freq = nfreq;
  no_change_box = 1;
  time_depend = 1;
  dim = domain->dimension;

  need_momentum = 0;
  need_density = 0;
  need_vgrad = 0;
  calculate_pair = 0;
  calculate_grad = 0;
  index_density = -1;
  for (int a = 0; a < 3; a++) {
    index_momentum[a] = -1;
    for (int b = 0; b < 3; b++) {
      index_vgrad[a][b] = -1;
    }
  }

  char *mygroup = arg[1];

  // parse values until one isn't recognized

  int iarg = 9;
  values.clear();
  labels.clear();
  while (iarg < narg) {
    if (strcmp(arg[iarg], "density") == 0) {
      values.push_back(std::make_pair(DENSITY, -1));
      labels.push_back("density");
      index_density = values.size() - 1;
    } else if (strcmp(arg[iarg], "volume/fraction") == 0) {
      values.push_back(std::make_pair(VOLFRAC, -1));
      labels.push_back("volume/fraction");
    } else if (utils::strmatch(arg[iarg], "^momentum/.$")) {
      add_vector_component(arg[iarg], MOMENTUM);
    } else if (utils::strmatch(arg[iarg], "^velocity/.$")) {
      add_vector_component(arg[iarg], VELOCITY);
      need_density = 1;
      need_momentum = 1;
    } else if (utils::strmatch(arg[iarg], "^momentum/grad/")) {
      add_tensor_component(arg[iarg], MGRAD);
      printf("momentum grad\n");
      need_momentum = 1;
      calculate_grad = 1;
    } else if (utils::strmatch(arg[iarg], "^velocity/grad/")) {
      add_tensor_component(arg[iarg], VGRAD);
      printf("velocity grad\n");
      need_density = 1;
      need_momentum = 1;
      calculate_grad = 1;
    } else if (utils::strmatch(arg[iarg], "^strain/rate/")) {
      add_tensor_component(arg[iarg], STRAINRATE);
      need_density = 1;
      need_momentum = 1;
      need_vgrad = 1;
      calculate_grad = 1;
    } else if (utils::strmatch(arg[iarg], "^stress/.$") ||
               utils::strmatch(arg[iarg], "^stress/..$")) {
      add_tensor_component(arg[iarg], STRESS);
      calculate_pair = 1;
    } else if (utils::strmatch(arg[iarg], "^stress/ke/")) {
      add_tensor_component(arg[iarg], STRESSKE);
    } else if (utils::strmatch(arg[iarg], "^stress/contacts/")) {
      add_tensor_component(arg[iarg], STRESSCON);
      calculate_pair = 1;
    } else if (utils::strmatch(arg[iarg], "^boundary/force")) {
      add_vector_component(arg[iarg], IFD);
      calculate_pair = 1;
    } else if (utils::strmatch(arg[iarg], "^fabric/")) {
      add_tensor_component(arg[iarg], FABRIC);
      calculate_pair = 1;
    } else {
      break;
    }
    iarg++;
  }

  // Add any necessary intermediate values, won't be printed
  nskip = 0;
  if (need_density && index_density == -1) {
    values.push_back(std::make_pair(DENSITY, -1));
    index_density = values.size() - 1;
    nskip += 1;
  }

  if (need_momentum) {
    for (int a = 0; a < dim; a++) {
      if (index_momentum[a] == -1) {
        values.push_back(std::make_pair(MOMENTUM, a));
        index_momentum[a] = values.size() - 1;
        nskip += 1;
      }
    }
  }

  if (need_vgrad) {
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        if (dim == 2 && (b == 2 || a == 2)) continue;
        if (index_vgrad[a][b] == -1) {
          values.push_back(std::make_pair(VGRAD, a * 3 + b));
          index_vgrad[a][b] = values.size() - 1;
          nskip += 1;
        }
      }
    }
  }

  nvalues = values.size();
  if (nvalues == 0) error->all(FLERR, "No values in fix continuum/chunk command");

  // optional args

  boundaryflag = BOUNDARY_NONE;
  boundary_group_flag = 0;
  ave = ONE;
  nwindow = 0;
  overwrite = 0;
  format_user = nullptr;
  format = (char *) " %g";
  char *title1 = nullptr;
  char *title2 = nullptr;
  char *title3 = nullptr;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "boundary/fix") == 0) {
      if (boundaryflag == BOUNDARY_ATOM)
        boundaryflag = BOUNDARY_BOTH;
      else
        boundaryflag = BOUNDARY_FIX;
      iarg += 1;
    } else if (strcmp(arg[iarg], "boundary/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix continuum/chunk ") + arg[iarg], error);
      if (boundaryflag == BOUNDARY_FIX)
        boundaryflag = BOUNDARY_BOTH;
      else
        boundaryflag = BOUNDARY_ATOM;
      boundary_group_flag = 1;
      boundary_groupbit = group->get_bitmask_by_id(FLERR, arg[iarg + 1], "fix continuum/chunk");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ave") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix continuum/chunk ave", error);
      if (strcmp(arg[iarg + 1], "one") == 0) ave = ONE;
      else if (strcmp(arg[iarg + 1], "running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg + 1], "window") == 0) ave = WINDOW;
      else error->all(FLERR, "Unknown fix continuum/chunk ave mode: {}", arg[iarg + 1]);
      if (ave == WINDOW) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix continuum/chunk ave window", error);
        nwindow = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
        if (nwindow <= 0) error->all(FLERR, "Illegal fix continuum/chunk number of windows: {}", nwindow);
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;


    } else if ((strcmp(arg[iarg], "file") == 0) || (strcmp(arg[iarg], "append") == 0)) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix continuum/chunk ") + arg[iarg], error);
      if (comm->me == 0) {
        if (strcmp(arg[iarg], "file") == 0) fp = fopen(arg[iarg + 1], "w");
        else fp = fopen(arg[iarg + 1], "a");
        if (fp == nullptr)
          error->one(FLERR, "Cannot open fix continuum/chunk file {}: {}",
                     arg[iarg + 1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg], "format") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix continuum/chunk format", error);
      delete[] format_user;
      format_user = utils::strdup(arg[iarg + 1]);
      format = format_user;
      iarg += 2;
    } else if (strcmp(arg[iarg], "title1") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix continuum/chunk title1", error);
      delete[] title1;
      title1 = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "title2") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix continuum/chunk title2", error);
      delete[] title2;
      title2 = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "title3") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix continuum/chunk title3", error);
      delete[] title3;
      title3 = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else error->all(FLERR, "Unknown fix continuum/chunk keyword: {}", arg[iarg]);
  }

  // setup and error check

  if (nevery <= 0) error->all(FLERR, "Illegal fix continuum/chunk nevery value: {}", nevery);
  if (nrepeat <= 0) error->all(FLERR, "Illegal fix continuum/chunk nrepeat value: {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR, "Illegal fix continuum/chunk nfreq value: {}", nfreq);
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR, "Inconsistent fix continuum/chunk nevery/nrepeat/nfreq values");
  if (ave != RUNNING && overwrite)
    error->all(FLERR, "Fix continuum/chunk overwrite keyword requires ave running setting");
  if (!boundaryflag)
    for (auto &val : values)
      if (val.first == IFD)
        error->all(FLERR, "Must specify type of boundary corrections, atom and/or fix, to compute boundary/force");

  // increment lock counter in compute chunk/atom
  // only if nrepeat > 1 or ave = RUNNING/WINDOW,
  //   so that locking spans multiple timesteps

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR, "Chunk/atom compute {} does not exist or is "
               "incorrect style for fix continuum/chunk", idchunk);

  int which = cchunk->get_which();
  if (which == ArgInfo::BIN1D) {
    bin_dim = 1;
  } else if (which == ArgInfo::BIN2D) {
    bin_dim = 2;
  } else if (which == ArgInfo::BIN3D) {
    bin_dim = 3;
  } else {
    error->all(FLERR, "Can only use bin chunk/atom styles with fix continuum/chunk");
  }

  w_cut_sq = w_cut * w_cut;
  w_sd_sq = w_sd * w_sd;

  // Normalization factor for truncated Gaussian
  double exp_cut = exp(-w_cut_sq / (2.0 * w_sd_sq));
  if (bin_dim == 1) {
    w_scale = sqrt(2.0 * MY_PI) * w_sd * erf(w_cut / (MY_SQRT2 * w_sd));
    w_scale -= 2.0 * w_cut * exp_cut;
    w_scale = 1.0 / w_scale;
  } else if (bin_dim == 2) {
    w_scale = -0.5 * w_cut_sq * exp_cut + w_sd_sq * (1.0 - exp_cut);
    w_scale = 1.0 / (2.0 * MY_PI * w_scale);
  } else {
    w_scale = -THIRD * w_cut * exp_cut * (w_cut_sq + 3.0 * w_sd_sq);
    w_scale += sqrt(0.5 * MY_PI) * w_sd_sq * w_sd * erf(w_cut / (MY_SQRT2 * w_sd));
    w_scale = 1.0 / (4.0 * MY_PI * w_scale);
  }
  w_offset = w_scale * exp(-w_cut_sq / (2 * w_sd_sq));

  if ((nrepeat > 1) || (ave == RUNNING) || (ave == WINDOW)) cchunk->lockcount++;
  lockforever = 0;

  // print file comment lines

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp, "%s\n",title1);
    else fprintf(fp, "# Chunk-averaged data for fix %s and group %s\n", id, mygroup);
    if (title2) fprintf(fp, "%s\n",title2);
    else fprintf(fp, "# Timestep Number-of-chunks Total-count\n");
    if (title3) fprintf(fp, "%s\n",title3);
    else {
      int compress = cchunk->compress;
      int ncoord = cchunk->ncoord;
      if (!compress) {
        if (ncoord == 0) fprintf(fp, "# Chunk Ncount");
        else if (ncoord == 1) fprintf(fp, "# Chunk Coord1 Ncount");
        else if (ncoord == 2) fprintf(fp, "# Chunk Coord1 Coord2 Ncount");
        else if (ncoord == 3)
          fprintf(fp, "# Chunk Coord1 Coord2 Coord3 Ncount");
      } else {
        if (ncoord == 0) fprintf(fp, "# Chunk OrigID Ncount");
        else if (ncoord == 1) fprintf(fp, "# Chunk OrigID Coord1 Ncount");
        else if (ncoord == 2) fprintf(fp, "# Chunk OrigID Coord1 Coord2 Ncount");
        else if (ncoord == 3)
          fprintf(fp, "# Chunk OrigID Coord1 Coord2 Coord3 Ncount");
      }
      for (int i = 0; i < (nvalues - nskip); i++)
        fprintf(fp, " %s", labels[i].c_str());
      fprintf(fp, "\n");
    }
    if (ferror(fp))
      error->one(FLERR, "Error writing file header");

    filepos = platform::ftell(fp);
  }

  delete[] title1;
  delete[] title2;
  delete[] title3;

  // this fix produces a global array
  // size_array_rows is variable and set by allocate()

  int compress = cchunk->compress;
  int ncoord = cchunk->ncoord;
  colextra = compress + ncoord;

  array_flag = 1;
  size_array_cols = colextra + 1 + nvalues - nskip;
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

  if (lmp->citeme) {
    lmp->citeme->add(cite_continuum);
    if (boundaryflag)
      lmp->citeme->add(cite_boundary);
  }
}

/* ---------------------------------------------------------------------- */

FixContinuumChunk::~FixContinuumChunk()
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

  memory->destroy(density_one);
  memory->destroy(density_sum_now);
  memory->destroy(momentum_one);
  memory->destroy(momentum_sum_now);

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

int FixContinuumChunk::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContinuumChunk::init()
{
  // set indices and check validity of all computes,fixes,variables
  // check that fix frequency is acceptable

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR, "Chunk/atom compute {} does not exist or is "
               "incorrect style for fix continuum/chunk",idchunk);

  if (boundaryflag == BOUNDARY_FIX || boundaryflag == BOUNDARY_BOTH) {
    auto wall_fixes = modify->get_fix_by_style("wall/gran");
    if (wall_fixes.size() == 0)
      error->all(FLERR, "Could not find any instances of fix wall/gran for boundary corrections");
    for (auto fix : wall_fixes)
      if (!fix->peratom_flag)
        error->all(FLERR, "Must use contacts keyword in fix wall/gran {} for boundary corrections", fix->id);
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }

  // to calculate the stress, need pair->single()

  if (calculate_pair) {
    if (force->pair == nullptr) error->all(FLERR, "No pair style is defined for fix continuum chunk stress calculation");
    if (force->pair->single_enable == 0)
      error->all(FLERR, "Pair style does not support stress calculation");

    // Find if granular or gran, need to include tangential forces

    pstyle = OTHER;
    if (force->pair_match("^granular", 0) || force->pair_match("^gran/", 0)) pstyle = GRANULAR;

    // need an occasional full neighbor list
    // set size to same value as request made by force->pair
    // should be able to derive full list from pair list

    auto *pairrequest = neighbor->find_request(force->pair);
    if (pairrequest && pairrequest->get_size())
      neighbor->add_request(this, NeighConst::REQ_SIZE |   NeighConst::REQ_OCCASIONAL | REQ_FULL);
    else
      neighbor->add_request(this, NeighConst::REQ_OCCASIONAL | REQ_FULL);
  }
}

/* ---------------------------------------------------------------------- */

void FixContinuumChunk::init_list(int /*id*/, NeighList *ptr)
{
  if (calculate_pair) list = ptr;
}

/* ----------------------------------------------------------------------
   only does averaging if nvalid = current timestep
   do not call setup_chunks(), even though fix continuum/chunk called setup_bins()
   b/c could cause nchunk to change if Nfreq epoch crosses 2 runs
   does mean that if change_box is used between runs to change box size,
     that nchunk may not track it
------------------------------------------------------------------------- */

void FixContinuumChunk::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixContinuumChunk::end_of_step()
{
  int i, j,m;

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
      modify->addstep_compute(ntimestep + nevery);
      modify->addstep_compute(ntimestep + nfreq);
    }
    allocate();
    if (nrepeat > 1 && ave == ONE)
      cchunk->lock(this, ntimestep, ntimestep + ((bigint)nrepeat - 1) * nevery);
    else if ((ave == RUNNING || ave == WINDOW) && !lockforever) {
      cchunk->lock(this,update->ntimestep,-1);
      lockforever = 1;
    }
    for (m = 0; m < nchunk; m++) {
      count_many[m] = count_sum[m] = 0.0;
      for (i = 0; i < nvalues; i++) values_many[m][i] = 0.0;
    }
  }

  // invoke setup_chunks() on each sampling step and grab relevant values
  // geometry could change, e.g. for NPT simulation

  cchunk->setup_chunks();
  int ncoord = cchunk->ncoord;
  double **coord = cchunk->coord;
  int reducedflag = cchunk->get_reducedflag();
  int *cdim = cchunk->get_dim();
  double *delta = cchunk->get_delta();

  double width;
  for (m = 0; m < ncoord; m++) {
    width = delta[m];
    if (reducedflag)
      width *= domain->prd[cdim[m]];
    if (0.5 * width < w_cut)
      error->all(FLERR, "Chunk half width {} smaller than specified cutoff {}", 0.5 * width, w_cut);
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

  if (cchunk->computeflag) modify->addstep_compute(ntimestep + nevery);

  // perform the computation for one sample
  // count # of atoms in each bin
  // accumulate results of attributes to local copy
  // sum within each chunk, only include atoms in fix group
  // compute/fix/variable may invoke computes so wrap with clear/add

  int a, b, itype, style, component, field_index;
  double w, wc, mi, voli, r, rsq_bin, rsq_pair, rsq_wall, rbin_dot_r, f_norm, w_int_tmp, wc_int_tmp;
  double xbin[3], xcont[3], dx_bin[3], dx_pair[3], dx_cont[3], f_pair[3], f_wall[3], dx_wall[3];

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int jj, jnum;
  int *jlist, *numneigh, **firstneigh;

  double **array_atom_fix;

  if (calculate_pair) {
    neighbor->build_one(list);
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  auto wall_fixes = modify->get_fix_by_style("wall/gran");

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && ichunk[i] > 0)
      count_one[ichunk[i]-1]++;

  modify->clearstep_compute();

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && ichunk[i] > 0) {
      m = ichunk[i] - 1;

      // x[i] is default so won't contribute unless binned in that coord
      MathExtra::copy3(x[i], xbin);
      for (a = 0; a < ncoord; a++) {
        xbin[cdim[a]] = coord[m][a];
        if (reducedflag) {
          xbin[cdim[a]] *= domain->prd[cdim[a]];
          xbin[cdim[a]] += domain->boxlo[cdim[a]];
        }
      }

      itype = type[i];
      if (rmass) mi = rmass[i];
      else mi = mass[itype];
      voli = MY_PI * radius[i] * radius[i];
      if (dim == 3)
        voli *= 4.0 * THIRD * radius[i];

      MathExtra::sub3(x[i], xbin, dx_bin);
      rsq_bin = MathExtra::lensq3(dx_bin);

      if (rsq_bin > w_cut_sq) continue;
      w = calc_w(sqrt(rsq_bin));

      // contributions from single atoms (excluding boundary)

      if (!(boundary_group_flag && (mask[i] & boundary_groupbit))) {
        field_index = 0;
        for (auto &val : values) {
          style = val.first;
          component = val.second;

          a = component % 3;
          b = (component - a) / 3;

          if (style == DENSITY) {
            values_one[m][field_index] += mi * w;
          } else if (style == VOLFRAC) {
            values_one[m][field_index] += voli * w;
          } else if (style == MOMENTUM) {
            values_one[m][field_index] += mi * v[i][component] * w;
          } else if (style == STRESS || style == STRESSKE) {
            values_one[m][field_index] -= mi * v[i][a] * v[i][b] * w;
          }

          // Boundary corrections from Weinhart et al. 2012
          if (boundaryflag && (style == STRESS || style == STRESSCON)) {
            for (auto wall_fix : wall_fixes) {
              array_atom_fix = wall_fix->array_atom;

              // Skip if not in contact
              if (array_atom_fix[i][0] != 1.0) continue;
              f_wall[0] = array_atom_fix[i][1];
              f_wall[1] = array_atom_fix[i][2];
              f_wall[2] = array_atom_fix[i][3];
              dx_wall[0] = x[i][0] - array_atom_fix[i][4];
              dx_wall[1] = x[i][1] - array_atom_fix[i][5];
              dx_wall[2] = x[i][2] - array_atom_fix[i][6];

              rsq_wall = MathExtra::lensq3(dx_wall);
              rbin_dot_r = MathExtra::dot3(dx_bin, dx_pair);
              w_int_tmp = calc_w_int(rsq_bin, rbin_dot_r, rsq_wall);

              values_one[m][field_index] -= f_wall[a] * dx_wall[b] * w_int_tmp;
            }
          }

          field_index++;
        }
      }

      // contributions from pairs of atoms

      if (calculate_pair) {
        jlist = firstneigh[i];
        jnum = numneigh[i];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          if (!mask[j] & groupbit) continue;

          MathExtra::sub3(x[i], x[j], dx_pair);
          rsq_pair = MathExtra::lensq3(dx_pair);
          pair->single(i, j, itype, type[j], rsq_pair, 1.0, 1.0, f_norm);

          MathExtra::scale3(f_norm / sqrt(rsq_pair), dx_pair, f_pair);
          if (pstyle == GRANULAR) {
            // add tangential forces
            f_pair[0] += force->pair->svector[0];
            f_pair[1] += force->pair->svector[1];
            f_pair[2] += force->pair->svector[2];
          }

          if (MathExtra::lensq3(f_pair) == 0.0)
            continue;

          rbin_dot_r = MathExtra::dot3(dx_bin, dx_pair);
          w_int_tmp = calc_w_int(rsq_bin, rbin_dot_r, rsq_pair);
          if (boundary_group_flag) {
            // Skip if both are boundary particles
            if ((mask[i] & boundary_groupbit) && (mask[j] & boundary_groupbit))
              continue;

            MathExtra::scale3(0.5, dx_pair, xcont);
            MathExtra::sub3(x[i], xcont, xcont);

            // Check which one is a boundary
            if (mask[i] & boundary_groupbit) {
              MathExtra::sub3(x[i], xcont, dx_cont);
            } else {
              MathExtra::sub3(xcont, x[j], dx_cont);
            }

            // Calculate alternate kernel metrics (could be more selective if only one needed)
            rsq_pair = MathExtra::lensq3(dx_cont); // reuse temp variables
            rbin_dot_r = MathExtra::dot3(dx_bin, dx_cont);
            wc_int_tmp = calc_w_int(rsq_bin, rbin_dot_r, rsq_pair);
            wc = calc_w(MathExtra::len3(dx_cont));
          }

          field_index = 0;
          for (auto &val : values) {
            style = val.first;
            component = val.second;

            a = component % 3;
            b = (component - a) / 3;

            if (style == STRESS || style == STRESSCON) {
              if (boundary_group_flag && (mask[i] & boundary_groupbit)) {
                values_one[m][field_index] -= f_pair[a] * dx_cont[b] * wc_int_tmp;
              } else {
                values_one[m][field_index] -= f_pair[a] * dx_pair[b] * w_int_tmp;
              }
            } else if (style == IFD) {
              if (!(boundary_group_flag && (mask[i] & boundary_groupbit))) continue;
              values_one[m][field_index] -= f_pair[a] * dx_cont[b] * wc;
            } else if (style == FABRIC) {
              if (boundary_group_flag && (mask[i] & boundary_groupbit)) continue;
              values_one[m][field_index] += voli * dx_pair[a] * dx_pair[b] * w_int_tmp / rsq_pair;
            }

            field_index++;
          }
        }
      }
    }
  }

  if (calculate_grad) {

    // Copy intermediate values and sum across processors (will repeat later)
    for (m = 0; m < nchunk; m++) {
      density_one[m] = values_one[m][index_density];
      for (a = 0; a < 3; a++) {
        if (a < dim) {
          momentum_one[m][a] = values_one[m][index_momentum[a]];
        } else {
          momentum_one[m][a] = 0.0;
        }
      }
    }

    MPI_Allreduce(&density_one[0], &density_sum_now[0], nchunk,
                MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&momentum_one[0][0], &momentum_sum_now[0][0], nchunk * 3,
                MPI_DOUBLE, MPI_SUM, world);

    double dw, vbin;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && ichunk[i] > 0) {
        if (boundary_group_flag && (mask[i] & boundary_groupbit)) continue;

        itype = type[i];
        if (rmass) mi = rmass[i];
        else mi = mass[itype];
        voli = MY_PI * radius[i] * radius[i];
        if (dim == 3)
        voli *= 4.0 * THIRD * radius[i];

        m = ichunk[i] - 1;

        MathExtra::copy3(x[i], xbin);
        for (a = 0; a < ncoord; a++) {
          xbin[cdim[a]] = coord[m][a];
          if (reducedflag) {
            xbin[cdim[a]] *= domain->prd[cdim[a]];
            xbin[cdim[a]] += domain->boxlo[cdim[a]];
          }
        }

        itype = type[i];
        if (rmass) mi = rmass[i];
        else mi = mass[itype];

        MathExtra::sub3(x[i], xbin, dx_bin);
        rsq_bin = MathExtra::lensq3(dx_bin);

        if (rsq_bin > w_cut_sq) continue;
        dw = calc_dw(rsq_bin); // sans dx factor

        field_index = 0;
        for (auto &val : values) {
          style = val.first;
          component = val.second;

          a = component % 3;
          b = (component - a) / 3;

          if (style == MGRAD) {
            values_one[m][field_index] += voli * (momentum_sum_now[m][a] - mi * v[i][a]) * dx_bin[b] * dw;
          } else if (style == VGRAD) {
            if (density_sum_now[m] != 0.0) {
              vbin = momentum_sum_now[m][a] / density_sum_now[m];
            } else{
              vbin = 0.0;
            }

            values_one[m][field_index] += voli * (vbin - v[i][a]) * dx_bin[b] * dw;
          }

          field_index++;
        }
      }
    }
  }

  // process the current sample, one = value/count, accumulate one to many

  MPI_Allreduce(count_one, count_many, nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (m = 0; m < nchunk; m++) {
    for (j = 0; j < nvalues; j++)
      values_many[m][j] += values_one[m][j];
    count_sum[m] += count_many[m];
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
  nvalid = ntimestep + nfreq - ((bigint) nrepeat - 1) * nevery;
  modify->addstep_compute(nvalid);

  // unlock compute chunk/atom at end of Nfreq epoch
  // do not unlock if ave = RUNNING or WINDOW

  if (nrepeat > 1 && ave == ONE) cchunk->unlock(this);

  // time average across samples, final is sum of ave / repeat

  double repeat = nrepeat;

  MPI_Allreduce(&values_many[0][0], &values_sum[0][0], nchunk * nvalues,
                MPI_DOUBLE, MPI_SUM, world);
  for (m = 0; m < nchunk; m++) {
    for (j = 0; j < nvalues; j++) values_sum[m][j] /= repeat;
    count_sum[m] /= repeat;
  }


  // Calculate trivially derived values

  double dtemp, mtemp, mtemp2;
  for (m = 0; m < nchunk; m++) {
    field_index = 0;
    for (auto &val : values) {
      style = val.first;
      component = val.second;

      a = component % 3;
      b = (component - a) / 3;

      if (style == VELOCITY) {
        dtemp = values_sum[m][index_density];
        mtemp = values_sum[m][index_momentum[component]];
        if (dtemp != 0.0)
          values_sum[m][field_index] = mtemp / dtemp;
      } else if (style == STRAINRATE) {
        mtemp = values_sum[m][index_vgrad[a][b]];
        mtemp2 = values_sum[m][index_vgrad[b][a]];
        values_sum[m][field_index] = 0.5 * (mtemp + mtemp2);
      }

      field_index++;
    }
  }

  // Normalize by any unused dimensions

  if (bin_dim != dim) {
    int unused_dim[3] = {1, 1, 1};
    for (a = 0; a < ncoord; a++)
      unused_dim[cdim[a]] = 0;

    for (a = 0; a < dim; a++)
      if (unused_dim[a])
        for (m = 0; m < nchunk; m++)
          for (int n = 0; n < nvalues; n++)
            values_sum[m][n] /= domain->prd[a];
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

    if (overwrite) (void) platform::fseek(fp,filepos);
    double count = 0.0;
    for (m = 0; m < nchunk; m++) count += count_total[m];
    fmt::print(fp, "{} {} {}\n", ntimestep, nchunk, count);

    int compress = cchunk->compress;
    int *chunkID = cchunk->chunkID;
    int ncoord = cchunk->ncoord;
    double **coord = cchunk->coord;

    if (!compress) {
      if (ncoord == 0) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp, "  %d %g", m + 1, count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      } else if (ncoord == 1) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp, "  %d %g %g", m + 1, coord[m][0],
                  count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      } else if (ncoord == 2) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp, "  %d %g %g %g", m + 1, coord[m][0], coord[m][1],
                  count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      } else if (ncoord == 3) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp, "  %d %g %g %g %g", m + 1,
                  coord[m][0], coord[m][1], coord[m][2], count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      }
    } else {
      if (ncoord == 0) {
        for (m = 0; m < nchunk; m++) {
          fprintf(fp, "  %d %d %g", m + 1, chunkID[m], count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      } else if (ncoord == 1) {
        for (m = 0; m < nchunk; m++) {
          j = chunkID[m];
          fprintf(fp, "  %d %d %g %g", m + 1, j, coord[j - 1][0],
                  count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      } else if (ncoord == 2) {
        for (m = 0; m < nchunk; m++) {
          j = chunkID[m];
          fprintf(fp, "  %d %d %g %g %g", m + 1, j, coord[j - 1][0], coord[j - 1][1],
                  count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      } else if (ncoord == 3) {
        for (m = 0; m < nchunk; m++) {
          j = chunkID[m];
          fprintf(fp, "  %d %d %g %g %g %g", m + 1, j, coord[j - 1][0],
                  coord[j - 1][1], coord[j - 1][2], count_total[m] / normcount);
          for (i = 0; i < (nvalues - nskip); i++)
            fprintf(fp,format,values_total[m][i] / normcount);
          fprintf(fp, "\n");
        }
      }
    }
    if (ferror(fp))
      error->one(FLERR, "Error writing averaged chunk data");

    fflush(fp);

    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR, "Error while tuncating output: {}", utils::getsyserror());
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all per-chunk vectors
------------------------------------------------------------------------- */

void FixContinuumChunk::allocate()
{
  size_array_rows = nchunk;

  // reallocate chunk arrays if needed

  if (nchunk > maxchunk) {
    maxchunk = nchunk;
    memory->grow(count_one, nchunk, "continuum/chunk:count_one");
    memory->grow(count_many, nchunk, "continuum/chunk:count_many");
    memory->grow(count_sum, nchunk, "continuum/chunk:count_sum");
    memory->grow(count_total, nchunk, "continuum/chunk:count_total");

    memory->grow(values_one, nchunk, nvalues, "continuum/chunk:values_one");
    memory->grow(values_many, nchunk, nvalues, "continuum/chunk:values_many");
    memory->grow(values_sum, nchunk, nvalues, "continuum/chunk:values_sum");
    memory->grow(values_total, nchunk, nvalues, "continuum/chunk:values_total");

    if (calculate_grad) {
      memory->grow(momentum_one, nchunk, 3, "continuum/chunk:momentum_one");
      memory->grow(momentum_sum_now, nchunk, 3, "continuum/chunk:momentum_sum_now");
      memory->grow(density_one, nchunk, "continuum/chunk:density_one");
      memory->grow(density_sum_now, nchunk, "continuum/chunk:density_sum_now");
    }

    // only allocate count and values list for ave = WINDOW

    if (ave == WINDOW) {
      memory->create(count_list, nwindow, nchunk, "continuum/chunk:count_list");
      memory->create(values_list, nwindow, nchunk, nvalues, "continuum/chunk:values_list");
    }

    // reinitialize regrown count/values total since they accumulate

    int i, m;
    for (m = 0; m < nchunk; m++) {
      for (i = 0; i < nvalues; i++) values_total[m][i] = 0.0;
      count_total[m] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   return I, j array value
   if I exceeds current nchunks, return 0.0 instead of generating an error
   columns 1 to colextra = chunkID + ncoord
   next column = count, remaining columns = Nvalues
------------------------------------------------------------------------- */

double FixContinuumChunk::compute_array(int i, int j)
{
  if (values_total == nullptr) return 0.0;
  if (i >= nchunk) return 0.0;
  if (j < colextra) {
    if (cchunk->compress) {
      if (j == 0) return (double) cchunk->chunkID[i];
      return cchunk->coord[i][j - 1];
    } else return cchunk->coord[i][j];
  }
  j -= colextra + 1;
  if (!normcount) return 0.0;
  if (j < 0) return count_total[i] / normcount;
  return values_total[i][j] / normcount;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixContinuumChunk::nextvalid()
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

double FixContinuumChunk::memory_usage()
{
  double bytes = (double)maxvar * sizeof(double);         // varatom
  bytes += (double)4 * maxchunk * sizeof(double);           // count one,many,sum,total
  bytes += (double)nvalues * maxchunk * sizeof(double);     // values one,many,sum,total
  bytes += (double)nwindow * maxchunk * sizeof(double);          // count_list
  bytes += (double)nwindow * maxchunk*nvalues * sizeof(double);  // values_list
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixContinuumChunk::add_tensor_component(char *option, int variable)
{
  if (((std::string) option).back() == '*') {
    std::vector<std::string> suffices = {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
    std::string trimmed_option = std::string(option);
    trimmed_option = trimmed_option.substr(0, trimmed_option.length() - 1);
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        if (dim == 2 && (b == 2 || a == 2)) continue;
        values.push_back(std::make_pair(variable, a * 3 + b));
        labels.push_back(trimmed_option + suffices[a * 3 + b]);
        if (variable == VGRAD)
          index_vgrad[a][b] = values.size() - 1;
      }
    }
  } else {
    int index = -1;
    int dim_error = 0;

    if (utils::strmatch(option, "xx$")) {
      index = 0;
    } else if (utils::strmatch(option, "xy$")) {
      index = 1;
    } else if (utils::strmatch(option, "xz$")) {
      index = 2;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "yx$")) {
      index = 3;
    } else if (utils::strmatch(option, "yy$")) {
      index = 4;
    } else if (utils::strmatch(option, "yz$")) {
      index = 5;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zx$")) {
      index = 6;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zy$")) {
      index = 7;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zz$")) {
      index = 8;
      if (dim == 2) dim_error = 1;
    } else {
      error->all(FLERR, "Invalid fix continuum/chunk property {}", option);
    }

    if (dim_error)
      error->all(FLERR, "Invalid fix continuum/chunk property {} in 2D", option);

    values.push_back(std::make_pair(variable, index));
    labels.push_back(option);
    if (variable == VGRAD) {
      int a = index % 3;
      int b = (index - a) / 3;
      index_vgrad[a][b] = values.size() - 1;
    }
  }

  return;
}

/* ---------------------------------------------------------------------- */

void FixContinuumChunk::add_vector_component(char *option, int variable)
{
  if (((std::string) option).back() == '*') {
    std::vector<std::string> suffices = {"x", "y", "z"};
    std::string trimmed_option = std::string(option);
    trimmed_option = trimmed_option.substr(0, trimmed_option.length() - 1);
    for (int a = 0; a < dim; a++) {
      values.push_back(std::make_pair(variable, a));
      labels.push_back(trimmed_option + suffices[a]);
      if (variable == MOMENTUM)
        index_momentum[a] = values.size() - 1;
    }
  } else {
    int index = -1;
    if (utils::strmatch(option, "x$")) {
      index = 0;
    } else if (utils::strmatch(option, "y$")) {
      index = 1;
    } else if (utils::strmatch(option, "z$")) {
      if (dim == 2)
        error->all(FLERR, "Invalid fix continuum/chunk property {} in 2D", option);
      index = 2;
    } else {
      error->all(FLERR, "Invalid fix continuum/chunk property {}", option);
    }

    values.push_back(std::make_pair(variable, index));
    labels.push_back(option);
    if (variable == MOMENTUM)
        index_momentum[index] = values.size() - 1;
  }
}

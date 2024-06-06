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

#include "fix_ave_histo.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { ONE, RUNNING };
enum { SCALAR, VECTOR, WINDOW };
enum { DEFAULT, GLOBAL, PERATOM, LOCAL };
enum { IGNORE, END, EXTRA };

static constexpr double BIG = 1.0e20;
/* ---------------------------------------------------------------------- */

FixAveHisto::FixAveHisto(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), fp(nullptr), stats_list(nullptr), bin(nullptr),
    bin_total(nullptr), bin_all(nullptr), bin_list(nullptr), coord(nullptr), vector(nullptr)
{
  auto mycmd = fmt::format("fix {}", style);
  if (narg < 10) utils::missing_cmd_args(FLERR, mycmd, error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[5], false, lmp);

  global_freq = nfreq;
  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;
  dynamic_group_allow = 1;
  time_depend = 1;

  lo = utils::numeric(FLERR, arg[6], false, lmp);
  hi = utils::numeric(FLERR, arg[7], false, lmp);
  nbins = utils::inumeric(FLERR, arg[8], false, lmp);

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0 ||
        strcmp(arg[iarg],"y") == 0 ||
        strcmp(arg[iarg],"z") == 0 ||
        strcmp(arg[iarg],"vx") == 0 ||
        strcmp(arg[iarg],"vy") == 0 ||
        strcmp(arg[iarg],"vz") == 0 ||
        strcmp(arg[iarg],"fx") == 0 ||
        strcmp(arg[iarg],"fy") == 0 ||
        strcmp(arg[iarg],"fz") == 0 ||
        strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0 ||
        strncmp(arg[iarg],"v_",2) == 0) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in {} command", mycmd);

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = utils::expand_args(FLERR, nvalues, &arg[9], mode, earg, lmp);

  if (earg != &arg[9]) expand = 1;
  arg = earg;

  // parse values

  values.clear();

  for (int i = 0; i < nvalues; i++) {
    value_t val;
    val.id = "";
    val.val.c = nullptr;

    if (strcmp(arg[i],"x") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 0;
    } else if (strcmp(arg[i],"y") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 1;
    } else if (strcmp(arg[i],"z") == 0) {
      val.which = ArgInfo::X;
      val.argindex = 2;

    } else if (strcmp(arg[i],"vx") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 0;
    } else if (strcmp(arg[i],"vy") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 1;
    } else if (strcmp(arg[i],"vz") == 0) {
      val.which = ArgInfo::V;
      val.argindex = 2;

    } else if (strcmp(arg[i],"fx") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 0;
    } else if (strcmp(arg[i],"fy") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 1;
    } else if (strcmp(arg[i],"fz") == 0) {
      val.which = ArgInfo::F;
      val.argindex = 2;

    } else {
      ArgInfo argi(arg[i]);

      if (argi.get_type() == ArgInfo::NONE) break;
      if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR,"Invalid {} argument: {}", mycmd, arg[i]);

      val.which = argi.get_type();
      val.argindex = argi.get_index1();
      val.id = argi.get_name();
    }
    values.push_back(val);
  }
  if (nvalues != (int)values.size())
    error->all(FLERR, "Could not parse value data consistently for {}", mycmd);

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  // check input args for kind consistency
  // inputs must all be all either global, per-atom, or local

  if (nevery <= 0)
    error->all(FLERR,"Illegal {} nevery value: {}", mycmd, nevery);
  if (nrepeat <= 0)
    error->all(FLERR,"Illegal {} nrepeat value: {}", mycmd, nrepeat);
  if (nfreq <= 0)
    error->all(FLERR,"Illegal {} nfreq value: {}", mycmd, nfreq);
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Inconsistent {} nevery/nrepeat/nfreq values", mycmd);
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"{} overwrite keyword requires ave running setting", mycmd);

  int kindglobal,kindperatom,kindlocal;
  for (auto &val : values) {
    kindglobal = kindperatom = kindlocal = 0;

    if ((val.which == ArgInfo::X) || (val.which == ArgInfo::V) || (val.which == ArgInfo::F)) {
      kindperatom = 1;

    } else if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR,"Compute ID {} for {} does not exist", val.id, mycmd);
      // computes can produce multiple kinds of output
      if (val.val.c->scalar_flag || val.val.c->vector_flag || val.val.c->array_flag)
        kindglobal = 1;
      if (val.val.c->peratom_flag) kindperatom = 1;
      if (val.val.c->local_flag) kindlocal = 1;

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for {} does not exist", val.id, mycmd);
      // fixes can produce multiple kinds of output
      if (val.val.f->scalar_flag || val.val.f->vector_flag || val.val.f->array_flag)
        kindglobal = 1;
      if (val.val.f->peratom_flag) kindperatom = 1;
      if (val.val.f->local_flag) kindlocal = 1;

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for {} does not exist", val.id, mycmd);
      // variables only produce one kind of output
      if (input->variable->equalstyle(val.val.v) || input->variable->vectorstyle(val.val.v))
          kindglobal = 1;
      else if (input->variable->atomstyle(val.val.v)) kindperatom = 1;
      else error->all(FLERR,"{} variable {} is incompatible style", mycmd, val.id);
    }

    if (kind == DEFAULT) {
      if (kindglobal + kindperatom + kindlocal > 1)
        error->all(FLERR,"{} input kind is ambiguous", mycmd);
      if (kindglobal) kind = GLOBAL;
      if (kindperatom) kind = PERATOM;
      if (kindlocal) kind = LOCAL;
    } else if (kind == GLOBAL) {
      if (!kindglobal)
        error->all(FLERR,"{} input kind is not global", mycmd);
    } else if (kind == PERATOM) {
      if (!kindperatom)
        error->all(FLERR,"{} input kind is not peratom", mycmd);
    } else if (kind == LOCAL) {
      if (!kindlocal)
        error->all(FLERR,"{} input kind is not local", mycmd);
    }
  }

  // more error checks
  // for fix inputs, check that fix frequency is acceptable

  if (kind == PERATOM && mode == SCALAR)
    error->all(FLERR, "{} cannot process per-atom values in scalar mode", mycmd);
  if (kind == LOCAL && mode == SCALAR)
    error->all(FLERR,"{} cannot process local values in scalar mode", mycmd);

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE && kind == GLOBAL && mode == SCALAR) {
      if (val.argindex == 0 && val.val.c->scalar_flag == 0)
        error->all(FLERR, "{} compute {} does not calculate a global scalar", mycmd, val.id);
      if (val.argindex && val.val.c->vector_flag == 0)
        error->all(FLERR, "{} compute {} does not calculate a global vector", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.c->size_vector)
        error->all(FLERR, "{} compute {} vector is accessed out-of-range", mycmd, val.id);

    } else if (val.which == ArgInfo::COMPUTE && kind == GLOBAL && mode == VECTOR) {
      if (val.argindex == 0 && val.val.c->vector_flag == 0)
        error->all(FLERR, "{} compute {} does not calculate a global vector", mycmd, val.id);
      if (val.argindex && val.val.c->array_flag == 0)
        error->all(FLERR, "{} compute {} does not calculate a global array", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.c->size_array_cols)
        error->all(FLERR, "{} compute {} array is accessed out-of-range", mycmd, val.id);

    } else if (val.which == ArgInfo::COMPUTE && kind == PERATOM) {
      if (val.val.c->peratom_flag == 0)
        error->all(FLERR, "{} compute {} does not calculate per-atom values", mycmd, val.id);
      if (val.argindex == 0 && val.val.c->size_peratom_cols != 0)
        error->all(FLERR, "{} compute {} does not calculate a per-atom vector", mycmd, val.id);
      if (val.argindex && val.val.c->size_peratom_cols == 0)
        error->all(FLERR, "{} compute {} does not calculate a per-atom array", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.c->size_peratom_cols)
        error->all(FLERR, "{} compute {} array is accessed out-of-range", mycmd, val.id);

    } else if (val.which == ArgInfo::COMPUTE && kind == LOCAL) {
      if (val.val.c->local_flag == 0)
        error->all(FLERR, "{} compute {} does not calculate local values", mycmd, val.id);
      if (val.argindex == 0 && val.val.c->size_local_cols != 0)
        error->all(FLERR, "{} compute {} does not calculate a local vector", mycmd, val.id);
      if (val.argindex && val.val.c->size_local_cols == 0)
        error->all(FLERR, "{} compute {} does not calculate a local array", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.c->size_local_cols)
        error->all(FLERR, "{} compute {} array is accessed out-of-range", mycmd, val.id);

    } else if (val.which == ArgInfo::FIX && kind == GLOBAL && mode == SCALAR) {
      if (val.argindex == 0 && val.val.f->scalar_flag == 0)
        error->all(FLERR, "{} fix {} does not calculate a global scalar", mycmd, val.id);
      if (val.argindex && val.val.f->vector_flag == 0)
        error->all(FLERR, "{} fix {} does not calculate a global vector", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.f->size_vector)
        error->all(FLERR, "{} fix {} vector is accessed out-of-range", mycmd, val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for {} not computed at compatible time", val.id, mycmd);

    } else if (val.which == ArgInfo::FIX && kind == GLOBAL && mode == VECTOR) {
      if (val.argindex == 0 && val.val.f->vector_flag == 0)
        error->all(FLERR, "{} fix {} does not calculate a global vector", mycmd, val.id);
      if (val.argindex && val.val.f->array_flag == 0)
        error->all(FLERR, "{} fix {} does not calculate a global array", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.f->size_array_cols)
        error->all(FLERR, "{} fix {} array is accessed out-of-range", mycmd, val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for {} not computed at compatible time", val.id, mycmd);

    } else if (val.which == ArgInfo::FIX && kind == PERATOM) {
      if (val.val.f->peratom_flag == 0)
        error->all(FLERR, "{} fix {} does not calculate per-atom values", mycmd, val.id);
      if (val.argindex == 0 && val.val.f->size_peratom_cols != 0)
        error->all(FLERR," {} fix {} does not calculate a per-atom vector", mycmd, val.id);
      if (val.argindex && val.val.f->size_peratom_cols == 0)
        error->all(FLERR, "{} fix {} does not ""calculate a per-atom array", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.f->size_peratom_cols)
        error->all(FLERR, "{} fix {} array is accessed out-of-range", mycmd, val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for {} not computed at compatible time", val.id, mycmd);

    } else if (val.which == ArgInfo::FIX && kind == LOCAL) {
      if (val.val.f->local_flag == 0)
        error->all(FLERR, "{} fix {} does not calculate local values", mycmd, val.id);
      if (val.argindex == 0 && val.val.f->size_local_cols != 0)
        error->all(FLERR, "{} fix {} does not calculate a local vector", mycmd, val.id);
      if (val.argindex && val.val.f->size_local_cols == 0)
        error->all(FLERR, "{} fix does not calculate a local array", mycmd, val.id);
      if (val.argindex && val.argindex > val.val.f->size_local_cols)
        error->all(FLERR, "{} fix {} array is accessed out-of-range", mycmd, val.id);
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, "Fix {} for {} not computed at compatible time", val.id, mycmd);

    } else if (val.which == ArgInfo::VARIABLE && kind == GLOBAL && mode == SCALAR) {
      if (val.argindex == 0 && input->variable->equalstyle(val.val.v) == 0)
        error->all(FLERR,"{} variable {} is not equal-style variable", mycmd, val.id);
      if (val.argindex && input->variable->vectorstyle(val.val.v) == 0)
        error->all(FLERR,"{} variable {} is not vector-style variable" , mycmd, val.id);

    } else if (val.which == ArgInfo::VARIABLE && kind == GLOBAL && mode == VECTOR) {
      if (val.argindex == 0 && input->variable->vectorstyle(val.val.v) == 0)
        error->all(FLERR,"{} variable {} is not vector-style variable", mycmd, val.id);
      if (val.argindex) error->all(FLERR,"{} variable {} cannot be indexed", mycmd, val.id);

    } else if (val.which == ArgInfo::VARIABLE && kind == PERATOM) {
      if (val.argindex == 0 && input->variable->atomstyle(val.val.v) == 0)
        error->all(FLERR,"{} variable {} is not atom-style variable", mycmd, val.id);
      if (val.argindex) error->all(FLERR,"{} variable {} cannot be indexed", mycmd, val.id);
    }
  }

  // print file comment lines

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Histogrammed data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# TimeStep Number-of-bins "
                 "Total-counts Missing-counts Min-value Max-value\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else fprintf(fp,"# Bin Coord Count Count/Total\n");

    if (ferror(fp)) error->one(FLERR,"Error writing file header: {}", utils::getsyserror());
    filepos = platform::ftell(fp);
  }

  delete[] title1;
  delete[] title2;
  delete[] title3;

  // allocate and initialize memory for averaging

  if (beyond == EXTRA) nbins += 2;
  size_array_rows = nbins;

  bin = new double[nbins];
  bin_total = new double[nbins];
  bin_all = new double[nbins];
  coord = new double[nbins];

  stats_list = nullptr;
  bin_list = nullptr;
  vector = nullptr;
  maxatom = 0;

  if (ave == WINDOW) {
    memory->create(stats_list,nwindow,4,"ave/histo:stats_list");
    memory->create(bin_list,nwindow,nbins,"ave/histo:bin_list");
  }

  // initializations
  // set coord to bin centers

  if (beyond == EXTRA) {
    binsize = (hi-lo)/(nbins-2);
    bininv = 1.0/binsize;
  } else {
    binsize = (hi-lo)/nbins;
    bininv = 1.0/binsize;
  }

  if (beyond == EXTRA) {
    coord[0] = lo;
    coord[nbins-1] = hi;
    for (int i = 1; i < nbins-1; i++)
      coord[i] = lo + (i-1+0.5)*binsize;
  } else {
    for (int i = 0; i < nbins; i++)
      coord[i] = lo + (i+0.5)*binsize;
  }

  irepeat = 0;
  iwindow = window_limit = 0;

  stats_total[0] = stats_total[1] = stats_total[2] = stats_total[3] = 0.0;
  for (int i = 0; i < nbins; i++) bin_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveHisto::~FixAveHisto()
{
  if (fp && comm->me == 0) fclose(fp);

  delete[] bin;
  delete[] bin_total;
  delete[] bin_all;
  delete[] coord;
  memory->destroy(stats_list);
  memory->destroy(bin_list);
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

int FixAveHisto::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveHisto::init()
{
  auto mycmd = fmt::format("fix {}", style);

  // update indices/pointers for all computes,fixes,variables

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c) error->all(FLERR,"Compute ID {} for {} does not exist", val.id, mycmd);
    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f) error->all(FLERR,"Fix ID {} for {} does not exist", val.id, mycmd);
    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0) error->all(FLERR,"Variable name {} for {} does not exist", val.id, mycmd);
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
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveHisto::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveHisto::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  if (irepeat == 0) {
    stats[0] = stats[1] = 0.0;
    stats[2] = BIG;
    stats[3] = -BIG;
    for (int i = 0; i < nbins; i++) bin[i] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (auto &val : values) {
    int j = val.argindex;

    // atom attributes

    if (val.which == ArgInfo::X)
      bin_atoms(&atom->x[0][j],3);
    else if (val.which == ArgInfo::V)
      bin_atoms(&atom->v[0][j],3);
    else if (val.which == ArgInfo::F)
      bin_atoms(&atom->f[0][j],3);

    // invoke compute if not previously invoked

    if (val.which == ArgInfo::COMPUTE) {

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) {
          if (!(val.val.c->invoked_flag & Compute::INVOKED_SCALAR)) {
            val.val.c->compute_scalar();
            val.val.c->invoked_flag |= Compute::INVOKED_SCALAR;
          }
          bin_one(val.val.c->scalar);
        } else {
          if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
            val.val.c->compute_vector();
            val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
          }
          bin_one(val.val.c->vector[j-1]);
        }
      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
            val.val.c->compute_vector();
            val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
          }
          bin_vector(val.val.c->size_vector,val.val.c->vector,1);
        } else {
          if (!(val.val.c->invoked_flag & Compute::INVOKED_ARRAY)) {
            val.val.c->compute_array();
            val.val.c->invoked_flag |= Compute::INVOKED_ARRAY;
          }
          if (val.val.c->array)
            bin_vector(val.val.c->size_array_rows,&val.val.c->array[0][j-1],
                       val.val.c->size_array_cols);
        }

      } else if (kind == PERATOM) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
          val.val.c->compute_peratom();
          val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
        }
        if (j == 0)
          bin_atoms(val.val.c->vector_atom,1);
        else if (val.val.c->array_atom)
          bin_atoms(&val.val.c->array_atom[0][j-1],val.val.c->size_peratom_cols);

      } else if (kind == LOCAL) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_LOCAL)) {
          val.val.c->compute_local();
          val.val.c->invoked_flag |= Compute::INVOKED_LOCAL;
        }
        if (j == 0)
          bin_vector(val.val.c->size_local_rows,val.val.c->vector_local,1);
        else if (val.val.c->array_local)
          bin_vector(val.val.c->size_local_rows,&val.val.c->array_local[0][j-1],
                     val.val.c->size_local_cols);
      }

      // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) bin_one(val.val.f->compute_scalar());
        else bin_one(val.val.f->compute_vector(j-1));

      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          int n = val.val.f->size_vector;
          for (int i = 0; i < n; i++) bin_one(val.val.f->compute_vector(i));
        } else {
          int n = val.val.f->size_vector;
          for (int i = 0; i < n; i++) bin_one(val.val.f->compute_array(i,j-1));
        }

      } else if (kind == PERATOM) {
        if (j == 0) bin_atoms(val.val.f->vector_atom,1);
        else if (val.val.f->array_atom)
          bin_atoms(&val.val.f->array_atom[0][j-1],val.val.f->size_peratom_cols);

      } else if (kind == LOCAL) {
        if (j == 0) bin_vector(val.val.f->size_local_rows,val.val.f->vector_local,1);
        else if (val.val.f->array_local)
          bin_vector(val.val.f->size_local_rows,&val.val.f->array_local[0][j-1],
                     val.val.f->size_local_cols);
      }

    // evaluate equal-style or vector-style or atom-style variable
    // if index exceeds vector length, use a zero value
    //   this can be useful if vector length is not known a priori

    } else if (val.which == ArgInfo::VARIABLE) {
      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) bin_one(input->variable->compute_equal(val.val.v));
        else {
          double *varvec;
          int nvec = input->variable->compute_vector(val.val.v,&varvec);
          if (j > nvec) bin_one(0.0);
          else bin_one(varvec[j-1]);
        }

      } else if (kind == GLOBAL && mode == VECTOR) {
        double *varvec;
        int nvec = input->variable->compute_vector(val.val.v,&varvec);
        bin_vector(nvec,varvec,1);

      } else if (val.which == ArgInfo::VARIABLE && kind == PERATOM) {
        if (atom->nmax > maxatom) {
          memory->destroy(vector);
          maxatom = atom->nmax;
          memory->create(vector,maxatom,"ave/histo:vector");
        }
        input->variable->compute_atom(val.val.v,igroup,vector,1,0);
        bin_atoms(vector,1);
      }
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
  nvalid = ntimestep + nfreq - static_cast<bigint>(nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // merge histogram stats across procs if necessary

  if (kind == PERATOM || kind == LOCAL) {
    MPI_Allreduce(stats,stats_all,2,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&stats[2],&stats_all[2],1,MPI_DOUBLE,MPI_MIN,world);
    MPI_Allreduce(&stats[3],&stats_all[3],1,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(bin,bin_all,nbins,MPI_DOUBLE,MPI_SUM,world);

    stats[0] = stats_all[0];
    stats[1] = stats_all[1];
    stats[2] = stats_all[2];
    stats[3] = stats_all[3];
    for (int i = 0; i < nbins; i++) bin[i] = bin_all[i];
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    stats_total[0] = stats[0];
    stats_total[1] = stats[1];
    stats_total[2] = stats[2];
    stats_total[3] = stats[3];
    for (int i = 0; i < nbins; i++) bin_total[i] = bin[i];

  } else if (ave == RUNNING) {
    stats_total[0] += stats[0];
    stats_total[1] += stats[1];
    stats_total[2] = MIN(stats_total[2],stats[2]);
    stats_total[3] = MAX(stats_total[3],stats[3]);
    for (int i = 0; i < nbins; i++) bin_total[i] += bin[i];

  } else if (ave == WINDOW) {
    int m;
    stats_total[0] += stats[0];
    if (window_limit) stats_total[0] -= stats_list[iwindow][0];
    stats_list[iwindow][0] = stats[0];
    stats_total[1] += stats[1];
    if (window_limit) stats_total[1] -= stats_list[iwindow][1];
    stats_list[iwindow][1] = stats[1];

    if (window_limit) m = nwindow;
    else m = iwindow+1;

    stats_list[iwindow][2] = stats[2];
    stats_total[2] = stats_list[0][2];
    for (int i = 1; i < m; i++)
      stats_total[2] = MIN(stats_total[2],stats_list[i][2]);
    stats_list[iwindow][3] = stats[3];
    stats_total[3] = stats_list[0][3];
    for (int i = 1; i < m; i++)
      stats_total[3] = MAX(stats_total[3],stats_list[i][3]);

    for (int i = 0; i < nbins; i++) {
      bin_total[i] += bin[i];
      if (window_limit) bin_total[i] -= bin_list[iwindow][i];
      bin_list[iwindow][i] = bin[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
  }

  // output result to file

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (overwrite) platform::fseek(fp,filepos);
    fmt::print(fp,"{} {} {} {} {} {}\n",ntimestep,nbins,
            stats_total[0],stats_total[1],stats_total[2],stats_total[3]);
    if (stats_total[0] != 0.0)
      for (int i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",
                i+1,coord[i],bin_total[i],bin_total[i]/stats_total[0]);
    else
      for (int i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",i+1,coord[i],0.0,0.0);

    if (ferror(fp))
      error->one(FLERR,"Error writing out histogram data");

    fflush(fp);
    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR,"Error while tuncating output: {}",utils::getsyserror());
    }
  }
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveHisto::compute_vector(int i)
{
  return stats_total[i];
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveHisto::compute_array(int i, int j)
{
  if (j == 0) return coord[i];
  else if (j == 1) return bin_total[i];
  else if (stats_total[0] != 0.0) return bin_total[i]/stats_total[0];
  return 0.0;
}

/* ----------------------------------------------------------------------
   bin a single value
------------------------------------------------------------------------- */

void FixAveHisto::bin_one(double value)
{
  stats[2] = MIN(stats[2],value);
  stats[3] = MAX(stats[3],value);

  if (value < lo) {
    if (beyond == IGNORE) {
      stats[1] += 1.0;
      return;
    } else bin[0] += 1.0;
  } else if (value > hi) {
    if (beyond == IGNORE) {
      stats[1] += 1.0;
      return;
    } else bin[nbins-1] += 1.0;
  } else {
    int ibin = static_cast<int> ((value-lo)*bininv);
    ibin = MIN(ibin,nbins-1);
    if (beyond == EXTRA) ibin++;
    bin[ibin] += 1.0;
  }

  stats[0] += 1.0;
}

/* ----------------------------------------------------------------------
   bin a vector of values with stride
------------------------------------------------------------------------- */

void FixAveHisto::bin_vector(int n, double *values, int stride)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    bin_one(values[m]);
    m += stride;
  }
}

/* ----------------------------------------------------------------------
   bin a per-atom vector of values with stride
   only bin if atom is in group
------------------------------------------------------------------------- */

void FixAveHisto::bin_atoms(double *values, int stride)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) bin_one(values[m]);
    m += stride;
  }
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveHisto::options(int iarg, int narg, char **arg)
{
  // option defaults

  fp = nullptr;
  kind = DEFAULT;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  beyond = IGNORE;
  overwrite = 0;
  title1 = nullptr;
  title2 = nullptr;
  title3 = nullptr;

  // optional args
  auto mycmd = fmt::format("fix {}", style);

  while (iarg < narg) {
    if ((strcmp(arg[iarg],"file") == 0) || (strcmp(arg[iarg],"append") == 0)) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ave/histo ")+arg[iarg], error);
      if (comm->me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == nullptr)
          error->one(FLERR, "Cannot open fix ave/histo file {}: {}",
                     arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"kind") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " kind", error);
      if (strcmp(arg[iarg+1],"global") == 0) kind = GLOBAL;
      else if (strcmp(arg[iarg+1],"peratom") == 0) kind = PERATOM;
      else if (strcmp(arg[iarg+1],"local") == 0) kind = LOCAL;
      else error->all(FLERR,"Unknown fix ave/histo kind option: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " ave", error);
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Unknown fix ave/histo ave option: {}", arg[iarg+1]);
      if (ave == WINDOW) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, mycmd + " ave window", error);
        nwindow = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/histo ave window size: {}", nwindow);
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " start", error);
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " mode", error);
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Unknown fix ave/histo mode option: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"beyond") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " beyond", error);
      if (strcmp(arg[iarg+1],"ignore") == 0) beyond = IGNORE;
      else if (strcmp(arg[iarg+1],"end") == 0) beyond = END;
      else if (strcmp(arg[iarg+1],"extra") == 0) beyond = EXTRA;
      else error->all(FLERR,"Unknown fix ave/histo beyond option: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " title1", error);
      delete[] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " title2", error);
      delete[] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, mycmd + " title3", error);
      delete[] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Unknown {} option: {}", mycmd, arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveHisto::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= static_cast<bigint>(nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

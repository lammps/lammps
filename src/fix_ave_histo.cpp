// clang-format off
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

#include "fix_ave_histo.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ONE,RUNNING};
enum{SCALAR,VECTOR,WINDOW};
enum{DEFAULT,GLOBAL,PERATOM,LOCAL};
enum{IGNORE,END,EXTRA};


#define BIG 1.0e20
/* ---------------------------------------------------------------------- */

FixAveHisto::FixAveHisto(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0), which(nullptr), argindex(nullptr), value2index(nullptr),
  ids(nullptr), fp(nullptr), stats_list(nullptr),
  bin(nullptr), bin_total(nullptr), bin_all(nullptr), bin_list(nullptr),
  coord(nullptr), vector(nullptr)
{
  if (narg < 10) error->all(FLERR,"Illegal fix ave/histo command");

  MPI_Comm_rank(world,&me);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  nrepeat = utils::inumeric(FLERR,arg[4],false,lmp);
  nfreq = utils::inumeric(FLERR,arg[5],false,lmp);

  global_freq = nfreq;
  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;
  dynamic_group_allow = 1;

  lo = utils::numeric(FLERR,arg[6],false,lmp);
  hi = utils::numeric(FLERR,arg[7],false,lmp);
  nbins = utils::inumeric(FLERR,arg[8],false,lmp);

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

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/histo command");

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = utils::expand_args(FLERR,nvalues,&arg[9],mode,earg,lmp);

  if (earg != &arg[9]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nvalues];
  argindex = new int[nvalues];
  value2index = new int[nvalues];
  ids = new char*[nvalues];

  for (int i = 0; i < nvalues; i++) {
    if (strcmp(arg[i],"x") == 0) {
      which[i] = ArgInfo::X;
      argindex[i] = 0;
      ids[i] = nullptr;
    } else if (strcmp(arg[i],"y") == 0) {
      which[i] = ArgInfo::X;
      argindex[i] = 1;
      ids[i] = nullptr;
    } else if (strcmp(arg[i],"z") == 0) {
      which[i] = ArgInfo::X;
      argindex[i] = 2;
      ids[i] = nullptr;

    } else if (strcmp(arg[i],"vx") == 0) {
      which[i] = ArgInfo::V;
      argindex[i] = 0;
      ids[i] = nullptr;
    } else if (strcmp(arg[i],"vy") == 0) {
      which[i] = ArgInfo::V;
      argindex[i] = 1;
      ids[i] = nullptr;
    } else if (strcmp(arg[i],"vz") == 0) {
      which[i] = ArgInfo::V;
      argindex[i] = 2;
      ids[i] = nullptr;

    } else if (strcmp(arg[i],"fx") == 0) {
      which[i] = ArgInfo::F;
      argindex[i] = 0;
      ids[i] = nullptr;
    } else if (strcmp(arg[i],"fy") == 0) {
      which[i] = ArgInfo::F;
      argindex[i] = 1;
      ids[i] = nullptr;
    } else if (strcmp(arg[i],"fz") == 0) {
      which[i] = ArgInfo::F;
      argindex[i] = 2;
      ids[i] = nullptr;

    } else {
      ArgInfo argi(arg[i]);

      if (argi.get_type() == ArgInfo::NONE) break;
      if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR,"Invalid fix ave/histo command");

      which[i] = argi.get_type();
      argindex[i] = argi.get_index1();
      ids[i] = argi.copy_name();
    }
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // check input args for kind consistency
  // all inputs must all be global, per-atom, or local

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/histo command");
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Illegal fix ave/histo command");
  if (lo >= hi) error->all(FLERR,"Illegal fix ave/histo command");
  if (nbins <= 0) error->all(FLERR,"Illegal fix ave/histo command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/histo command");

  int kindglobal,kindperatom,kindlocal;

  for (int i = 0; i < nvalues; i++) {
    kindglobal = kindperatom = kindlocal = 0;

    if ((which[i] == ArgInfo::X) || (which[i] == ArgInfo::V)
        || (which[i] == ArgInfo::F)) {
      kindperatom = 1;

    } else if (which[i] == ArgInfo::COMPUTE) {
      int c_id = modify->find_compute(ids[i]);
      if (c_id < 0) error->all(FLERR,"Fix ave/histo input is invalid compute");
      Compute *compute = modify->compute[c_id];
      // computes can produce multiple kinds of output
      if (compute->scalar_flag || compute->vector_flag || compute->array_flag)
        kindglobal = 1;
      if (compute->peratom_flag) kindperatom = 1;
      if (compute->local_flag) kindlocal = 1;

    } else if (which[i] == ArgInfo::FIX) {
      int f_id = modify->find_fix(ids[i]);
      if (f_id < 0) error->all(FLERR,"Fix ave/histo input is invalid fix");
      Fix *fix = modify->fix[f_id];
      // fixes can produce multiple kinds of output
      if (fix->scalar_flag || fix->vector_flag || fix->array_flag)
        kindglobal = 1;
      if (fix->peratom_flag) kindperatom = 1;
      if (fix->local_flag) kindlocal = 1;

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Fix ave/histo input is invalid variable");
      // variables only produce one kind of output
      if (input->variable->equalstyle(ivariable)) kindglobal = 1;
      else if (input->variable->atomstyle(ivariable)) kindperatom = 1;
      else error->all(FLERR,"Fix ave/histo input is invalid kind of variable");
    }

    if (kind == DEFAULT) {
      if (kindglobal + kindperatom + kindlocal > 1)
        error->all(FLERR,"Fix ave/histo input kind is ambiguous");
      if (kindglobal) kind = GLOBAL;
      if (kindperatom) kind = PERATOM;
      if (kindlocal) kind = LOCAL;
    } else if (kind == GLOBAL) {
      if (!kindglobal)
        error->all(FLERR,"Fix ave/histo input kind is invalid");
    } else if (kind == PERATOM) {
      if (!kindperatom)
        error->all(FLERR,"Fix ave/histo input kind is invalid");
    } else if (kind == LOCAL) {
      if (!kindlocal)
        error->all(FLERR,"Fix ave/histo input kind is invalid");
    }
  }

  // more error checks
  // for fix inputs, check that fix frequency is acceptable

  if (kind == PERATOM && mode == SCALAR)
    error->all(FLERR,
               "Fix ave/histo cannot input per-atom values in scalar mode");
  if (kind == LOCAL && mode == SCALAR)
    error->all(FLERR,"Fix ave/histo cannot input local values in scalar mode");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE && kind == GLOBAL && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
        error->all(FLERR,
                   "Fix ave/histo compute vector is accessed out-of-range");

    } else if (which[i] == ArgInfo::COMPUTE && kind == GLOBAL && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_array_cols)
        error->all(FLERR,
                   "Fix ave/histo compute array is accessed out-of-range");

    } else if (which[i] == ArgInfo::COMPUTE && kind == PERATOM) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,
                   "Fix ave/histo compute array is accessed out-of-range");

    } else if (which[i] == ArgInfo::COMPUTE && kind == LOCAL) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      if (modify->compute[icompute]->local_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate local values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_local_cols != 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a local vector");
      if (argindex[i] && modify->compute[icompute]->size_local_cols == 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a local array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_local_cols)
        error->all(FLERR,
                   "Fix ave/histo compute array is accessed out-of-range");

    } else if (which[i] == ArgInfo::FIX && kind == GLOBAL && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate a global scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate a global vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
        error->all(FLERR,"Fix ave/histo fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == ArgInfo::FIX && kind == GLOBAL && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate a global vector");
      if (argindex[i] && modify->fix[ifix]->array_flag == 0)
        error->all(FLERR,"Fix ave/histo fix does not calculate a global array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_array_cols)
        error->all(FLERR,"Fix ave/histo fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == ArgInfo::FIX && kind == PERATOM) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,"Fix ave/histo fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == ArgInfo::FIX && kind == LOCAL) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      if (modify->fix[ifix]->local_flag == 0)
        error->all(FLERR,"Fix ave/histo fix does not calculate local values");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_local_cols != 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a local vector");
      if (argindex[i] && modify->fix[ifix]->size_local_cols == 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a local array");
      if (argindex[i] &&
          argindex[i] > modify->fix[ifix]->size_local_cols)
        error->all(FLERR,"Fix ave/histo fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == ArgInfo::VARIABLE && kind == GLOBAL && mode == SCALAR) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
      if (argindex[i] == 0 && input->variable->equalstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/histo variable is not equal-style variable");
      if (argindex[i] && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/histo variable is not vector-style variable");

    } else if (which[i] == ArgInfo::VARIABLE && kind == GLOBAL && mode == VECTOR) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
      if (argindex[i] == 0 && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/histo variable is not vector-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/histo variable cannot be indexed");

    } else if (which[i] == ArgInfo::VARIABLE && kind == PERATOM) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
      if (argindex[i] == 0 && input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/histo variable is not atom-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/histo variable cannot be indexed");
    }
  }

  // print file comment lines

  if (fp && me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Histogrammed data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# TimeStep Number-of-bins "
                 "Total-counts Missing-counts Min-value Max-value\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else fprintf(fp,"# Bin Coord Count Count/Total\n");

    if (ferror(fp))
      error->one(FLERR,"Error writing file header");

    filepos = ftell(fp);
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;

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
  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  if (fp && me == 0) fclose(fp);

  delete [] bin;
  delete [] bin_total;
  delete [] bin_all;
  delete [] coord;
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
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      value2index[i] = icompute;

    } else if (which[i] == ArgInfo::FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      value2index[i] = ifix;

    } else if (which[i] == ArgInfo::VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
      value2index[i] = ivariable;
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
  int i,j,m;

  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/histo");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  if (irepeat == 0) {
    stats[0] = stats[1] = 0.0;
    stats[2] = BIG;
    stats[3] = -BIG;
    for (i = 0; i < nbins; i++) bin[i] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (i = 0; i < nvalues; i++) {
    m = value2index[i];
    j = argindex[i];

    // atom attributes

    if (which[i] == ArgInfo::X)
      bin_atoms(&atom->x[0][j],3);
    else if (which[i] == ArgInfo::V)
      bin_atoms(&atom->v[0][j],3);
    else if (which[i] == ArgInfo::F)
      bin_atoms(&atom->f[0][j],3);

    // invoke compute if not previously invoked

    if (which[i] == ArgInfo::COMPUTE) {
      Compute *compute = modify->compute[m];

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) {
          if (!(compute->invoked_flag & Compute::INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= Compute::INVOKED_SCALAR;
          }
          bin_one(compute->scalar);
        } else {
          if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= Compute::INVOKED_VECTOR;
          }
          bin_one(compute->vector[j-1]);
        }
      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= Compute::INVOKED_VECTOR;
          }
          bin_vector(compute->size_vector,compute->vector,1);
        } else {
          if (!(compute->invoked_flag & Compute::INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= Compute::INVOKED_ARRAY;
          }
          if (compute->array)
            bin_vector(compute->size_array_rows,&compute->array[0][j-1],
                       compute->size_array_cols);
        }

      } else if (kind == PERATOM) {
        if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
          compute->compute_peratom();
          compute->invoked_flag |= Compute::INVOKED_PERATOM;
        }
        if (j == 0)
          bin_atoms(compute->vector_atom,1);
        else if (compute->array_atom)
          bin_atoms(&compute->array_atom[0][j-1],compute->size_peratom_cols);

      } else if (kind == LOCAL) {
        if (!(compute->invoked_flag & Compute::INVOKED_LOCAL)) {
          compute->compute_local();
          compute->invoked_flag |= Compute::INVOKED_LOCAL;
        }
        if (j == 0)
          bin_vector(compute->size_local_rows,compute->vector_local,1);
        else if (compute->array_local)
          bin_vector(compute->size_local_rows,&compute->array_local[0][j-1],
                     compute->size_local_cols);
      }

      // access fix fields, guaranteed to be ready

    } else if (which[i] == ArgInfo::FIX) {

      Fix *fix = modify->fix[m];

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) bin_one(fix->compute_scalar());
        else bin_one(fix->compute_vector(j-1));

      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          int n = fix->size_vector;
          for (i = 0; i < n; i++) bin_one(fix->compute_vector(i));
        } else {
          int n = fix->size_vector;
          for (i = 0; i < n; i++) bin_one(fix->compute_array(i,j-1));
        }

      } else if (kind == PERATOM) {
        if (j == 0) bin_atoms(fix->vector_atom,1);
        else if (fix->array_atom)
          bin_atoms(fix->array_atom[j-1],fix->size_peratom_cols);

      } else if (kind == LOCAL) {
        if (j == 0) bin_vector(fix->size_local_rows,fix->vector_local,1);
        else if (fix->array_local)
          bin_vector(fix->size_local_rows,&fix->array_local[0][j-1],
                     fix->size_local_cols);
      }

    // evaluate equal-style or vector-style or atom-style variable

    } else if (which[i] == ArgInfo::VARIABLE) {
      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) bin_one(input->variable->compute_equal(m));
        else {
          double *varvec;
          int nvec = input->variable->compute_vector(m,&varvec);
          if (nvec < j) bin_one(0.0);
          else bin_one(varvec[j-1]);
        }

      } else if (kind == GLOBAL && mode == VECTOR) {
        double *varvec;
        int nvec = input->variable->compute_vector(m,&varvec);
        bin_vector(nvec,varvec,1);

      } else if (which[i] == ArgInfo::VARIABLE && kind == PERATOM) {
        if (atom->nmax > maxatom) {
          memory->destroy(vector);
          maxatom = atom->nmax;
          memory->create(vector,maxatom,"ave/histo:vector");
        }
        input->variable->compute_atom(m,igroup,vector,1,0);
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
    for (i = 0; i < nbins; i++) bin[i] = bin_all[i];
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    stats_total[0] = stats[0];
    stats_total[1] = stats[1];
    stats_total[2] = stats[2];
    stats_total[3] = stats[3];
    for (i = 0; i < nbins; i++) bin_total[i] = bin[i];

  } else if (ave == RUNNING) {
    stats_total[0] += stats[0];
    stats_total[1] += stats[1];
    stats_total[2] = MIN(stats_total[2],stats[2]);
    stats_total[3] = MAX(stats_total[3],stats[3]);
    for (i = 0; i < nbins; i++) bin_total[i] += bin[i];

  } else if (ave == WINDOW) {
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
    for (i = 1; i < m; i++)
      stats_total[2] = MIN(stats_total[2],stats_list[i][2]);
    stats_list[iwindow][3] = stats[3];
    stats_total[3] = stats_list[0][3];
    for (i = 1; i < m; i++)
      stats_total[3] = MAX(stats_total[3],stats_list[i][3]);

    for (i = 0; i < nbins; i++) {
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

  if (fp && me == 0) {
    clearerr(fp);
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d %g %g %g %g\n",ntimestep,nbins,
            stats_total[0],stats_total[1],stats_total[2],stats_total[3]);
    if (stats_total[0] != 0.0)
      for (i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",
                i+1,coord[i],bin_total[i],bin_total[i]/stats_total[0]);
    else
      for (i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",i+1,coord[i],0.0,0.0);

    if (ferror(fp))
      error->one(FLERR,"Error writing out histogram data");

    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      if ((fileend > 0) && (ftruncate(fileno(fp),fileend)))
        perror("Error while tuncating output");
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

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == nullptr)
          error->one(FLERR,"Cannot open fix ave/histo file {}: {}",
                                       arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"kind") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"global") == 0) kind = GLOBAL;
      else if (strcmp(arg[iarg+1],"peratom") == 0) kind = PERATOM;
      else if (strcmp(arg[iarg+1],"local") == 0) kind = LOCAL;
      else error->all(FLERR,"Illegal fix ave/histo command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/histo command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/histo command");
        nwindow = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/histo command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Illegal fix ave/histo command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"beyond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"ignore") == 0) beyond = IGNORE;
      else if (strcmp(arg[iarg+1],"end") == 0) beyond = END;
      else if (strcmp(arg[iarg+1],"extra") == 0) beyond = EXTRA;
      else error->all(FLERR,"Illegal fix ave/histo command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      delete [] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      delete [] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      delete [] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/histo command");
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

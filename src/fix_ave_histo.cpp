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

#include "stdlib.h"
#include "string.h"
#include "fix_ave_histo.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};
enum{SCALAR,VECTOR,WINDOW};
enum{GLOBAL,PERATOM,LOCAL};
enum{IGNORE,END,EXTRA};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

#define BIG 1.0e20
/* ---------------------------------------------------------------------- */

FixAveHisto::FixAveHisto(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix ave/histo command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);

  global_freq = nfreq;
  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;

  lo = force->numeric(FLERR,arg[6]);
  hi = force->numeric(FLERR,arg[7]);
  nbins = force->inumeric(FLERR,arg[8]);

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

  options(narg,arg);

  // parse values until one isn't recognized
  // if mode = VECTOR and value is a global array:
  //   expand it as if columns listed one by one
  //   adjust nvalues accordingly via maxvalues

  which = argindex = value2index = NULL;
  ids = NULL;
  int maxvalues = nvalues;
  allocate_values(maxvalues);
  nvalues = 0;

  iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      which[nvalues] = X;
      argindex[nvalues] = 0;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;
    } else if (strcmp(arg[iarg],"y") == 0) {
      which[nvalues] = X;
      argindex[nvalues] = 1;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;
    } else if (strcmp(arg[iarg],"z") == 0) {
      which[nvalues] = X;
      argindex[nvalues] = 2;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = V;
      argindex[nvalues] = 0;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = V;
      argindex[nvalues] = 1;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = V;
      argindex[nvalues] = 2;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = F;
      argindex[nvalues] = 0;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = F;
      argindex[nvalues] = 1;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = F;
      argindex[nvalues] = 2;
      ids[nvalues] = NULL;
      nvalues++;
      iarg++;

    } else if ((strncmp(arg[iarg],"c_",2) == 0) ||
        (strncmp(arg[iarg],"f_",2) == 0) ||
        (strncmp(arg[iarg],"v_",2) == 0)) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ave/histo command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      delete [] suffix;

      if (mode == VECTOR && which[nvalues] == COMPUTE &&
          argindex[nvalues] == 0) {
        int icompute = modify->find_compute(ids[nvalues]);
        if (icompute < 0)
          error->all(FLERR,"Compute ID for fix ave/histo does not exist");
        if (modify->compute[icompute]->array_flag) {
          int ncols = modify->compute[icompute]->size_array_cols;
          maxvalues += ncols-1;
          allocate_values(maxvalues);
          argindex[nvalues] = 1;
          for (int icol = 1; icol < ncols; icol++) {
            which[nvalues+icol] = which[nvalues];
            argindex[nvalues+icol] = icol+1;
            n = strlen(ids[nvalues]) + 1;
            ids[nvalues+icol] = new char[n];
            strcpy(ids[nvalues+icol],ids[nvalues]);
          }
          nvalues += ncols-1;
        }

      } else if (mode == VECTOR && which[nvalues] == FIX &&
                 argindex[nvalues] == 0) {
        int ifix = modify->find_fix(ids[nvalues]);
        if (ifix < 0)
          error->all(FLERR,"Fix ID for fix ave/histo does not exist");
        if (modify->fix[ifix]->array_flag) {
          int ncols = modify->fix[ifix]->size_array_cols;
          maxvalues += ncols-1;
          allocate_values(maxvalues);
          argindex[nvalues] = 1;
          for (int icol = 1; icol < ncols; icol++) {
            which[nvalues+icol] = which[nvalues];
            argindex[nvalues+icol] = icol+1;
            n = strlen(ids[nvalues]) + 1;
            ids[nvalues+icol] = new char[n];
            strcpy(ids[nvalues+icol],ids[nvalues]);
          }
          nvalues += ncols-1;
        }
      }

      nvalues++;
      iarg++;
    } else break;
  }

  // setup and error check
  // kind = inputs are all global, or all per-atom, or all local
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/histo command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix ave/histo command");
  if (lo >= hi) error->all(FLERR,"Illegal fix ave/histo command");
  if (nbins <= 0) error->all(FLERR,"Illegal fix ave/histo command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/histo command");

  int kindflag;
  for (int i = 0; i < nvalues; i++) {
    if (which[i] == X || which[i] == V || which[i] == F) kindflag = PERATOM;
    else if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[modify->find_compute(ids[0])];
      if (compute->scalar_flag || compute->vector_flag || compute->array_flag)
        kindflag = GLOBAL;
      else if (compute->peratom_flag) kindflag = PERATOM;
      else if (compute->local_flag) kindflag = LOCAL;
      else error->all(FLERR,"Fix ave/histo input is invalid compute");
    } else if (which[i] == FIX) {
      Fix *fix = modify->fix[modify->find_fix(ids[0])];
      if (fix->scalar_flag || fix->vector_flag || fix->array_flag)
        kindflag = GLOBAL;
      else if (fix->peratom_flag) kindflag = PERATOM;
      else if (fix->local_flag) kindflag = LOCAL;
      else error->all(FLERR,"Fix ave/histo input is invalid fix");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (input->variable->equalstyle(ivariable)) kindflag = GLOBAL;
      else if (input->variable->atomstyle(ivariable)) kindflag = PERATOM;
      else error->all(FLERR,"Fix ave/histo input is invalid variable");
    }
    if (i == 0) kind = kindflag;
    else if (kindflag != kind)
      error->all(FLERR,
                 "Fix ave/histo inputs are not all global, peratom, or local");
  }

  if (kind == PERATOM && mode == SCALAR)
    error->all(FLERR,
               "Fix ave/histo cannot input per-atom values in scalar mode");
  if (kind == LOCAL && mode == SCALAR)
    error->all(FLERR,"Fix ave/histo cannot input local values in scalar mode");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE && kind == GLOBAL && mode == SCALAR) {
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

    } else if (which[i] == COMPUTE && kind == GLOBAL && mode == VECTOR) {
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

    } else if (which[i] == COMPUTE && kind == PERATOM) {
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

    } else if (which[i] == COMPUTE && kind == LOCAL) {
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

    } else if (which[i] == FIX && kind == GLOBAL && mode == SCALAR) {
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

    } else if (which[i] == FIX && kind == GLOBAL && mode == VECTOR) {
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

    } else if (which[i] == FIX && kind == PERATOM) {
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

    } else if (which[i] == FIX && kind == LOCAL) {
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

    } else if (which[i] == VARIABLE && kind == GLOBAL) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");

    } else if (which[i] == VARIABLE && kind == PERATOM) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
    }
  }

  // print file comment lines

  if (fp && me == 0) {
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Histogrammed data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# TimeStep Number-of-bins "
                 "Total-counts Missing-counts Min-value Max-value\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else fprintf(fp,"# Bin Coord Count Count/Total\n");
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

  stats_list = NULL;
  bin_list = NULL;
  vector = NULL;
  maxatom = 0;
  if (ave == WINDOW) {
    memory->create(stats_list,nwindow,4,"histo:stats_list");
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

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveHisto::~FixAveHisto()
{
  memory->destroy(which);
  memory->destroy(argindex);
  memory->destroy(value2index);
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  memory->sfree(ids);

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
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      value2index[i] = ifix;

    } else if (which[i] == VARIABLE) {
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

void FixAveHisto::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveHisto::end_of_step()
{
  int i,j,m;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

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

    if (which[i] == X)
      bin_atoms(&atom->x[0][j],3);
    else if (which[i] == V)
      bin_atoms(&atom->v[0][j],3);
    else if (which[i] == F)
      bin_atoms(&atom->f[0][j],3);

    // invoke compute if not previously invoked

    if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) {
          if (!(compute->invoked_flag & INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= INVOKED_SCALAR;
          }
          bin_one(compute->scalar);
        } else {
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }
          bin_one(compute->vector[j-1]);
        }
      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }
          bin_vector(compute->size_vector,compute->vector,1);
        } else {
          if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }
          if (compute->array)
            bin_vector(compute->size_array_rows,&compute->array[0][j-1],
                       compute->size_array_cols);
        }

      } else if (kind == PERATOM) {
        if (!(compute->invoked_flag & INVOKED_PERATOM)) {
          compute->compute_peratom();
          compute->invoked_flag |= INVOKED_PERATOM;
        }
        if (j == 0)
          bin_atoms(compute->vector_atom,1);
        else if (compute->array_atom)
          bin_atoms(&compute->array_atom[0][j-1],compute->size_peratom_cols);

      } else if (kind == LOCAL) {
        if (!(compute->invoked_flag & INVOKED_LOCAL)) {
          compute->compute_local();
          compute->invoked_flag |= INVOKED_LOCAL;
        }
        if (j == 0)
          bin_vector(compute->size_local_rows,compute->vector_local,1);
        else if (compute->array_local)
          bin_vector(compute->size_local_rows,&compute->array_local[0][j-1],
                     compute->size_local_cols);
      }

      // access fix fields, guaranteed to be ready

    } else if (which[i] == FIX) {

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

      // evaluate equal-style variable

    } else if (which[i] == VARIABLE && kind == GLOBAL) {
      bin_one(input->variable->compute_equal(m));

    } else if (which[i] == VARIABLE && kind == PERATOM) {
      if (atom->nlocal > maxatom) {
        memory->destroy(vector);
        maxatom = atom->nmax;
        memory->create(vector,maxatom,"ave/histo:vector");
      }
      input->variable->compute_atom(m,igroup,vector,1,0);
      bin_atoms(vector,1);
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
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;
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
    fflush(fp);
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

void FixAveHisto::options(int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  beyond = IGNORE;
  overwrite = 0;
  title1 = NULL;
  title2 = NULL;
  title3 = NULL;

  // optional args

  int iarg = 9 + nvalues;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix ave/histo file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/histo command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/histo command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/histo command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
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
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/histo command");
  }
}

/* ----------------------------------------------------------------------
   reallocate vectors for each input value, of length N
------------------------------------------------------------------------- */

void FixAveHisto::allocate_values(int n)
{
  memory->grow(which,n,"ave/time:which");
  memory->grow(argindex,n,"ave/time:argindex");
  memory->grow(value2index,n,"ave/time:value2index");
  ids = (char **) memory->srealloc(ids,n*sizeof(char *),"ave/time:ids");
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
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixAveHisto::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix ave/histo missed timestep");
}

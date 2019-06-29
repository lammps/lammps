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

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdio>    // IWYU pragma: keep
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "fix_ave_time.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING,WINDOW};
enum{SCALAR,VECTOR};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0), which(NULL), argindex(NULL), value2index(NULL),
  offcol(NULL), varlen(NULL), ids(NULL),
  fp(NULL), offlist(NULL), format(NULL), format_user(NULL),
  vector(NULL), vector_total(NULL), vector_list(NULL),
  column(NULL), array(NULL), array_total(NULL), array_list(NULL)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/time command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);

  global_freq = nfreq;

  dynamic_group_allow = 1;

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;

  int iarg = 6;
  while (iarg < narg) {
    if ((strncmp(arg[iarg],"c_",2) == 0) ||
        (strncmp(arg[iarg],"f_",2) == 0) ||
        (strncmp(arg[iarg],"v_",2) == 0)) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/time command");

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = input->expand_args(nvalues,&arg[6],mode,earg);

  if (earg != &arg[6]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nvalues];
  argindex = new int[nvalues];
  value2index = new int[nvalues];
  offcol = new int[nvalues];
  varlen = new int[nvalues];
  ids = new char*[nvalues];

  for (int i = 0; i < nvalues; i++) {
    if (arg[i][0] == 'c') which[i] = COMPUTE;
    else if (arg[i][0] == 'f') which[i] = FIX;
    else if (arg[i][0] == 'v') which[i] = VARIABLE;

    int n = strlen(arg[i]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[i][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal fix ave/time command");
      argindex[i] = atoi(ptr+1);
      *ptr = '\0';
    } else argindex[i] = 0;

    n = strlen(suffix) + 1;
    ids[i] = new char[n];
    strcpy(ids[i],suffix);
    delete [] suffix;
  }

  // set off columns now that nvalues is finalized

  for (int i = 0; i < nvalues; i++) offcol[i] = 0;
  for (int i = 0; i < noff; i++) {
    if (offlist[i] < 1 || offlist[i] > nvalues)
      error->all(FLERR,"Invalid fix ave/time off column");
    offcol[offlist[i]-1] = 1;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable
  // set variable_length if any compute is variable length

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/time command");
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Illegal fix ave/time command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/time command");

  for (int i = 0; i < nvalues; i++) {
    varlen[i] = 0;

    if (which[i] == COMPUTE && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector &&
          modify->compute[icompute]->size_vector_variable == 0)
        error->all(FLERR,
                   "Fix ave/time compute vector is accessed out-of-range");
      if (argindex[i] && modify->compute[icompute]->size_vector_variable)
        varlen[i] = 1;

    } else if (which[i] == COMPUTE && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate a vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate an array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_array_cols)
        error->all(FLERR,"Fix ave/time compute array is accessed out-of-range");
      if (argindex[i] == 0 && modify->compute[icompute]->size_vector_variable)
        varlen[i] = 1;
      if (argindex[i] && modify->compute[icompute]->size_array_rows_variable)
        varlen[i] = 1;

    } else if (which[i] == FIX && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
        error->all(FLERR,"Fix ave/time fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time fix does not calculate a vector");
      if (argindex[i] && modify->fix[ifix]->size_vector_variable)
        error->all(FLERR,"Fix ave/time fix vector cannot be variable length");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
        error->all(FLERR,"Fix ave/time fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == FIX && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time fix does not calculate a vector");
      if (argindex[i] && modify->fix[ifix]->array_flag == 0)
        error->all(FLERR,"Fix ave/time fix does not calculate an array");
      if (argindex[i] && modify->fix[ifix]->size_array_rows_variable)
        error->all(FLERR,"Fix ave/time fix array cannot be variable length");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_array_cols)
        error->all(FLERR,"Fix ave/time fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == VARIABLE && mode == SCALAR) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/time does not exist");
      if (argindex[i] == 0 && input->variable->equalstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/time variable is not equal-style variable");
      if (argindex[i] && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/time variable is not vector-style variable");

    } else if (which[i] == VARIABLE && mode == VECTOR) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/time does not exist");
      if (argindex[i] == 0 && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/time variable is not vector-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/time mode vector variable cannot be indexed");
      varlen[i] = 1;
    }
  }

  // all_variable_length = 1 if all values are variable length
  // any_variable_length = 1 if any values are variable length

  all_variable_length = 1;
  any_variable_length = 0;
  for (int i = 0; i < nvalues; i++) {
    if (varlen[i] == 0) all_variable_length = 0;
    if (varlen[i]) any_variable_length = 1;
  }

  // if VECTOR mode, check that all columns are same length
  // nrows = # of rows in output array
  // if all columns are variable length, just set nrows = 1 for now

  column = NULL;
  if (mode == VECTOR) {
    if (all_variable_length == 0) nrows = column_length(0);
    else nrows = 1;
    memory->create(column,nrows,"ave/time:column");
  }

  // enable locking of row count by this fix for computes of variable length
  // only if nrepeat > 1 or ave = RUNNING/WINDOW,
  //   so that locking spans multiple timesteps

  if (any_variable_length &&
      (nrepeat > 1 || ave == RUNNING || ave == WINDOW)) {
    for (int i = 0; i < nvalues; i++)
      if (varlen[i] && which[i] == COMPUTE) {
        int icompute = modify->find_compute(ids[i]);
        modify->compute[icompute]->lock_enable();
      }
    lockforever = 0;
  }

  // print file comment lines
  // for mode = VECTOR, cannot use arg to print
  // since array args may have been expanded to multiple vectors

  if (fp && me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Time-averaged data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else if (mode == SCALAR) {
      fprintf(fp,"# TimeStep");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    } else fprintf(fp,"# TimeStep Number-of-rows\n");
    if (title3 && mode == VECTOR) fprintf(fp,"%s\n",title3);
    else if (mode == VECTOR) {
      fprintf(fp,"# Row");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    }
    if (ferror(fp))
      error->one(FLERR,"Error writing file header");

    filepos = ftell(fp);
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // allocate memory for averaging

  vector = vector_total = NULL;
  vector_list = NULL;
  array = array_total = NULL;
  array_list = NULL;

  if (mode == SCALAR) {
    vector = new double[nvalues];
    vector_total = new double[nvalues];
    if (ave == WINDOW)
      memory->create(vector_list,nwindow,nvalues,"ave/time:vector_list");
  } else allocate_arrays();

  // this fix produces either a global scalar or vector or array
  // SCALAR mode produces either a scalar or vector
  // VECTOR mode produces either a vector or array
  // intensive/extensive flags set by compute,fix,variable that produces value

  extlist = NULL;

  if (mode == SCALAR) {
    if (nvalues == 1) {
      scalar_flag = 1;
      if (which[0] == COMPUTE) {
        Compute *compute = modify->compute[modify->find_compute(ids[0])];
        if (argindex[0] == 0) extscalar = compute->extscalar;
        else if (compute->extvector >= 0) extscalar = compute->extvector;
        else extscalar = compute->extlist[argindex[0]-1];
      } else if (which[0] == FIX) {
        Fix *fix = modify->fix[modify->find_fix(ids[0])];
        if (argindex[0] == 0) extscalar = fix->extscalar;
        else if (fix->extvector >= 0) extscalar = fix->extvector;
        else extscalar = fix->extlist[argindex[0]-1];
      } else if (which[0] == VARIABLE) {
        extscalar = 0;
      }

    } else {
      vector_flag = 1;
      size_vector = nrows = nvalues;
      extvector = -1;
      extlist = new int[nvalues];
      for (int i = 0; i < nvalues; i++) {
        if (which[i] == COMPUTE) {
          Compute *compute = modify->compute[modify->find_compute(ids[i])];
          if (argindex[i] == 0) extlist[i] = compute->extscalar;
          else if (compute->extvector >= 0) extlist[i] = compute->extvector;
          else extlist[i] = compute->extlist[argindex[i]-1];
        } else if (which[i] == FIX) {
          Fix *fix = modify->fix[modify->find_fix(ids[i])];
          if (argindex[i] == 0) extlist[i] = fix->extscalar;
          else if (fix->extvector >= 0) extlist[i] = fix->extvector;
          else extlist[i] = fix->extlist[argindex[i]-1];
        } else if (which[i] == VARIABLE) {
          extlist[i] = 0;
        }
      }
    }

  } else {
    if (nvalues == 1) {
      vector_flag = 1;
      size_vector = nrows;
      if (all_variable_length) size_vector_variable = 1;
      if (which[0] == COMPUTE) {
        Compute *compute = modify->compute[modify->find_compute(ids[0])];
        if (argindex[0] == 0) {
          extvector = compute->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = compute->extlist[i];
          }
        } else extvector = compute->extarray;
      } else if (which[0] == FIX) {
        Fix *fix = modify->fix[modify->find_fix(ids[0])];
        if (argindex[0] == 0) {
          extvector = fix->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = fix->extlist[i];
          }
        } else extvector = fix->extarray;
      } else if (which[0] == VARIABLE) {
        extlist = new int[nrows];
        for (int i = 0; i < nrows; i++) extlist[i] = 0;
      }

    } else {
      array_flag = 1;
      size_array_rows = nrows;
      size_array_cols = nvalues;
      if (all_variable_length) size_array_rows_variable = 1;
      int value;
      for (int i = 0; i < nvalues; i++) {
        if (which[i] == COMPUTE) {
          Compute *compute = modify->compute[modify->find_compute(ids[i])];
          if (argindex[i] == 0) value = compute->extvector;
          else value = compute->extarray;
        } else if (which[i] == FIX) {
          Fix *fix = modify->fix[modify->find_fix(ids[i])];
          if (argindex[i] == 0) value = fix->extvector;
          else value = fix->extarray;
        } else if (which[i] == VARIABLE) {
          value = 0;
        }
        if (value == -1)
          error->all(FLERR,"Fix ave/time cannot set output array "
                     "intensive/extensive from these inputs");
        if (i == 0) extarray = value;
        else if (value != extarray)
          error->all(FLERR,"Fix ave/time cannot set output array "
                     "intensive/extensive from these inputs");
      }
    }
  }

  // initializations
  // set vector_total to zero since it accumulates
  // array_total already zeroed in allocate_arrays

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  if (mode == SCALAR)
    for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  // decrement lock counter in compute chunk/atom, it if still exists

  if (any_variable_length &&
      (nrepeat > 1 || ave == RUNNING || ave == WINDOW)) {
    for (int i = 0; i < nvalues; i++)
      if (varlen[i]) {
        int icompute = modify->find_compute(ids[i]);
        if (icompute >= 0) {
          if (ave == RUNNING || ave == WINDOW)
            modify->compute[icompute]->unlock(this);
          modify->compute[icompute]->lock_disable();
        }
      }
  }

  delete [] format_user;

  delete [] which;
  delete [] argindex;
  delete [] value2index;
  delete [] offcol;
  delete [] varlen;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  delete [] extlist;

  if (fp && me == 0) fclose(fp);

  memory->destroy(column);

  delete [] vector;
  delete [] vector_total;
  memory->destroy(array);
  memory->destroy(array_total);
  memory->destroy(array_list);
}

/* ---------------------------------------------------------------------- */

int FixAveTime::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveTime::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/time does not exist");
      value2index[i] = icompute;
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/time does not exist");
      value2index[i] = ifix;
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/time does not exist");
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

void FixAveTime::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveTime::end_of_step()
{
  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/time");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  if (mode == SCALAR) invoke_scalar(ntimestep);
  else invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_scalar(bigint ntimestep)
{
  int i,m;
  double scalar;

  // zero if first sample within single Nfreq epoch
  // if any input is variable length, initialize current length
  // check for exceeding length is done below

  if (irepeat == 0) {
    if (any_variable_length) {
      modify->clearstep_compute();
      column_length(1);
      modify->addstep_compute(ntimestep+nevery);
      modify->addstep_compute(ntimestep+nfreq);
    }
    for (i = 0; i < nvalues; i++) vector[i] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (i = 0; i < nvalues; i++) {
    m = value2index[i];

    // invoke compute if not previously invoked
    // insure no out-of-range access to variable-length compute vector

    if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[i] == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        scalar = compute->scalar;
      } else {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        if (varlen[i] && compute->size_vector < argindex[i]) scalar = 0.0;
        else scalar = compute->vector[argindex[i]-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == FIX) {
      if (argindex[i] == 0)
        scalar = modify->fix[m]->compute_scalar();
      else
        scalar = modify->fix[m]->compute_vector(argindex[i]-1);

    // evaluate equal-style or vector-style variable
    // insure no out-of-range access to vector-style variable

    } else if (which[i] == VARIABLE) {
      if (argindex[i] == 0)
        scalar = input->variable->compute_equal(m);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(m,&varvec);
        if (nvec < argindex[i]) scalar = 0.0;
        else scalar = varvec[argindex[i]-1];
      }
    }

    // add value to vector or just set directly if offcol is set

    if (offcol[i]) vector[i] = scalar;
    else vector[i] += scalar;
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

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (i = 0; i < nvalues; i++)
    if (offcol[i] == 0) vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (i = 0; i < nvalues; i++) vector_total[i] = vector[i];
    norm = 1;

  } else if (ave == RUNNING) {
    for (i = 0; i < nvalues; i++) vector_total[i] += vector[i];
    norm++;

  } else if (ave == WINDOW) {
    for (i = 0; i < nvalues; i++) {
      vector_total[i] += vector[i];
      if (window_limit) vector_total[i] -= vector_list[iwindow][i];
      vector_list[iwindow][i] = vector[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // insure any columns with offcol set are effectively set to last value

  for (i = 0; i < nvalues; i++)
    if (offcol[i]) vector_total[i] = norm*vector[i];

  // output result to file

  if (fp && me == 0) {
    clearerr(fp);
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT,ntimestep);
    for (i = 0; i < nvalues; i++) fprintf(fp,format,vector_total[i]/norm);
    fprintf(fp,"\n");
    if (ferror(fp))
      error->one(FLERR,"Error writing out time averaged data");

    fflush(fp);

    if (overwrite) {
      long fileend = ftell(fp);
      if (fileend > 0) ftruncate(fileno(fp),fileend);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_vector(bigint ntimestep)
{
  int i,j,m;

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
    if (any_variable_length) {
      modify->clearstep_compute();
      int nrows_new = column_length(1);
      modify->addstep_compute(ntimestep+nevery);
      modify->addstep_compute(ntimestep+nfreq);

      if (all_variable_length && nrows_new != nrows) {
        nrows = nrows_new;
        memory->destroy(column);
        memory->create(column,nrows,"ave/time:column");
        allocate_arrays();
      }

      bigint ntimestep = update->ntimestep;
      int lockforever_flag = 0;
      for (i = 0; i < nvalues; i++) {
        if (!varlen[i] || which[i] != COMPUTE) continue;
        if (nrepeat > 1 && ave == ONE) {
          Compute *compute = modify->compute[value2index[i]];
          compute->lock(this,ntimestep,ntimestep+(nrepeat-1)*nevery);
        } else if ((ave == RUNNING || ave == WINDOW) && !lockforever) {
          Compute *compute = modify->compute[value2index[i]];
          compute->lock(this,update->ntimestep,-1);
          lockforever_flag = 1;
        }
      }
      if (lockforever_flag) lockforever = 1;
    }

    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array[i][j] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (j = 0; j < nvalues; j++) {
    m = value2index[j];

    // invoke compute if not previously invoked

    if (which[j] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[j] == 0) {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        double *cvector = compute->vector;
        for (i = 0; i < nrows; i++)
          column[i] = cvector[i];

      } else {
        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        double **carray = compute->array;
        int icol = argindex[j]-1;
        for (i = 0; i < nrows; i++)
          column[i] = carray[i][icol];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[j] == FIX) {
      Fix *fix = modify->fix[m];
      if (argindex[j] == 0)
        for (i = 0; i < nrows; i++)
          column[i] = fix->compute_vector(i);
      else {
        int icol = argindex[j]-1;
        for (i = 0; i < nrows; i++)
          column[i] = fix->compute_array(i,icol);
      }

    // evaluate vector-style variable
    // insure nvec = nrows, else error
    // could be different on this timestep than when column_length(1) set nrows

    } else if (which[j] == VARIABLE) {
      double *varvec;
      int nvec = input->variable->compute_vector(m,&varvec);
      if (nvec != nrows)
        error->all(FLERR,"Fix ave/time vector-style variable changed length");
      for (i = 0; i < nrows; i++)
        column[i] = varvec[i];
    }

    // add columns of values to array or just set directly if offcol is set

    if (offcol[j]) {
      for (i = 0; i < nrows; i++)
        array[i][j] = column[i];
    } else {
      for (i = 0; i < nrows; i++)
        array[i][j] += column[i];
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
  nvalid = ntimestep+nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // unlock any variable length computes at end of Nfreq epoch
  // do not unlock if ave = RUNNING or WINDOW

  if (any_variable_length && nrepeat > 1 && ave == ONE) {
    for (i = 0; i < nvalues; i++) {
      if (!varlen[i]) continue;
      Compute *compute = modify->compute[value2index[i]];
      compute->unlock(this);
    }
  }

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j] == 0) array[i][j] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array_total[i][j] = array[i][j];
    norm = 1;

  } else if (ave == RUNNING) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array_total[i][j] += array[i][j];
    norm++;

  } else if (ave == WINDOW) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) {
        array_total[i][j] += array[i][j];
        if (window_limit) array_total[i][j] -= array_list[iwindow][i][j];
        array_list[iwindow][i][j] = array[i][j];
      }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // insure any columns with offcol set are effectively set to last value

  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j]) array_total[i][j] = norm*array[i][j];

  // output result to file

  if (fp && me == 0) {
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d\n",ntimestep,nrows);
    for (i = 0; i < nrows; i++) {
      fprintf(fp,"%d",i+1);
      for (j = 0; j < nvalues; j++) fprintf(fp,format,array_total[i][j]/norm);
      fprintf(fp,"\n");
    }
    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      if (fileend > 0) ftruncate(fileno(fp),fileend);
    }
  }
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

int FixAveTime::column_length(int dynamic)
{
  int m,length,lengthone;

  // determine nrows for static values

  if (!dynamic) {
    length = 0;
    for (int i = 0; i < nvalues; i++) {
      if (varlen[i]) continue;
      if (which[i] == COMPUTE) {
        int icompute = modify->find_compute(ids[i]);
        if (argindex[i] == 0)
          lengthone = modify->compute[icompute]->size_vector;
        else lengthone = modify->compute[icompute]->size_array_rows;
      } else if (which[i] == FIX) {
        int ifix = modify->find_fix(ids[i]);
        if (argindex[i] == 0) lengthone = modify->fix[ifix]->size_vector;
        else lengthone = modify->fix[ifix]->size_array_rows;
      } else if (which[i] == VARIABLE) {
        // variables are always varlen = 1, so dynamic
      }
      if (length == 0) length = lengthone;
      else if (lengthone != length)
        error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
    }
  }

  // determine new nrows for dynamic values
  // either all must be the same
  // or must match other static values
  // don't need to check if not MODE = VECTOR, just invoke lock_length()

  if (dynamic) {
    length = 0;
    for (int i = 0; i < nvalues; i++) {
      if (varlen[i] == 0) continue;
      m = value2index[i];
      if (which[i] == COMPUTE) {
        Compute *compute = modify->compute[m];
        lengthone = compute->lock_length();
      } else if (which[i] == VARIABLE) {
        double *varvec;
        lengthone = input->variable->compute_vector(m,&varvec);
      }
      if (mode == SCALAR) continue;
      if (all_variable_length) {
        if (length == 0) length = lengthone;
        else if (lengthone != length)
          error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
      } else {
        if (lengthone != nrows)
          error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
      }
    }
  }

  return length;
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveTime::compute_scalar()
{
  if (norm) return vector_total[0]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveTime::compute_vector(int i)
{
  if (i >= nrows) return 0.0;
  if (norm) {
    if (mode == SCALAR) return vector_total[i]/norm;
    if (mode == VECTOR) return array_total[i][0]/norm;
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveTime::compute_array(int i, int j)
{
  if (i >= nrows) return 0.0;
  if (norm) return array_total[i][j]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveTime::options(int iarg, int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  noff = 0;
  offlist = NULL;
  overwrite = 0;
  format_user = NULL;
  format = (char *) " %g";
  title1 = NULL;
  title2 = NULL;
  title3 = NULL;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          snprintf(str,128,"Cannot open fix ave/time file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/time command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/time command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/time command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Illegal fix ave/time command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"off") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      memory->grow(offlist,noff+1,"ave/time:offlist");
      offlist[noff++] = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      delete [] format_user;
      int n = strlen(arg[iarg+1]) + 2;
      format_user = new char[n];
      sprintf(format_user," %s",arg[iarg+1]);
      format = format_user;
      iarg += 2;
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
    } else error->all(FLERR,"Illegal fix ave/time command");
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays for mode = VECTOR of size Nrows x Nvalues
------------------------------------------------------------------------- */

void FixAveTime::allocate_arrays()
{
  memory->destroy(array);
  memory->destroy(array_total);
  memory->create(array,nrows,nvalues,"ave/time:array");
  memory->create(array_total,nrows,nvalues,"ave/time:array_total");
  if (ave == WINDOW) {
    memory->destroy(array_list);
    memory->create(array_list,nwindow,nrows,nvalues,"ave/time:array_list");
  }

  // reinitialize regrown array_total since it accumulates

  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nvalues; j++) array_total[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveTime::nextvalid()
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

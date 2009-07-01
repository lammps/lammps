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

#include "stdlib.h"
#include "string.h"
#include "fix_ave_time.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING,WINDOW};
enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix ave/time command");

  time_depend = 1;

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  nfreq = atoi(arg[5]);

  // parse values until one isn't recognized

  which = new int[narg-6];
  argindex = new int[narg-6];
  ids = new char*[narg-6];
  value2index = new int[narg-6];
  nvalues = 0;

  int iarg = 6;
  while (iarg < narg) {
    if ((strncmp(arg[iarg],"c_",2) == 0) || 
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
	  error->all("Illegal fix ave/time command");
	argindex[nvalues] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else break;

    iarg++;
  }

  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/time command");
      if (me == 0) {
	fp = fopen(arg[iarg+1],"w");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open fix ave/time file %s",arg[iarg+1]);
	  error->one(str);
	}
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all("Illegal fix ave/time command");
      if (ave == WINDOW) {
	if (iarg+3 > narg) error->all("Illegal fix ave/time command");
	nwindow = atoi(arg[iarg+2]);
	if (nwindow <= 0) error->all("Illegal fix ave/time command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/time command");
      startstep = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal fix ave/time command");
  }

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/time command");
  if (nfreq < nevery || nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all("Illegal fix ave/time command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
	error->all("Fix ave/time compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
	error->all("Fix ave/time compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
	error->all("Fix ave/time compute vector is accessed out-of-range");
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
	error->all("Fix ave/time fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
	error->all("Fix ave/time fix does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
	error->all("Fix ave/time fix vector is accessed out-of-range");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all("Variable name for fix ave/time does not exist");
      if (input->variable->equalstyle(ivariable) == 0)
	error->all("Fix ave/time variable is not equal-style variable");
    }
  }

  // print header into file

  if (fp && me == 0) {
    fprintf(fp,"# Time-averaged data for fix %s\n",id);
    fprintf(fp,"# TimeStep");
    for (int i = 0; i < nvalues; i++)
      if (which[i] == COMPUTE) fprintf(fp," c_%s",ids[i]);
      else if (which[i] == FIX) fprintf(fp," f_%s",ids[i]);
      else if (which[i] == VARIABLE) fprintf(fp," v_%s",ids[i]);
    fprintf(fp,"\n");
  }

  // allocate and initialize memory for averaging

  vector = new double[nvalues];
  vector_total = new double[nvalues];
  for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;

  vector_list = NULL;
  if (ave == WINDOW)
    vector_list = memory->create_2d_double_array(nwindow,nvalues,
						 "ave/time:vector_list");

  // this fix produces either a global scalar or vector
  // intensive/extensive flags set by compute,fix,variable that produces value

  scalar_vector_freq = nfreq;
  extlist = NULL;

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
    } else if (which[0] == VARIABLE)
      extscalar = 0;

  } else {
    vector_flag = 1;
    size_vector = nvalues;
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
      } else if (which[i] == VARIABLE)
	extlist[i] = 0;
    }
  }

  // initializations
  // set vector total to zero since it accumulates

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;
  for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // can be this timestep if multiple of nfreq and nrepeat = 1
  // else backup from next multiple of nfreq

  nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;

  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  delete [] which;
  delete [] argindex;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] value2index;
  delete [] extlist;

  if (fp && me == 0) fclose(fp);
  delete [] vector;
  delete [] vector_total;
  memory->destroy_2d_double_array(vector_list);
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
  // set indices and check validity of all computes,fixes,variables
  // check that fix frequency is acceptable

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for fix ave/time does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0) 
	error->all("Fix ID for fix ave/time does not exist");
      value2index[i] = ifix;

      if (nevery % modify->fix[ifix]->scalar_vector_freq)
	error->all("Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0) 
	error->all("Variable name for fix ave/time does not exist");
      value2index[i] = ivariable;
    }
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveTime::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveTime::end_of_step()
{
  int i,m;

  // skip if not step which requires doing something

  int ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero if first step

  if (irepeat == 0)
    for (i = 0; i < nvalues; i++) vector[i] = 0.0;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (i = 0; i < nvalues; i++) {
    m = value2index[i];

    // invoke compute if not previously invoked

    if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[i] == 0) {
	if (!(compute->invoked_flag & INVOKED_SCALAR)) {
	  compute->compute_scalar();
	  compute->invoked_flag |= INVOKED_SCALAR;
	}
	vector[i] += compute->scalar;
      } else {
	if (!(compute->invoked_flag & INVOKED_VECTOR)) {
	  compute->compute_vector();
	  compute->invoked_flag |= INVOKED_VECTOR;
	}
	vector[i] += compute->vector[argindex[i]-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == FIX) {
      if (argindex[i] == 0) 
	vector[i] += modify->fix[m]->compute_scalar();
      else
	vector[i] += modify->fix[m]->compute_vector(argindex[i]-1);

    // evaluate equal-style variable

    } else if (which[i] == VARIABLE)
      vector[i] += input->variable->compute_equal(m);
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

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (i = 0; i < nvalues; i++) vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ntimestep >= startstep) {
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
  }

  // output result to file

  if (fp && me == 0) {
    fprintf(fp,"%d",ntimestep);
    for (i = 0; i < nvalues; i++) fprintf(fp," %g",vector_total[i]/norm);
    fprintf(fp,"\n");
    fflush(fp);
  }
}

/* ----------------------------------------------------------------------
   return scalar value
   could be ONE, RUNNING, or WINDOW value
------------------------------------------------------------------------- */

double FixAveTime::compute_scalar()
{
  if (norm) return vector_total[0]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Nth vector value
   could be ONE, RUNNING, or WINDOW value
------------------------------------------------------------------------- */

double FixAveTime::compute_vector(int n)
{
  if (norm) return vector_total[n]/norm;
  return 0.0;
}

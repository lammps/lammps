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
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{COMPUTE,FIX};
enum{SCALAR,VECTOR,BOTH};
enum{ONE,RUNNING,WINDOW};

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all("Illegal fix ave/time command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  nfreq = atoi(arg[5]);

  if (strcmp(arg[6],"compute") == 0) which = COMPUTE;
  else if (strcmp(arg[6],"fix") == 0) which = FIX;
  else error->all("Illegal fix ave/time command");

  int n = strlen(arg[7]) + 1;
  id = new char[n];
  strcpy(id,arg[7]);

  // option defaults

  sflag = 1;
  vflag = 0;
  fp = NULL;
  ave = ONE;

  // optional args

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"scalar") == 0) {
	sflag = 1;
	vflag = 0;
      } else if (strcmp(arg[iarg+1],"vector") == 0) {
	sflag = 0;
	vflag = 1;
      } else if (strcmp(arg[iarg+1],"both") == 0) sflag = vflag = 1;
      else error->all("Illegal fix ave/time command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
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
    } else error->all("Illegal fix ave/time command");
  }

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/time command");
  if (nfreq < nevery || nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all("Illegal fix ave/time command");

  int icompute,ifix;
  if (which == COMPUTE) {
    icompute = modify->find_compute(id);
    if (icompute < 0) error->all("Compute ID for fix ave/time does not exist");
  } else {
    ifix = modify->find_fix(id);
    if (ifix < 0) error->all("Fix ID for fix ave/time does not exist");
  }

  if (which == COMPUTE) {
    if (sflag && modify->compute[icompute]->scalar_flag == 0)
      error->all("Fix ave/time compute does not calculate a scalar");
    if (vflag && modify->compute[icompute]->vector_flag == 0)
      error->all("Fix ave/time compute does not calculate a vector");
  } else {
    if (sflag && modify->fix[ifix]->scalar_flag == 0)
      error->all("Fix ave/time fix does not calculate a scalar");
    if (vflag && modify->fix[ifix]->vector_flag == 0)
      error->all("Fix ave/time fix does not calculate a vector");
  }

  if (which == COMPUTE &&
      modify->compute[icompute]->pressflag) pressure_every = nevery;

  // setup list of computes to call, including pre-computes

  compute = NULL;
  if (which == COMPUTE) {
    ncompute = 1 + modify->compute[icompute]->npre;
    compute = new Compute*[ncompute];
  } else ncompute = 0;

  // print header into file

  if (fp && me == 0) {
    if (which == COMPUTE)
      fprintf(fp,"Time-averaged data for fix %s, group %s, and compute %s\n",
	      id,group->names[modify->compute[icompute]->igroup],id);
    else
      fprintf(fp,"Time-averaged data for fix %s, group %s, and fix %s\n",
	      id,group->names[modify->fix[ifix]->igroup],id);
    if (sflag and !vflag)
      fprintf(fp,"TimeStep Value\n");
    else if (!sflag and vflag)
      fprintf(fp,"TimeStep Vector-values\n");
    else if (!sflag and vflag)
      fprintf(fp,"TimeStep Scalar-value Vector-values\n");
  }

  // allocate memory for averaging

  vector = vector_total = NULL;
  if (vflag) {
    if (which == COMPUTE) size_vector = modify->compute[icompute]->size_vector;
    else size_vector = modify->fix[ifix]->size_vector;
    vector = new double[size_vector];
    vector_total = new double[size_vector];
    for (int i = 0; i < size_vector; i++) vector_total[i] = 0.0;
  }

  scalar_list = NULL;
  vector_list = NULL;
  if (sflag && ave == WINDOW) scalar_list = new double[nwindow];
  if (vflag && ave == WINDOW)
    vector_list = memory->create_2d_double_array(nwindow,size_vector,
						 "fix ave/time:vector_list");

  // enable this fix to produce a global scalar and/or vector

  if (sflag) scalar_flag = 1;
  if (vflag) vector_flag = 1;
  scalar_vector_freq = nfreq;
  if (which == COMPUTE) extensive = modify->compute[icompute]->extensive;
  else extensive = modify->fix[ifix]->extensive;

  // initializations

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  // nvalid = next step on which end_of_step does something
  // can be this timestep if multiple of nfreq and nrepeat = 1
  // else backup from next multiple of nfreq

  nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;

  if (nvalid < update->ntimestep)
    error->all("Fix ave/time cannot be started on this timestep");
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  delete [] id;
  if (fp && me == 0) fclose(fp);
  delete [] compute;
  delete [] vector;
  delete [] vector_total;
  delete [] scalar_list;
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
  // set ptrs to one or more computes called each end-of-step

  if (which == COMPUTE) {
    int icompute = modify->find_compute(id);
    if (icompute < 0)
      error->all("Compute ID for fix ave/time does not exist");
  
    ncompute = 0;
    if (modify->compute[icompute]->npre)
      for (int i = 0; i < modify->compute[icompute]->npre; i++) {
	int ic = modify->find_compute(modify->compute[icompute]->id_pre[i]);
	if (ic < 0)
	  error->all("Precompute ID for fix ave/time does not exist");
	compute[ncompute++] = modify->compute[ic];
      }
    
    compute[ncompute++] = modify->compute[icompute];
  }

  // set ptr to fix ID
  // check that fix frequency is acceptable

  if (which == FIX) {
    int ifix = modify->find_fix(id);
    if (ifix < 0) 
      error->all("Fix ID for fix ave/time does not exist");
    fix = modify->fix[ifix];
    if (nevery % fix->scalar_vector_freq)
      error->all("Fix ave/time and fix not computed at compatible times");
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveTime::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveTime::end_of_step()
{
  int i;

  // skip if not step which requires doing something

  if (update->ntimestep != nvalid) return;

  // zero if first step

  if (irepeat == 0) {
    scalar = 0.0;
    if (vflag) 
      for (i = 0; i < size_vector; i++) vector[i] = 0.0;
  }

  // accumulate results of compute or fix to local copy
  
  if (which == COMPUTE) {
    if (sflag) {
      double value;
      for (i = 0; i < ncompute; i++) value = compute[i]->compute_scalar();
      scalar += value;
    }
    if (vflag) {
      for (i = 0; i < ncompute; i++) compute[i]->compute_vector();
      double *cvector = compute[ncompute-1]->vector;
      for (i = 0; i < size_vector; i++) vector[i] += cvector[i];
    }

  } else {
    if (sflag) scalar += fix->compute_scalar();
    if (vflag)
      for (i = 0; i < size_vector; i++)
	vector[i] += fix->compute_vector(i);
  }

  // done if irepeat < nrepeat

  irepeat++;
  nvalid += nevery;
  if (irepeat < nrepeat) return;

  // average the final result for the Nfreq timestep
  // reset irepeat and nvalid

  double repeat = nrepeat;
  if (sflag) scalar /= repeat;
  if (vflag) for (i = 0; i < size_vector; i++) vector[i] /= repeat;

  irepeat = 0;
  nvalid = update->ntimestep+nfreq - (nrepeat-1)*nevery;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, comine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    if (sflag) scalar_total = scalar;
    if (vflag) for (i = 0; i < size_vector; i++) vector_total[i] = vector[i];
    norm = 1;

  } else if (ave == RUNNING) {
    if (sflag) scalar_total += scalar;
    if (vflag) for (i = 0; i < size_vector; i++) vector_total[i] += vector[i];
    norm++;

  } else if (ave == WINDOW) {
    if (sflag) {
      scalar_total += scalar;
      if (window_limit) scalar_total -= scalar_list[iwindow];
      scalar_list[iwindow] = scalar;
    }
    if (vflag) {
      for (i = 0; i < size_vector; i++) {
	vector_total[i] += vector[i];
	if (window_limit) vector_total[i] -= vector_list[iwindow][i];
	vector_list[iwindow][i] = vector[i];
      }
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // output result to file

  if (fp && me == 0) {
    fprintf(fp,"%d",update->ntimestep);
    if (sflag) fprintf(fp," %g",scalar_total/norm);
    if (vflag)
      for (i = 0; i < size_vector; i++) fprintf(fp," %g",vector_total[i]/norm);
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
  if (norm) return scalar_total/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Nth vector value
   could be ONE, RUNNING, or WINDOW value
------------------------------------------------------------------------- */

double FixAveTime::compute_vector(int n)
{
  if (norm) return vector_total[n]/norm;
  else return 0.0;
}

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
enum{SCALAR,VECTOR};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix ave/time command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  nfreq = atoi(arg[5]);

  global_freq = nfreq;
  time_depend = 1;

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

  options(narg,arg);

  // parse values until one isn't recognized

  which = new int[nvalues];
  argindex = new int[nvalues];
  ids = new char*[nvalues];
  value2index = new int[nvalues];
  nvalues = 0;

  iarg = 6;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 || 
	strncmp(arg[iarg],"f_",2) == 0 || 
	strncmp(arg[iarg],"v_",2) == 0) {
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

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0) error->all("Illegal fix ave/time command");
  if (nfreq < nevery || nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all("Illegal fix ave/time command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
	error->all("Fix ave/time compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
	error->all("Fix ave/time compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
	error->all("Fix ave/time compute vector is accessed out-of-range");

    } else if (which[i] == COMPUTE && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
	error->all("Fix ave/time compute does not calculate a vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
	error->all("Fix ave/time compute does not calculate a array");
      if (argindex[i] && 
	  argindex[i] > modify->compute[icompute]->size_array_cols)
	error->all("Fix ave/time compute array is accessed out-of-range");

    } else if (which[i] == FIX && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
	error->all("Fix ave/time fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
	error->all("Fix ave/time fix does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
	error->all("Fix ave/time fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
	error->all("Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == FIX && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
	error->all("Fix ave/time fix does not calculate a vector");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
	error->all("Fix ave/time fix does not calculate a array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_array_cols)
	error->all("Fix ave/time fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
	error->all("Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all("Variable name for fix ave/time does not exist");
      if (input->variable->equalstyle(ivariable) == 0)
	error->all("Fix ave/time variable is not equal-style variable");
      if (mode == VECTOR)
	error->all("Fix ave/time cannot use variable with vector mode");
    }
  }

  // if VECTOR mode, check that all columns are same length
  // nrows = # of rows in output array

  if (mode == VECTOR) {
    int length;

    for (int i = 0; i < nvalues; i++) {
      if (which[i] == COMPUTE) {
	int icompute = modify->find_compute(ids[i]);
	if (argindex[i] == 0) length = modify->compute[icompute]->size_vector;
	else length = modify->compute[icompute]->size_array_rows;
      } else if (which[i] == FIX) {
	int ifix = modify->find_fix(ids[i]);
	if (argindex[i] == 0) length = modify->fix[ifix]->size_vector;
	else length = modify->fix[ifix]->size_array_rows;
      }
      if (i == 0) nrows = length;
      else if (length != nrows)
	error->all("Fix ave/time columns are inconsistent lengths");
    }

    column = new double[nrows];
  } else column = NULL;

  // print file comment lines

  if (fp && me == 0) {
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Time-averaged data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else if (mode == SCALAR) {
      fprintf(fp,"# TimeStep");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",arg[6+i]);
      fprintf(fp,"\n");
    } else fprintf(fp,"# TimeStep Number-of-rows\n");
    if (title3 && mode == VECTOR) fprintf(fp,"%s\n",title3);
    else if (mode == VECTOR) {
      fprintf(fp,"# Row");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",arg[6+i]);
      fprintf(fp,"\n");
    }
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;

  // allocate and initialize memory for averaging

  vector = vector_total = NULL;
  vector_list = NULL;
  array = array_total = NULL;
  array_list = NULL;
  
  if (mode == SCALAR) {
    vector = new double[nvalues];
    vector_total = new double[nvalues];
    if (ave == WINDOW)
      vector_list = memory->create_2d_double_array(nwindow,nvalues,
						   "ave/time:vector_list");
  } else {
    array = memory->create_2d_double_array(nrows,nvalues,"ave/time:array");
    array_total = memory->create_2d_double_array(nrows,nvalues,
						 "ave/time:array_total");
    if (ave == WINDOW)
      array_list = memory->create_3d_double_array(nwindow,nrows,nvalues,
						  "ave/time:array_list");
  }

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

  } else {
    if (nvalues == 1) {
      vector_flag = 1;
      size_vector = nrows;
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
      }

    } else {
      array_flag = 1;
      size_local_rows = nrows;
      size_local_cols = nvalues;
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
	}
	if (value == -1) 
	  error->all("Fix ave/time cannot set output array "
		     "intensive/extensive from these inputs");
	if (i == 0) extarray = value;
	else if (value != extarray)
	  error->all("Fix ave/time cannot set output array "
		     "intensive/extensive from these inputs");
      }
    }
  }

  // initializations
  // set vector_total/array_total to zero since it accumulates

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  if (mode == SCALAR)
    for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;
  else
    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < nvalues; j++) array_total[i][j] = 0.0;

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
  delete [] offcol;

  delete [] which;
  delete [] argindex;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] value2index;
  delete [] extlist;

  if (fp && me == 0) fclose(fp);

  delete [] vector;
  delete [] vector_total;
  delete [] column;
  memory->destroy_2d_double_array(array);
  memory->destroy_2d_double_array(array_total);
  memory->destroy_3d_double_array(array_list);
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
	error->all("Compute ID for fix ave/time does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0) 
	error->all("Fix ID for fix ave/time does not exist");
      value2index[i] = ifix;

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
  // skip if not step which requires doing something

  int ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  if (mode == SCALAR) invoke_scalar(ntimestep);
  else invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_scalar(int ntimestep)
{
  int i,m;
  double scalar;

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
	scalar = compute->scalar;
      } else {
	if (!(compute->invoked_flag & INVOKED_VECTOR)) {
	  compute->compute_vector();
	  compute->invoked_flag |= INVOKED_VECTOR;
	}
	scalar = compute->vector[argindex[i]-1];
      }
      
      // access fix fields, guaranteed to be ready
      
    } else if (which[i] == FIX) {
      if (argindex[i] == 0) 
	scalar = modify->fix[m]->compute_scalar();
      else
	scalar = modify->fix[m]->compute_vector(argindex[i]-1);
      
      // evaluate equal-style variable
      
    } else if (which[i] == VARIABLE)
      scalar = input->variable->compute_equal(m);
    
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

  // insure any columns with offcol set are effectively set to last value

  for (i = 0; i < nvalues; i++)
    if (offcol[i]) vector_total[i] = norm*vector[i];

  // output result to file

  if (fp && me == 0) {
    fprintf(fp,"%d",ntimestep);
    for (i = 0; i < nvalues; i++) fprintf(fp," %g",vector_total[i]/norm);
    fprintf(fp,"\n");
    fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_vector(int ntimestep)
{
  int i,j,m;
  
  // zero if first step

  if (irepeat == 0)
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array[i][j] = 0.0;
  
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
      else
	for (i = 0; i < nrows; i++)
	  column[i] = fix->compute_array(i,argindex[j]);
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
  
  // average the final result for the Nfreq timestep
  
  double repeat = nrepeat;
  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j] == 0) array[i][j] /= repeat;
  
  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values
  
  if (ntimestep >= startstep) {
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
  }
  
  // insure any columns with offcol set are effectively set to last value
  
  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j]) array_total[i][j] = norm*array[i][j];
  
  // output result to file
      
  if (fp && me == 0) {
    fprintf(fp,"%d %d\n",ntimestep,nrows);
    for (i = 0; i < nrows; i++) {
      fprintf(fp,"%d",i+1);
      for (j = 0; j < nvalues; j++) fprintf(fp," %g",array_total[i][j]/norm);
      fprintf(fp,"\n");
    }
    fflush(fp);
  }
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
  if (norm) {
    if (mode == SCALAR) return vector_total[i]/norm;
    if (mode == VECTOR) return array_total[i][0];
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveTime::compute_array(int i, int j)
{
  if (norm) return array_total[i][j]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveTime::options(int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  offcol = new int[nvalues];
  for (int i = 0; i < nvalues; i++) offcol[i] = 0;
  title1 = NULL;
  title2 = NULL;
  title3 = NULL;

  // optional args

  int iarg = 6 + nvalues;
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
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all("Illegal fix ave/time command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"off") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/time command");
      int ncolumn = atoi(arg[iarg+1]);
      if (ncolumn <= 0 || ncolumn > nvalues)
	error->all("Invalid fix ave/time off column");
      offcol[ncolumn-1] = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal fix ave/time command");
  }
}

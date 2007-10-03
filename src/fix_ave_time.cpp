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
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 9) error->all("Illegal fix ave/time command");

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  nfreq = atoi(arg[5]);

  int n = strlen(arg[6]) + 1;
  id_compute = new char[n];
  strcpy(id_compute,arg[6]);

  int flag = atoi(arg[7]);

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[8],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix ave/time file %s",arg[8]);
      error->one(str);
    }
  }

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/time command");
  if (nfreq < nevery || nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all("Illegal fix ave/time command");

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all("Compute ID for fix ave/time does not exist");

  if (flag < 0 || flag > 2) error->all("Illegal fix ave/time command");
  sflag = vflag = 0;
  if (flag == 0 || flag == 2) sflag = 1;
  if (flag == 1 || flag == 2) vflag = 1;

  if (sflag && modify->compute[icompute]->scalar_flag == 0)
    error->all("Fix ave/time compute does not calculate a scalar");
  if (vflag && modify->compute[icompute]->vector_flag == 0)
    error->all("Fix ave/time compute does not calculate a vector");

  if (modify->compute[icompute]->pressflag) pressure_every = nevery;

  // setup list of computes to call, including pre-computes

  ncompute = 1 + modify->compute[icompute]->npre;
  compute = new Compute*[ncompute];

  // print header into file

  if (me == 0) {
    fprintf(fp,"Time-averaged data for fix %s, group %s, and compute %s\n",
	    id,group->names[modify->compute[icompute]->igroup],id_compute);
    if (sflag and !vflag)
      fprintf(fp,"TimeStep Value\n");
    else if (!sflag and vflag)
      fprintf(fp,"TimeStep Vector-values\n");
    else if (!sflag and vflag)
      fprintf(fp,"TimeStep Scalar-value Vector-values\n");
  }

  vector = NULL;
  if (vflag) {
    size_vector = modify->compute[icompute]->size_vector;
    vector = new double[size_vector];
  }

  // nvalid = next step on which end_of_step does something

  irepeat = 0;
  nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  nvalid -= (nrepeat-1)*nevery;
  if (nvalid <= update->ntimestep)
    error->all("Fix ave/time cannot be started on this timestep");
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  delete [] id_compute;
  if (me == 0) fclose(fp);
  delete [] compute;
  delete [] vector;
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

  int icompute = modify->find_compute(id_compute);
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

  // accumulate results of compute to local copy
  
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

  irepeat++;
  nvalid += nevery;

  // output the results
  // reset irepeat and nvalid

  if (irepeat == nrepeat) {
    double repeat = nrepeat;

    if (me == 0) {
      fprintf(fp,"%d",update->ntimestep);
      if (sflag) fprintf(fp," %g",scalar/repeat);
      if (vflag)
	for (i = 0; i < size_vector; i++) fprintf(fp," %g",vector[i]/repeat);
      fprintf(fp,"\n");
      fflush(fp);
    }

    irepeat = 0;
    nvalid = update->ntimestep+nfreq - (nrepeat-1)*nevery;
  }
}

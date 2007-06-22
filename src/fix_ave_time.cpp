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
  if (narg != 8) error->all("Illegal fix ave/time command");

  nevery = atoi(arg[3]);
  nfreq = atoi(arg[4]);

  int n = strlen(arg[5]) + 1;
  id_compute = new char[n];
  strcpy(id_compute,arg[5]);

  int flag = atoi(arg[6]);

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[7],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix ave/time file %s",arg[7]);
      error->one(str);
    }
  }

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/time command");
  if (nfreq < nevery || nfreq % nevery)
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

  nsum = 0;
  scalar = 0.0;
  vector = NULL;
  if (vflag) {
    size_vector = modify->compute[icompute]->size_vector;
    vector = new double[size_vector];
    for (int i = 0; i < size_vector; i++) vector[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  delete [] id_compute;
  if (me == 0) fclose(fp);
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
  // set ptrs to current compute and precompute

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all("Compute ID for fix ave/time does not exist");
  compute = modify->compute[icompute];

  if (compute->id_pre) {
    icompute = modify->find_compute(compute->id_pre);
    if (icompute < 0)
      error->all("Precompute ID for fix ave/time does not exist");
    precompute = modify->compute[icompute];
  } else precompute = NULL;
}

/* ---------------------------------------------------------------------- */

void FixAveTime::end_of_step()
{
  int i;

  if (precompute) {
    if (sflag) double tmp = precompute->compute_scalar();
    if (vflag) precompute->compute_vector();
  }

  nsum++;
  if (sflag) scalar += compute->compute_scalar();
  if (vflag) {
    compute->compute_vector();
    double *cvector = compute->vector;
    for (i = 0; i < size_vector; i++) vector[i] += cvector[i];
  }

  if (update->ntimestep % nfreq == 0) {
    if (me == 0) {
      fprintf(fp,"%d",update->ntimestep);
      if (sflag) fprintf(fp," %g",scalar/nsum);
      if (vflag)
	for (i = 0; i < size_vector; i++) fprintf(fp," %g",vector[i]/nsum);
      fprintf(fp,"\n");
      fflush(fp);
    }

    nsum = 0;
    scalar = 0.0;
    if (vflag) 
      for (i = 0; i < size_vector; i++) vector[i] = 0.0;
 }
}

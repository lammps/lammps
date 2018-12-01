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

#include <mpi.h>
#include <cstring>
#include <cctype>
#include "compute_pair.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{EPAIR,EVDWL,ECOUL};

/* ---------------------------------------------------------------------- */

ComputePair::ComputePair(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  pstyle(NULL), pair(NULL), one(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute pair command");

  scalar_flag = 1;
  extscalar = 1;
  peflag = 1;
  timeflag = 1;

  int n = strlen(arg[3]) + 1;
  if (lmp->suffix) n += strlen(lmp->suffix) + 1;
  pstyle = new char[n];
  strcpy(pstyle,arg[3]);

  int iarg = 4;
  nsub = 0;
  evalue = EPAIR;

  if (narg > iarg) {
    if (isdigit(arg[iarg][0])) {
      nsub = force->inumeric(FLERR,arg[iarg]);
      ++iarg;
      if (nsub <= 0)
        error->all(FLERR,"Illegal compute pair command");
    }
  }

  if  (narg > iarg) {
    if (strcmp(arg[iarg],"epair") == 0) evalue = EPAIR;
    else if (strcmp(arg[iarg],"evdwl") == 0) evalue = EVDWL;
    else if (strcmp(arg[iarg],"ecoul") == 0) evalue = ECOUL;
    else error->all(FLERR, "Illegal compute pair command");
    ++iarg;
  }

  // check if pair style with and without suffix exists

  pair = force->pair_match(pstyle,1,nsub);
  if (!pair && lmp->suffix) {
    strcat(pstyle,"/");
    strcat(pstyle,lmp->suffix);
    pair = force->pair_match(pstyle,1,nsub);
  }

  if (!pair)
    error->all(FLERR,"Unrecognized pair style in compute pair command");
  npair = pair->nextra;

  if (npair) {
    vector_flag = 1;
    size_vector = npair;
    extvector = 1;
    one = new double[npair];
    vector = new double[npair];
  } else one = vector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePair::~ComputePair()
{
  delete [] pstyle;
  delete [] one;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePair::init()
{
  // recheck for pair style in case it has been deleted

  pair = force->pair_match(pstyle,1,nsub);
  if (!pair)
    error->all(FLERR,"Unrecognized pair style in compute pair command");
}

/* ---------------------------------------------------------------------- */

double ComputePair::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  double eng;
  if (evalue == EPAIR) eng = pair->eng_vdwl + pair->eng_coul;
  else if (evalue == EVDWL) eng = pair->eng_vdwl;
  else if (evalue == ECOUL) eng = pair->eng_coul;

  MPI_Allreduce(&eng,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputePair::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_vector)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  for (int i = 0; i < npair; i++)
    one[i] = pair->pvector[i];
  MPI_Allreduce(one,vector,npair,MPI_DOUBLE,MPI_SUM,world);
}

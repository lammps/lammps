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

#include "mpi.h"
#include "string.h"
#include "compute_pair.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePair::ComputePair(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute pair command");
  if (igroup) error->all("Compute pair must use group all");

  pair = force->pair_match(arg[3],1);
  if (!pair) error->all("Unrecognized pair style in compute pair command");
  npair = pair->nextra;
  if (!npair) 
    error->all("Pair style in compute pair command stores no values");

  // settings

  vector_flag = 1;
  size_vector = npair;
  extvector = 1;
  peflag = 1;
  timeflag = 1;

  one = new double[npair];
  vector = new double[npair];
}

/* ---------------------------------------------------------------------- */

ComputePair::~ComputePair()
{
  delete [] one;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePair::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_vector)
    error->all("Energy was not tallied on needed timestep");

  for (int i = 0; i < npair; i++)
    one[i] = pair->pvector[i];
  MPI_Allreduce(one,vector,npair,MPI_DOUBLE,MPI_SUM,world);
}

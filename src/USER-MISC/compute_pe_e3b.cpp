/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   contributing author: Steven E Strong
   stevene.strong at gmail dot com
------------------------------------------------------------------------- */

#include "compute_pe_e3b.h"
#include "pair_e3b.h"

#include "pair.h"
#include "update.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePEE3B::ComputePEE3B(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), e3b(NULL)
{
  //        0  1   2
  //compute ID grp pe/e3b
  if (narg != 3) error->all(FLERR,"Illegal compute pe/e3b command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4; //etotA,etotB,etotC,etot2
  extvector = extscalar = 1;
  timeflag = 1;

  peflag = 1;                   // we need Pair::ev_tally() to be run

  invoked_vector = invoked_scalar = -1;
  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputePEE3B::~ComputePEE3B()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePEE3B::init() {
  Pair *pair = force->pair_match("e3b",false,0);
  if (pair==NULL)
    error->all(FLERR,"This compute must be used with pair_style e3b");

  e3b = (PairE3B *) pair;
  if (e3b==NULL)
    error->all(FLERR,"something went wrong");
}

/* ---------------------------------------------------------------------- */

void ComputePEE3B::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  // sum energies across procs
  MPI_Allreduce(e3b->etot,vector,4,MPI_DOUBLE,MPI_SUM,world);
}

double ComputePEE3B::compute_scalar() {
  invoked_scalar = update->ntimestep;
  if (invoked_scalar != invoked_vector)
    compute_vector();

  scalar = vector[0]+vector[1]+vector[2]+vector[3];
  return scalar;
}

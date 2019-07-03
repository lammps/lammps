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

#include "compute_improper.h"
#include <mpi.h>
#include "update.h"
#include "force.h"
#include "improper.h"
#include "improper_hybrid.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeImproper::ComputeImproper(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  emine(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute improper command");

  vector_flag = 1;
  extvector = 1;
  peflag = 1;
  timeflag = 1;

  // check if improper style hybrid exists

  improper = (ImproperHybrid *) force->improper_match("hybrid");
  if (!improper)
    error->all(FLERR,
               "Improper style for compute improper command is not hybrid");
  size_vector = nsub = improper->nstyles;

  emine = new double[nsub];
  vector = new double[nsub];
}

/* ---------------------------------------------------------------------- */

ComputeImproper::~ComputeImproper()
{
  delete [] emine;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeImproper::init()
{
  // recheck improper style in case it has been changed

  improper = (ImproperHybrid *) force->improper_match("hybrid");
  if (!improper)
    error->all(FLERR,
               "Improper style for compute improper command is not hybrid");
  if (improper->nstyles != nsub)
    error->all(FLERR,"Improper style for compute improper command has changed");
}

/* ---------------------------------------------------------------------- */

void ComputeImproper::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_vector)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  for (int i = 0; i < nsub; i++)
    emine[i] = improper->styles[i]->energy;

  MPI_Allreduce(emine,vector,nsub,MPI_DOUBLE,MPI_SUM,world);
}

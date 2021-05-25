// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_bond.h"

#include "bond.h"
#include "bond_hybrid.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBond::ComputeBond(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  emine(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal compute bond command");

  vector_flag = 1;
  extvector = 1;
  peflag = 1;
  timeflag = 1;

  // check if bond style hybrid exists

  bond = (BondHybrid *) force->bond_match("hybrid");
  if (!bond)
    error->all(FLERR,"Bond style for compute bond command is not hybrid");
  size_vector = nsub = bond->nstyles;

  emine = new double[nsub];
  vector = new double[nsub];
}

/* ---------------------------------------------------------------------- */

ComputeBond::~ComputeBond()
{
  delete [] emine;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeBond::init()
{
  // recheck bond style in case it has been changed

  bond = (BondHybrid *) force->bond_match("hybrid");
  if (!bond)
    error->all(FLERR,"Bond style for compute bond command is not hybrid");
  if (bond->nstyles != nsub)
    error->all(FLERR,"Bond style for compute bond command has changed");
}

/* ---------------------------------------------------------------------- */

void ComputeBond::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_vector)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  for (int i = 0; i < nsub; i++)
    emine[i] = bond->styles[i]->energy;

  MPI_Allreduce(emine,vector,nsub,MPI_DOUBLE,MPI_SUM,world);
}

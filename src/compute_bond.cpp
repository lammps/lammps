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
#include <string.h>
#include "compute_bond.h"
#include "update.h"
#include "force.h"
#include "bond.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBond::ComputeBond(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR,"Illegal compute bond command");
  if (igroup) error->all(FLERR,"Compute bond must use group all");

  scalar_flag = 1;
  extscalar = 1;
  peflag = 1;
  timeflag = 1;

  int n = strlen(arg[3]) + 1;
  if (lmp->suffix) n += strlen(lmp->suffix) + 1;
  bstyle = new char[n];
  strcpy(bstyle,arg[3]);

  // check if bond style with and without suffix exists

  bond = force->bond_match(bstyle);
  if (!bond && lmp->suffix) {
    strcat(bstyle,"/");
    strcat(bstyle,lmp->suffix);
    bond = force->bond_match(bstyle);
  }
  if (!bond)
    error->all(FLERR,"Unrecognized bond style in compute bond command");

  vector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBond::~ComputeBond()
{
  delete [] bstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeBond::init()
{
  // recheck for bond style in case it has been deleted

  bond = force->bond_match(bstyle);

  if (!bond)
    error->all(FLERR,"Unrecognized bond style in compute bond command");
}

/* ---------------------------------------------------------------------- */

double ComputeBond::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  double eng;
  eng = bond->energy;

  MPI_Allreduce(&eng,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  return scalar;
}


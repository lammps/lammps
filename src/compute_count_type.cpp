/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_count_type.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{ATOM,BOND};

/* ---------------------------------------------------------------------- */

ComputeCountType::ComputeCountType(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR, "Illegal compute count/type command");

  // process args

  if (strcmp(arg[3],"atom") == 0) mode = ATOM;
  else if (strcmp(arg[3],"bond") == 0) mode = BOND;
  else error->all(FLERR, "Invalid compute count/type keyword {}",arg[3]);

  if (mode == ATOM) {
    vector_flag = 1;
    size_vector = atom->ntypes;
    extvector = 1;
  } else if (mode == BOND) {
    scalar_flag = vector_flag = 1;
    size_vector = atom->nbondtypes;
    extscalar = 1;
    extvector = 1;
  }

  if (mode == BOND && !atom->avec->bonds_allow)
    error->all(FLERR,"Cannot use compute count/type bond command with no bonds allowed");

  vector = new double[size_vector];

  // work vectors

  count = new int[size_vector];
  bcount_me = new bigint[size_vector];
  bcount = new bigint[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeCountType::~ComputeCountType()
{
  delete[] vector;
}

/* ----------------------------------------------------------------------
   only invoked for mode = BOND to count broken bonds
   broken bonds have bond_type = 0
---------------------------------------------------------------------- */

double ComputeCountType::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int nlocal = atom->nlocal;

  int nbond;
  int count = 0;

  // NOTE: respect group setting
  
  for (int i = 0; i < nlocal; i++) {
    nbond = num_bond[i];
    for (int m = 0; m < nbond; m++)
      if (bond_type[i][m] == 0) count++;
  }

  // sum across procs as bigint, then convert to double
  // correct for double counting if newton_bond off
  
  bigint bcount_me = count;
  bigint bcount;
  MPI_Allreduce(&bcount_me, &bcount, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (force->newton_bond == 0) bcount /= 2;

  if (bcount > MAXDOUBLEINT)
    error->all(FLERR,"Compute count/type overflow");
  double scalar = bcount;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeCountType::compute_vector()
{
  invoked_vector = update->ntimestep;

  int n;

  // count atoms by type

  // NOTE: respect group setting

  if (mode == ATOM) {
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int ntypes = atom->ntypes;

    for (int m = 0; m < ntypes; m++) count[m] = 0;
    for (int i = 0; i < nlocal; i++) count[type[i]-1]++;

    n = ntypes;
  }
  
  // count bonds by type
  // skip type = 0 bonds, they are counted by compute_scalar
  // bond types can be negative for SHAKE

  // NOTE: respect group setting

  else if (mode == BOND) {
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;
    int nlocal = atom->nlocal;
    int nbondtypes = atom->nbondtypes;
    
    int nbond,itype;
    for (int m = 0; m < nbondtypes; m++) count[m] = 0;
    
    for (int i = 0; i < nlocal; i++) {
      nbond = num_bond[i];
      for (int m = 0; m < nbond; m++) {
        itype = bond_type[i][m];
        if (itype == 0) continue;
        if (itype > 0) count[itype-1]++;
        else count[-itype-1]++;
      }
    }

    n = nbondtypes;
  }

  // sum across procs as bigint, then convert to double
  // correct for double counting if newton_bond off

  for (int m = 0; m < n; m++) bcount_me[m] = count[m];
  MPI_Allreduce(&bcount_me, &bcount, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (force->newton_bond == 0) 
    for (int m = 0; m < n; m++) bcount[m] /= 2;

  for (int m = 0; m < n; m++)
    if (bcount[m] > MAXDOUBLEINT)
      error->all(FLERR,"Compute count/type overflow");
  for (int m = 0; m < n; m++) vector[m] *= bcount[m];
}

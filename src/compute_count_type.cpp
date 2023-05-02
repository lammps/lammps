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

enum { ATOM, BOND };

/* ---------------------------------------------------------------------- */

ComputeCountType::ComputeCountType(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR, "Illegal compute count/type command");

  // process args

  if (strcmp(arg[3], "atom") == 0)
    mode = ATOM;
  else if (strcmp(arg[3], "bond") == 0)
    mode = BOND;
  else
    error->all(FLERR, "Invalid compute count/type keyword {}", arg[3]);

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
    error->all(FLERR, "Cannot use compute count/type bond command with no bonds allowed");

  // output vector

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

  delete[] count;
  delete[] bcount_me;
  delete[] bcount;
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

  int m, nbond;
  int count = 0;

  // count broken bonds with bond_type = 0
  // ignore group setting since 2 atoms in a broken bond
  //   can be arbitrarily far apart

  for (int i = 0; i < nlocal; i++) {
    nbond = num_bond[i];
    for (m = 0; m < nbond; m++)
      if (bond_type[i][m] == 0) count++;
  }

  // sum across procs as bigint, then convert to double
  // correct for double counting if newton_bond off

  bigint bcount_me = count;
  bigint bcount;
  MPI_Allreduce(&bcount_me, &bcount, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (force->newton_bond == 0) bcount /= 2;

  if (bcount > MAXDOUBLEINT) error->all(FLERR, "Compute count/type overflow");
  scalar = bcount;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeCountType::compute_vector()
{
  invoked_vector = update->ntimestep;

  int nvec;

  // count atoms by type
  // atom must be in group to be counted

  if (mode == ATOM) {
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int ntypes = atom->ntypes;

    for (int m = 0; m < ntypes; m++) count[m] = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) count[type[i] - 1]++;

    nvec = ntypes;
  }

  // count bonds by type
  // both atoms in bond must be in group to be counted
  // skip type = 0 bonds, they are counted by compute_scalar
  // bond types can be negative for SHAKE

  else if (mode == BOND) {
    tagint **bond_atom = atom->bond_atom;
    int **bond_type = atom->bond_type;
    int *num_bond = atom->num_bond;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int nbondtypes = atom->nbondtypes;

    int j, m, nbond, itype;
    int flag = 0;
    for (int m = 0; m < nbondtypes; m++) count[m] = 0;

    for (int i = 0; i < nlocal; i++) {
      nbond = num_bond[i];
      for (m = 0; m < nbond; m++) {
        itype = bond_type[i][m];
        if (itype == 0) continue;

        j = atom->map(bond_atom[i][m]);
        if (j < 0) {
          flag = 1;
          continue;
        }

        if ((mask[i] & groupbit) && (mask[j] & groupbit)) {
          if (itype > 0)
            count[itype - 1]++;
          else
            count[-itype - 1]++;
        }
      }
    }

    int flagany;
    MPI_Allreduce(&flag, &flagany, 1, MPI_INT, MPI_SUM, world);
    if (flagany) error->all(FLERR, "Missing bond atom in compute count/type");

    nvec = nbondtypes;
  }

  // sum across procs as bigint, then convert to double
  // correct for double counting if newton_bond off

  for (int m = 0; m < nvec; m++) bcount_me[m] = count[m];
  MPI_Allreduce(bcount_me, bcount, nvec, MPI_LMP_BIGINT, MPI_SUM, world);
  if (force->newton_bond == 0)
    for (int m = 0; m < nvec; m++) bcount[m] /= 2;

  for (int m = 0; m < nvec; m++)
    if (bcount[m] > MAXDOUBLEINT) error->all(FLERR, "Compute count/type overflow");
  for (int m = 0; m < nvec; m++) vector[m] = bcount[m];
}

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
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "update.h"

using namespace LAMMPS_NS;

enum { ATOM, BOND, ANGLE, DIHEDRAL, IMPROPER };

/* ---------------------------------------------------------------------- */

ComputeCountType::ComputeCountType(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), count(nullptr), bcount_me(nullptr), bcount(nullptr)
{
  if (narg != 4) error->all(FLERR, "Incorrect number of args for compute count/type command");

  // process args

  if (strcmp(arg[3], "atom") == 0)
    mode = ATOM;
  else if (strcmp(arg[3], "bond") == 0)
    mode = BOND;
  else if (strcmp(arg[3], "angle") == 0)
    mode = ANGLE;
  else if (strcmp(arg[3], "dihedral") == 0)
    mode = DIHEDRAL;
  else if (strcmp(arg[3], "improper") == 0)
    mode = IMPROPER;
  else
    error->all(FLERR, "Invalid compute count/type keyword {}", arg[3]);

  // error check

  if (mode == BOND && !atom->nbondtypes)
    error->all(FLERR, "Compute count/type bond command with no bonds defined");
  if (mode == ANGLE && !atom->nangletypes)
    error->all(FLERR, "Compute count/type bond command with no angles defined");
  if (mode == DIHEDRAL && !atom->ndihedraltypes)
    error->all(FLERR, "Compute count/type dihedral command with no dihedrals defined");
  if (mode == IMPROPER && !atom->nimpropertypes)
    error->all(FLERR, "Compute count/type improper command with no impropers defined");

  // set vector lengths

  if (mode == ATOM) {
    vector_flag = 1;
    size_vector = atom->ntypes;
    extvector = 1;
  } else if (mode == BOND) {
    scalar_flag = vector_flag = 1;
    size_vector = atom->nbondtypes;
    extscalar = 1;
    extvector = 1;
  } else if (mode == ANGLE) {
    vector_flag = 1;
    size_vector = atom->nangletypes;
    extvector = 1;
  } else if (mode == DIHEDRAL) {
    vector_flag = 1;
    size_vector = atom->ndihedraltypes;
    extvector = 1;
  } else if (mode == IMPROPER) {
    vector_flag = 1;
    size_vector = atom->nimpropertypes;
    extvector = 1;
  }

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

  // count broken bonds with bond_type = 0
  // ignore group setting since 2 atoms in a broken bond
  //   can be arbitrarily far apart

  int count = 0;
  for (int i = 0; i < nlocal; i++) {
    int nbond = num_bond[i];
    for (int m = 0; m < nbond; m++)
      if (bond_type[i][m] == 0) count++;
  }

  // sum across procs as bigint, then convert to double
  // correct for double counting if newton_bond off

  bigint bcount = 0;
  bigint bcount_me = count;
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

  if (mode == ATOM)
    nvec = count_atoms();
  else if (mode == BOND)
    nvec = count_bonds();
  else if (mode == ANGLE)
    nvec = count_angles();
  else if (mode == DIHEDRAL)
    nvec = count_dihedrals();
  else if (mode == IMPROPER)
    nvec = count_impropers();
  else
    nvec = 0;

  // sum across procs as bigint, then convert to double
  // correct for multiple counting if newton_bond off

  for (int m = 0; m < nvec; m++) bcount_me[m] = count[m];
  MPI_Allreduce(bcount_me, bcount, nvec, MPI_LMP_BIGINT, MPI_SUM, world);

  if (force->newton_bond == 0) {
    if (mode == BOND)
      for (int m = 0; m < nvec; m++) bcount[m] /= 2;
    else if (mode == ANGLE)
      for (int m = 0; m < nvec; m++) bcount[m] /= 3;
    if (mode == DIHEDRAL || mode == IMPROPER)
      for (int m = 0; m < nvec; m++) bcount[m] /= 4;
  }

  for (int m = 0; m < nvec; m++)
    if (bcount[m] > MAXDOUBLEINT) error->all(FLERR, "Compute count/type overflow");
  for (int m = 0; m < nvec; m++) vector[m] = bcount[m];
}

/* ----------------------------------------------------------------------
   count atoms by type
   atom must be in group to be counted
---------------------------------------------------------------------- */

int ComputeCountType::count_atoms()
{
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  for (int m = 0; m < ntypes; m++) count[m] = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) count[type[i] - 1]++;

  return ntypes;
}

/* ----------------------------------------------------------------------
   count bonds by type
   both atoms in bond must be in group to be counted
   skip type = 0 bonds, they are counted by compute_scalar()
   bond types can be negative, count them as if positive
---------------------------------------------------------------------- */

int ComputeCountType::count_bonds()
{
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nbondtypes = atom->nbondtypes;

  int flag = 0;
  for (int m = 0; m < nbondtypes; m++) count[m] = 0;

  for (int i = 0; i < nlocal; i++) {
    int nbond = num_bond[i];
    for (int m = 0; m < nbond; m++) {
      int itype = bond_type[i][m];
      if (itype == 0) continue;

      int j = atom->map(bond_atom[i][m]);
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

  return nbondtypes;
}

/* ----------------------------------------------------------------------
   count angles by type
   all 3 atoms in angle must be in group to be counted
   angle types can be negative, count them as if positive
---------------------------------------------------------------------- */

int ComputeCountType::count_angles()
{
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int *num_angle = atom->num_angle;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nangletypes = atom->nangletypes;

  int flag = 0;
  for (int m = 0; m < nangletypes; m++) count[m] = 0;

  for (int i = 0; i < nlocal; i++) {
    int nangle = num_angle[i];
    for (int m = 0; m < nangle; m++) {
      int itype = angle_type[i][m];

      int j1 = atom->map(angle_atom1[i][m]);
      int j2 = atom->map(angle_atom2[i][m]);
      int j3 = atom->map(angle_atom3[i][m]);
      if (j1 < 0 || j2 < 0 || j3 < 0) {
        flag = 1;
        continue;
      }

      if ((mask[j1] & groupbit) && (mask[j2] & groupbit) && (mask[j3] & groupbit)) {
        if (itype > 0)
          count[itype - 1]++;
        else if (itype < 0)
          count[-itype - 1]++;
      }
    }
  }

  int flagany;
  MPI_Allreduce(&flag, &flagany, 1, MPI_INT, MPI_SUM, world);
  if (flagany) error->all(FLERR, "Missing angle atom in compute count/type");

  return nangletypes;
}

/* ----------------------------------------------------------------------
   count dihedrals by type
   all 4 atoms in dihedral must be in group to be counted
   dihedral types can be negative, count them as if positive
---------------------------------------------------------------------- */

int ComputeCountType::count_dihedrals()
{
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int *num_dihedral = atom->num_dihedral;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ndihedraltypes = atom->ndihedraltypes;

  int flag = 0;
  for (int m = 0; m < ndihedraltypes; m++) count[m] = 0;

  for (int i = 0; i < nlocal; i++) {
    int ndihedral = num_dihedral[i];
    for (int m = 0; m < ndihedral; m++) {
      int itype = dihedral_type[i][m];

      int j1 = atom->map(dihedral_atom1[i][m]);
      int j2 = atom->map(dihedral_atom2[i][m]);
      int j3 = atom->map(dihedral_atom3[i][m]);
      int j4 = atom->map(dihedral_atom4[i][m]);
      if (j1 < 0 || j2 < 0 || j3 < 0 || j4 < 0) {
        flag = 1;
        continue;
      }

      if ((mask[j1] & groupbit) && (mask[j2] & groupbit) && (mask[j3] & groupbit) &&
          (mask[j4] & groupbit)) {
        if (itype > 0)
          count[itype - 1]++;
        else if (itype < 0)
          count[-itype - 1]++;
      }
    }
  }

  int flagany;
  MPI_Allreduce(&flag, &flagany, 1, MPI_INT, MPI_SUM, world);
  if (flagany) error->all(FLERR, "Missing dihedral atom in compute count/type");

  return ndihedraltypes;
}

/* ----------------------------------------------------------------------
   count impropers by type
   all 4 atoms in improper must be in group to be counted
   improper types can be negative, count them as if positive
---------------------------------------------------------------------- */

int ComputeCountType::count_impropers()
{
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int *num_improper = atom->num_improper;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nimpropertypes = atom->nimpropertypes;

  int flag = 0;
  for (int m = 0; m < nimpropertypes; m++) count[m] = 0;

  for (int i = 0; i < nlocal; i++) {
    int nimproper = num_improper[i];
    for (int m = 0; m < nimproper; m++) {
      int itype = improper_type[i][m];

      int j1 = atom->map(improper_atom1[i][m]);
      int j2 = atom->map(improper_atom2[i][m]);
      int j3 = atom->map(improper_atom3[i][m]);
      int j4 = atom->map(improper_atom4[i][m]);
      if (j1 < 0 || j2 < 0 || j3 < 0 || j4 < 0) {
        flag = 1;
        continue;
      }

      if ((mask[j1] & groupbit) && (mask[j2] & groupbit) && (mask[j3] & groupbit) &&
          (mask[j4] & groupbit)) {
        if (itype > 0)
          count[itype - 1]++;
        else if (itype < 0)
          count[-itype - 1]++;
      }
    }
  }

  int flagany;
  MPI_Allreduce(&flag, &flagany, 1, MPI_INT, MPI_SUM, world);
  if (flagany) error->all(FLERR, "Missing improper atom in compute count/type");

  return nimpropertypes;
}

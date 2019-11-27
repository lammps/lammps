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

#include "atom_vec_molecular.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMolecular::AtomVecMolecular(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 1;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 1;
  mass_type = 1;

  atom->molecule_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) 
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "nspecial special";
  fields_copy = (char *)
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "nspecial special";
  fields_comm = NULL;
  fields_comm_vel = NULL;
  fields_reverse = NULL;
  fields_border = (char *) "molecule";
  fields_border_vel = (char *) "molecule";
  fields_exchange = (char *)
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "nspecial special";
  fields_restart = (char *) 
    "molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4";
  fields_create = (char *) 
    "molecule num_bond num_angle num_dihedral num_improper nspecial";
  fields_data_atom = (char *) "id molecule type x";
  fields_data_vel = NULL;

  setup_fields();

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  bond_negative = angle_negative = dihedral_negative = improper_negative = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecMolecular::~AtomVecMolecular()
{
  delete [] bond_negative;
  delete [] angle_negative;
  delete [] dihedral_negative;
  delete [] improper_negative;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_restart() to pack
------------------------------------------------------------------------- */

void AtomVecMolecular::pack_restart_pre(int i)
{
  // insure negative vectors are needed length

  if (bond_per_atom < atom->bond_per_atom) {
    delete [] bond_negative;
    bond_per_atom = atom->bond_per_atom;
    bond_negative = new int[bond_per_atom];
  }
  if (angle_per_atom < atom->angle_per_atom) {
    delete [] angle_negative;
    angle_per_atom = atom->angle_per_atom;
    angle_negative = new int[angle_per_atom];
  }
  if (dihedral_per_atom < atom->dihedral_per_atom) {
    delete [] dihedral_negative;
    dihedral_per_atom = atom->dihedral_per_atom;
    dihedral_negative = new int[dihedral_per_atom];
  }
  if (improper_per_atom < atom->improper_per_atom) {
    delete [] improper_negative;
    improper_per_atom = atom->improper_per_atom;
    improper_negative = new int[improper_per_atom];
  }

  // flip any negative types to positive and flag which ones

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_type = atom->dihedral_type;
  int *num_improper = atom->num_improper;
  int **improper_type = atom->improper_type;

  int any_bond_negative = 0;
  for (int m = 0; m < num_bond[i]; m++) {
    if (bond_type[i][m] < 0) {
      bond_negative[m] = 1;
      bond_type[i][m] = -bond_type[i][m];
      any_bond_negative = 1;
    } else bond_negative[m] = 0;
  }

  int any_angle_negative = 0;
  for (int m = 0; m < num_angle[i]; m++) {
    if (angle_type[i][m] < 0) {
      angle_negative[m] = 1;
      angle_type[i][m] = -angle_type[i][m];
      any_angle_negative = 1;
    } else angle_negative[m] = 0;
  }

  int any_dihedral_negative = 0;
  for (int m = 0; m < num_dihedral[i]; m++) {
    if (dihedral_type[i][m] < 0) {
      dihedral_negative[m] = 1;
      dihedral_type[i][m] = -dihedral_type[i][m];
      any_dihedral_negative = 1;
    } else dihedral_negative[m] = 0;
  }

  int any_improper_negative = 0;
  for (int m = 0; m < num_improper[i]; m++) {
    if (improper_type[i][m] < 0) {
      improper_negative[m] = 1;
      improper_type[i][m] = -improper_type[i][m];
      any_improper_negative = 1;
    } else improper_negative[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_restart()
------------------------------------------------------------------------- */

void AtomVecMolecular::pack_restart_post(int i)
{
  // restore the flagged types to their negative values

  if (any_bond_negative) {
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;
    for (int m = 0; m < num_bond[i]; m++)
      if (bond_negative[m]) bond_type[i][m] = -bond_type[i][m];
  }

  if (any_angle_negative) {
    int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;
    for (int m = 0; m < num_angle[i]; m++)
      if (angle_negative[m]) angle_type[i][m] = -angle_type[i][m];
  }

  if (any_dihedral_negative) {
    int *num_dihedral = atom->num_dihedral;
    int **dihedral_type = atom->dihedral_type;
    for (int m = 0; m < num_dihedral[i]; m++)
      if (dihedral_negative[m]) dihedral_type[i][m] = -dihedral_type[i][m];
  }

  if (any_improper_negative) {
    int *num_improper = atom->num_improper;
    int **improper_type = atom->improper_type;
    for (int m = 0; m < num_improper[i]; m++)
      if (improper_negative[m]) improper_type[i][m] = -improper_type[i][m];
  }
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecMolecular::unpack_restart_init(int ilocal)
{
  atom->nspecial[ilocal][0] = 0;
  atom->nspecial[ilocal][1] = 0;
  atom->nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecMolecular::data_atom_post(int ilocal)
{
  atom->num_bond[ilocal] = 0;
  atom->num_angle[ilocal] = 0;
  atom->num_dihedral[ilocal] = 0;
  atom->num_improper[ilocal] = 0;
  atom->nspecial[ilocal][0] = 0;
  atom->nspecial[ilocal][1] = 0;
  atom->nspecial[ilocal][2] = 0;
}

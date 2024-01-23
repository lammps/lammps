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

#include "atom_vec_molecular.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMolecular::AtomVecMolecular(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::MOLECULAR;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 1;
  mass_type = PER_TYPE;

  atom->molecule_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  // clang-format off
  fields_grow = {"molecule", "num_bond", "bond_type", "bond_atom", "num_angle", "angle_type",
    "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type", "dihedral_atom1",
    "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper", "improper_type",
    "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4", "nspecial", "special"};
  fields_copy = {"molecule", "num_bond", "bond_type", "bond_atom", "num_angle", "angle_type",
    "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type", "dihedral_atom1",
    "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper", "improper_type",
    "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4", "nspecial", "special"};
  fields_border = {"molecule"};
  fields_border_vel = {"molecule"};
  fields_exchange = {"molecule", "num_bond", "bond_type", "bond_atom", "num_angle", "angle_type",
    "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type", "dihedral_atom1",
    "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper", "improper_type",
    "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4", "nspecial", "special"};
  fields_restart = {"molecule", "num_bond", "bond_type", "bond_atom", "num_angle", "angle_type",
    "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type", "dihedral_atom1",
    "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper", "improper_type",
    "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4"};
  fields_create = {"molecule", "num_bond", "num_angle", "num_dihedral", "num_improper", "nspecial"};
  fields_data_atom = {"id", "molecule", "type", "x"};
  fields_data_vel = {"id", "v"};
  // clang-format on

  setup_fields();

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  bond_negative = angle_negative = dihedral_negative = improper_negative = nullptr;
}

/* ---------------------------------------------------------------------- */

AtomVecMolecular::~AtomVecMolecular()
{
  delete[] bond_negative;
  delete[] angle_negative;
  delete[] dihedral_negative;
  delete[] improper_negative;
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecMolecular::grow_pointers()
{
  num_bond = atom->num_bond;
  bond_type = atom->bond_type;
  num_angle = atom->num_angle;
  angle_type = atom->angle_type;
  num_dihedral = atom->num_dihedral;
  dihedral_type = atom->dihedral_type;
  num_improper = atom->num_improper;
  improper_type = atom->improper_type;
  nspecial = atom->nspecial;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_restart() to pack
------------------------------------------------------------------------- */

void AtomVecMolecular::pack_restart_pre(int ilocal)
{
  // ensure negative vectors are needed length

  if (bond_per_atom < atom->bond_per_atom) {
    delete[] bond_negative;
    bond_per_atom = atom->bond_per_atom;
    bond_negative = new int[bond_per_atom];
  }
  if (angle_per_atom < atom->angle_per_atom) {
    delete[] angle_negative;
    angle_per_atom = atom->angle_per_atom;
    angle_negative = new int[angle_per_atom];
  }
  if (dihedral_per_atom < atom->dihedral_per_atom) {
    delete[] dihedral_negative;
    dihedral_per_atom = atom->dihedral_per_atom;
    dihedral_negative = new int[dihedral_per_atom];
  }
  if (improper_per_atom < atom->improper_per_atom) {
    delete[] improper_negative;
    improper_per_atom = atom->improper_per_atom;
    improper_negative = new int[improper_per_atom];
  }

  // flip any negative types to positive and flag which ones

  any_bond_negative = 0;
  for (int m = 0; m < num_bond[ilocal]; m++) {
    if (bond_type[ilocal][m] < 0) {
      bond_negative[m] = 1;
      bond_type[ilocal][m] = -bond_type[ilocal][m];
      any_bond_negative = 1;
    } else
      bond_negative[m] = 0;
  }

  any_angle_negative = 0;
  for (int m = 0; m < num_angle[ilocal]; m++) {
    if (angle_type[ilocal][m] < 0) {
      angle_negative[m] = 1;
      angle_type[ilocal][m] = -angle_type[ilocal][m];
      any_angle_negative = 1;
    } else
      angle_negative[m] = 0;
  }

  any_dihedral_negative = 0;
  for (int m = 0; m < num_dihedral[ilocal]; m++) {
    if (dihedral_type[ilocal][m] < 0) {
      dihedral_negative[m] = 1;
      dihedral_type[ilocal][m] = -dihedral_type[ilocal][m];
      any_dihedral_negative = 1;
    } else
      dihedral_negative[m] = 0;
  }

  any_improper_negative = 0;
  for (int m = 0; m < num_improper[ilocal]; m++) {
    if (improper_type[ilocal][m] < 0) {
      improper_negative[m] = 1;
      improper_type[ilocal][m] = -improper_type[ilocal][m];
      any_improper_negative = 1;
    } else
      improper_negative[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_restart()
------------------------------------------------------------------------- */

void AtomVecMolecular::pack_restart_post(int ilocal)
{
  // restore the flagged types to their negative values

  if (any_bond_negative) {
    for (int m = 0; m < num_bond[ilocal]; m++)
      if (bond_negative[m]) bond_type[ilocal][m] = -bond_type[ilocal][m];
  }

  if (any_angle_negative) {
    for (int m = 0; m < num_angle[ilocal]; m++)
      if (angle_negative[m]) angle_type[ilocal][m] = -angle_type[ilocal][m];
  }

  if (any_dihedral_negative) {
    for (int m = 0; m < num_dihedral[ilocal]; m++)
      if (dihedral_negative[m]) dihedral_type[ilocal][m] = -dihedral_type[ilocal][m];
  }

  if (any_improper_negative) {
    for (int m = 0; m < num_improper[ilocal]; m++)
      if (improper_negative[m]) improper_type[ilocal][m] = -improper_type[ilocal][m];
  }
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecMolecular::unpack_restart_init(int ilocal)
{
  nspecial[ilocal][0] = 0;
  nspecial[ilocal][1] = 0;
  nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecMolecular::data_atom_post(int ilocal)
{
  num_bond[ilocal] = 0;
  num_angle[ilocal] = 0;
  num_dihedral[ilocal] = 0;
  num_improper[ilocal] = 0;
  nspecial[ilocal][0] = 0;
  nspecial[ilocal][1] = 0;
  nspecial[ilocal][2] = 0;
}
